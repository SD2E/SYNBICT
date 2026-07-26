#!/usr/bin/env python3
"""Compute a circuit's truth table from a gate netlist (via Yosys).

Input: a netlist dict produced by gate_assembler.GateAssembler:
    {primary_inputs:[...], primary_outputs:[...],
     gates:[{id,type,inputs:[promoter...],output,is_output}], wires:[...]}

Gate semantics (Cello repressor logic):
    repressor gate (NOT/NOR): output_promoter = ~(OR of input promoters)  -> NOR
    reporter gate (OUTPUT)  : reporter        =  (OR of input promoters)  -> OR

Pipeline: netlist -> structural Verilog -> Yosys `eval` over all 2^n input
assignments -> truth table. A pure-Python evaluation is also computed and the
two are cross-checked (Yosys is the authority; Python guards against a broken
Yosys call). Assumes a COMPLETE netlist (no gaps); on an incomplete netlist the
dangling intermediate promoters are simply reported as the outputs.
"""
import argparse
import json
import re
import subprocess
import tempfile
import os
import logging

logger = logging.getLogger('synbict')

_EVAL_RE = re.compile(r"Eval result:\s*\\?(\S+)\s*=\s*\d+'([01xz]+)")


def _san(name):
    """Sanitize a part name into a valid Verilog identifier."""
    s = re.sub(r'\W', '_', name)
    if s and s[0].isdigit():
        s = '_' + s
    return s


def netlist_to_verilog(netlist, top='top'):
    """Emit structural Verilog for the netlist. Returns (verilog, in_names, out_names)."""
    ins = [_san(n) for n in netlist['primary_inputs']]
    outs = [_san(n) for n in netlist['primary_outputs']]

    # every signal that is driven by a gate
    driven = {_san(g['output']) for g in netlist['gates']}
    intermediates = sorted(driven - set(outs))

    lines = [f'module {top}(']
    ports = [f'  input {n}' for n in ins] + [f'  output {n}' for n in outs]
    lines.append(',\n'.join(ports))
    lines.append(');')
    for w in intermediates:
        lines.append(f'  wire {w};')

    for g in netlist['gates']:
        o = _san(g['output'])
        terms = [_san(i) for i in g['inputs']] or ["1'b0"]
        or_expr = ' | '.join(terms)
        if g.get('is_output'):
            lines.append(f'  assign {o} = {or_expr};        // reporter (OR)')
        else:
            lines.append(f'  assign {o} = ~({or_expr});     // {g["type"]} (NOR)')
    lines.append('endmodule')
    return '\n'.join(lines), ins, outs


def _python_truth_table(netlist):
    """Reference evaluation in pure Python (topological, promoter activities)."""
    ins = [_san(n) for n in netlist['primary_inputs']]
    outs = [_san(n) for n in netlist['primary_outputs']]
    gates = [dict(output=_san(g['output']),
                  inputs=[_san(i) for i in g['inputs']],
                  is_output=g.get('is_output', False)) for g in netlist['gates']]

    rows = []
    for combo in range(2 ** len(ins)):
        val = {ins[i]: (combo >> i) & 1 for i in range(len(ins))}
        # iterate to a fixed point (handles arbitrary gate order; DAG converges)
        for _ in range(len(gates) + 1):
            changed = False
            for g in gates:
                terms = [val.get(s, 0) for s in g['inputs']]
                orv = 1 if any(terms) else 0
                out = orv if g['is_output'] else (0 if orv else 1)
                if val.get(g['output']) != out:
                    val[g['output']] = out
                    changed = True
            if not changed:
                break
        rows.append((tuple(val[i] for i in ins), tuple(val.get(o, 0) for o in outs)))
    return ins, outs, rows


def truth_table_yosys(netlist, yosys='yosys'):
    """Compute the truth table using Yosys `eval`. Returns (ins, outs, rows)."""
    verilog, ins, outs = netlist_to_verilog(netlist)
    with tempfile.NamedTemporaryFile('w', suffix='.v', delete=False) as f:
        f.write(verilog)
        vpath = f.name
    try:
        script = [f'read_verilog {vpath}', 'hierarchy -top top', 'proc', 'flatten', 'opt']
        combos = []
        for combo in range(2 ** len(ins)):
            sets = ' '.join(f'-set {ins[i]} {(combo >> i) & 1}' for i in range(len(ins)))
            shows = ' '.join(f'-show {o}' for o in outs)
            script.append(f'eval {sets} {shows}')
            combos.append(combo)
        proc = subprocess.run([yosys, '-p', '; '.join(script)],
                              capture_output=True, text=True, check=True)
    finally:
        os.unlink(vpath)

    # parse: each "Executing EVAL pass" starts a combo; following Eval result lines are outputs
    rows, cur = [], {}
    combo_iter = iter(combos)
    cur_combo = None
    for line in proc.stdout.splitlines():
        if 'Executing EVAL pass' in line:
            if cur_combo is not None:
                rows.append(_row(cur_combo, ins, outs, cur))
            cur_combo = next(combo_iter)
            cur = {}
        m = _EVAL_RE.search(line)
        if m:
            cur[m.group(1)] = m.group(2)
    if cur_combo is not None:
        rows.append(_row(cur_combo, ins, outs, cur))
    return ins, outs, rows


def _row(combo, ins, outs, vals):
    invec = tuple((combo >> i) & 1 for i in range(len(ins)))
    outvec = tuple(int(vals.get(o, 'x')) if vals.get(o, 'x') in '01' else 'x' for o in outs)
    return (invec, outvec)


def truth_table(netlist, yosys='yosys', check=True):
    """Yosys truth table, cross-checked against the Python reference."""
    ins, outs, rows = truth_table_yosys(netlist, yosys)
    if check:
        _, _, ref = _python_truth_table(netlist)
        ref_map = {iv: ov for iv, ov in ref}
        for iv, ov in rows:
            if ref_map.get(iv) != ov:
                logger.warning('Yosys/Python truth-table mismatch at %s: yosys=%s python=%s',
                               iv, ov, ref_map.get(iv))
    return ins, outs, rows


def format_table(ins, outs, rows):
    head = ' '.join(ins) + '  |  ' + ' '.join(outs)
    lines = [head, '-' * len(head)]
    for iv, ov in rows:
        lines.append(' '.join(map(str, iv)) + '  |  ' + ' '.join(map(str, ov)))
    return '\n'.join(lines)


def output_bitstring(ins, outs, rows):
    """Pack the (single-output) truth table into a hex string, MSB = all-ones input.
    Mirrors how Cello names designs (e.g. 0xEA)."""
    if len(outs) != 1:
        return None
    by_in = {iv: ov[0] for iv, ov in rows}
    n = len(ins)
    bits = ''
    for combo in range(2 ** n - 1, -1, -1):
        iv = tuple((combo >> i) & 1 for i in range(n))
        bits += str(by_in.get(iv, 'x'))
    try:
        return '0x%X' % int(bits, 2), bits
    except ValueError:
        return None, bits


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='gate netlist -> truth table (via Yosys)')
    ap.add_argument('netlist', help='netlist JSON from gate_assembler')
    ap.add_argument('--verilog', help='also write the generated structural Verilog here')
    ap.add_argument('--yosys', default='yosys', help='path to the yosys binary')
    args = ap.parse_args()

    with open(args.netlist) as f:
        nl = json.load(f)

    if nl.get('gaps'):
        print('WARNING: netlist has gaps (incomplete input); truth table reflects '
              'only the assembled part:')
        for g in nl['gaps']:
            print(f'  - {g}')
        print()

    if args.verilog:
        v, _, _ = netlist_to_verilog(nl)
        with open(args.verilog, 'w') as f:
            f.write(v)
        print(f'wrote Verilog: {args.verilog}\n')

    ins, outs, rows = truth_table(nl, yosys=args.yosys)
    print(format_table(ins, outs, rows))
    bs = output_bitstring(ins, outs, rows)
    if bs:
        print(f'\noutput bitstring: {bs[0]} ({bs[1]})')
