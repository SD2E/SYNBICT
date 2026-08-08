#!/usr/bin/env python3
"""Check a run of run_test.sh against the expected netlists and the references.

Checks per circuit:

  1. netlist  -- gates / wires / inputs / outputs vs <dataset>/expected/netlists/
  2. topology -- no broken edges. A finding is reported when a gate's input promoter is
                 neither produced by another gate nor a primary input; when a gate produces
                 a promoter no gate consumes; when a gate appears in no wire at all; or when
                 a logic gate's output could not be resolved to a promoter. OUTPUT/reporter
                 gates are exempt from the isolated-gate rule (an MD5 quorum-sensing output
                 gene, or an unused branch of a multi-output Cello circuit, legitimately has
                 no edges), and there is no separate "reporter is driven" rule. The MD5
                 dataset also switches off the unconsumed-promoter and no-wires rules: its
                 plasmids are PARTITIONS of one circuit spread across cells and wired by
                 quorum-sensing signals, so an unwired fragment is the expected shape.
  3. truth table -- our Yosys table compared ROW BY ROW against Cello's own table in
                 cello/reference/cello_logic/<C>_A000_*.txt, matching rows by the input
                 promoters' states.

Nothing about the input convention is hard-coded. The reference file states it per circuit:
its INPUT rows give each sensor promoter's waveform and its OUTPUT rows give the reporter's,
column by column, so a column IS an input state and the mapping falls out of the file. This
matters because the column order is NOT constant -- across the 56 reference files the three
sensors appear in 6 different orders, and only 12 of 32 hex circuits have an OUTPUT_OR string
that reads as the circuit's own name. Matching by state instead of by position sidesteps all
of it, and works for the 2-input circuits and the multi-output ones (whose reporter row is
selected by the reporter gene our netlist actually annotates).

Circuits with no reference file, and every MD5 plasmid, are checked on netlist + topology
only -- see TESTING.md for what MD5 validation still needs.

Usage: python check_results.py <output_dir> [--dataset cello|md5]
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SENSORS = {'pTac', 'pTet', 'pBAD'}   # only for the topology rules, not for the truth table

_REF_INPUT = re.compile(r'^INPUT\s+([01]+)\s+input_(\S+)', re.M)
_REF_OUTPUT = re.compile(r'^OUTPUT\S*\s+([01]+)\s+output_(\S+)', re.M)


def our_tables(path):
    """Parse the Yosys truth table out of a run log, one table per output.

    A multi-output circuit prints several output columns (`pTac pTet | PhlF AmeR YFP`),
    so this returns {output name -> {state -> bit}}, where a state is a tuple of
    (promoter, value) pairs sorted by promoter name. Returns (None, reason) on failure;
    an output whose column contains an `x` is dropped (a combinational loop leaves it
    undefined) and reported only if every output is.
    """
    if not os.path.exists(path):
        return None, 'no log'
    txt = open(path, errors='ignore').read()
    if 'Traceback' in txt:
        return None, 'crashed'
    hdr = re.search(r'^(\S+(?: \S+)*)\s*\|\s*(\S+(?: \S+)*)\s*$', txt, re.M)
    rows = re.findall(r'^([01](?: [01])*)\s*\|\s*([01x](?: [01x])*)\s*$', txt, re.M)
    if not hdr or not rows:
        return None, 'no truth table'
    ins, outs = hdr.group(1).split(), hdr.group(2).split()
    tables = {name: {} for name in outs}
    undefined = set()
    for vals, res in rows:
        state = tuple(sorted(zip(ins, (int(x) for x in vals.split()))))
        bits = res.split()
        if len(bits) != len(outs):
            return None, 'truth table row does not match its header'
        for name, b in zip(outs, bits):
            if b == 'x':
                undefined.add(name)
            else:
                tables[name][state] = int(b)
    tables = {k: v for k, v in tables.items() if k not in undefined}
    if not tables:
        return None, 'undefined (x) - combinational loop'
    return tables, None


def reference_table(circuit, reporter):
    """Cello's published table for a circuit, as {state -> output bit}.

    `reporter` is the reporter gene our netlist annotates (YFP, ...); a multi-output
    circuit publishes one OUTPUT row per reporter and only that one is comparable.
    """
    path = os.path.join(HERE, 'cello', 'reference', 'cello_logic',
                        f'{circuit}_A000_logic_circuit.txt')
    if not os.path.exists(path):
        return None, 'no reference file'
    txt = open(path, errors='ignore').read()
    inputs = [(name, wave) for wave, name in _REF_INPUT.findall(txt)]
    outputs = {name.upper(): wave for wave, name in _REF_OUTPUT.findall(txt)}
    if not inputs or not outputs:
        return None, 'reference file has no INPUT/OUTPUT rows'
    out_wave = outputs.get((reporter or '').upper())
    if out_wave is None:
        return None, f'reference has no output_{reporter} row (has {sorted(outputs)})'
    width = len(out_wave)
    if any(len(w) != width for _, w in inputs):
        return None, 'reference INPUT/OUTPUT rows differ in width'
    table = {}
    for col in range(width):
        state = tuple(sorted((name, int(wave[col])) for name, wave in inputs))
        table[state] = int(out_wave[col])
    return table, None


def compare_truth_table(log_path, circuit, netlist):
    """(verdict, detail) from comparing our table to Cello's, row by row."""
    tables, why = our_tables(log_path)
    if tables is None:
        # no table at all is a failure only if the reference has one to compare against
        ref, ref_why = reference_table(circuit, None)
        return (None, ref_why) if ref_why == 'no reference file' else (False, why)
    # A multi-output circuit yields several outputs, only some of which are reporters
    # the reference publishes a row for -- compare on one it knows.
    ref = None
    for reporter in tables:
        ref, why = reference_table(circuit, reporter)
        if ref is not None:
            ours = tables[reporter]
            break
    if ref is None:
        return None, why
    ref_ins = {n for st in ref for n, _ in st}
    our_ins = {n for st in ours for n, _ in st}
    if ref_ins != our_ins:
        return False, f'inputs {sorted(our_ins)} vs reference {sorted(ref_ins)}'
    wrong = [st for st, v in ref.items() if ours.get(st) != v]
    if wrong:
        st = sorted(wrong)[0]
        shown = ', '.join(f'{n}={v}' for n, v in st)
        return False, (f'{len(ref) - len(wrong)}/{len(ref)} rows, first wrong at {shown}: '
                       f'ours={ours.get(st)} ref={ref[st]}')
    return True, f'{len(ref)}/{len(ref)} rows'


def topology_problems(net, partitioned=False):
    gates, wires = net['gates'], net['wires']
    produced = {g['output'] for g in gates if not g.get('is_output')}
    consumed = {i for g in gates for i in g['inputs']}
    primary = set(net.get('primary_inputs', [])) | SENSORS
    touched = {x for w in wires for x in w}
    p = []
    if not wires and len(gates) > 1 and not partitioned:
        p.append('NO WIRES (discrete)')
    p += [f"gate output is not a promoter: {g['id']}:{g['cds']}"
          for g in gates if not g.get('is_output') and g['output'] == g['cds']]
    p += [f'input produced by nobody: {i}' for i in sorted(consumed - produced - primary)]
    if not partitioned:
        p += [f'produced but unconsumed: {o}' for o in sorted(produced - consumed)]
    # An isolated *reporter/output* gate is legitimate: an MD5 quorum-sensing output
    # gene driven straight off a sensor promoter has neither an in- nor an out-edge.
    # An isolated logic gate is not.
    if not partitioned:
        p += [f"isolated gate: {g['id']}:{g['cds']}"
              for g in gates
              if g['id'] not in touched and len(gates) > 1 and not g.get('is_output')]
    return p


def norm(net):
    """Netlist reduced to what must be reproducible (gate ids are positional)."""
    return {
        'gates': sorted((g['type'], g['cds'], tuple(sorted(g['inputs'])), g['output'])
                        for g in net['gates']),
        'primary_inputs': sorted(net['primary_inputs']),
        'primary_outputs': sorted(net['primary_outputs']),
        'n_wires': len(net['wires']),
    }


def main(out_dir, dataset):
    exp_dir = os.path.join(HERE, dataset, 'expected', 'netlists')
    ran = sorted(f.replace('_circuit_netlist.json', '')
                 for f in os.listdir(out_dir) if f.endswith('_circuit_netlist.json'))
    if not ran:
        print(f'no netlists in {out_dir} -- did run_test.sh produce anything?')
        return 1
    print(f'dataset: {dataset}   circuits run: {len(ran)}\n')
    print(f"{'circuit':<34}{'netlist':<10}{'topology':<10}{'truth table':<34}result")
    print('-' * 96)
    failed = skipped_tt = 0
    for c in ran:
        got = json.load(open(os.path.join(out_dir, f'{c}_circuit_netlist.json')))
        exp_p = os.path.join(exp_dir, f'{c}_circuit_netlist.json')
        net_ok = os.path.exists(exp_p) and norm(got) == norm(json.load(open(exp_p)))
        net_txt = 'match' if net_ok else ('DIFFERS' if os.path.exists(exp_p) else 'no ref')
        probs = topology_problems(got, partitioned=(dataset == 'md5'))
        if dataset == 'cello':
            # a repaired sequence is checked against the circuit it repairs
            verdict, tt_txt = compare_truth_table(os.path.join(out_dir, f'{c}.log'),
                                                  c.replace('_fixed', ''), got)
        else:
            verdict, tt_txt = None, 'n/a (no published table)'
        if verdict is None:          # nothing to compare against
            tt_ok = True
            skipped_tt += 1
            tt_txt = f'n/a ({tt_txt})' if not tt_txt.startswith('n/a') else tt_txt
        else:
            tt_ok = verdict
        verdict = 'PASS' if (net_ok and not probs and tt_ok) else 'FAIL'
        failed += verdict == 'FAIL'
        print(f'{c:<34}{net_txt:<10}{"clean" if not probs else "BROKEN":<10}{tt_txt:<34}{verdict}')
        if net_txt == 'DIFFERS':
            e, g = norm(json.load(open(exp_p))), norm(got)
            print(f'      expected {len(e["gates"])} gates / {e["n_wires"]} wires, '
                  f'got {len(g["gates"])} / {g["n_wires"]}')
            for x in sorted(set(map(str, e['gates'])) - set(map(str, g['gates']))):
                print(f'      missing gate: {x}')
            for x in sorted(set(map(str, g['gates'])) - set(map(str, e['gates']))):
                print(f'      unexpected  : {x}')
        for x in probs:
            print(f'      {x}')
    print('-' * 96)
    print(f'{len(ran) - failed}/{len(ran)} PASS'
          + (f'   ({skipped_tt} checked on netlist+topology only)' if skipped_tt else ''))
    return 1 if failed else 0


if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    ds = 'cello'
    if '--dataset' in sys.argv:
        ds = sys.argv[sys.argv.index('--dataset') + 1]
        args = [a for a in args if a != ds]
    if len(args) != 1:
        sys.exit(__doc__)
    sys.exit(main(args[0], ds))
