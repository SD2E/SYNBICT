#!/usr/bin/env python3
"""Check a run of run_test.sh against the expected netlists and the references.

Checks per circuit:

  1. netlist  -- gates / wires / inputs / outputs vs <dataset>/expected/netlists/
  2. topology -- no broken edges: every input promoter is produced by some gate or is a
                 primary input, every produced promoter is consumed, no isolated gate,
                 the reporter is driven
  3. truth table (Cello hex circuits only) -- the hex recomputed from the Yosys table
                 must equal the circuit name, which is the same string as OUTPUT_OR in
                 cello/reference/cello_logic/<C>_A000_*.txt

Input convention for check 3 (validated on 43 of the 44 published hex circuits): Yosys
drives the three sensor promoters directly and every input is ACTIVE-LOW, so the hex bit
index is j = 4*(1-pBAD) + 2*(1-pTet) + 1*(1-pTac); bit 0 = no inducer added.

Circuits whose name is not a truth table (named Cello circuits, all MD5 plasmids) are
checked on netlist + topology only -- see TESTING.md for what MD5 validation still needs.

Usage: python check_results.py <output_dir> [--dataset cello|md5]
"""
import json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
SENSORS = {'pTac', 'pTet', 'pBAD'}
HEX_NAME = re.compile(r'^0x[0-9A-Fa-f]{2}$')


def hex_from_log(path):
    """Recompute a Cello circuit's hex identity from the Yosys truth table in a run log."""
    if not os.path.exists(path):
        return None, 'no log'
    txt = open(path, errors='ignore').read()
    if 'Traceback' in txt:
        return None, 'crashed'
    hdr = re.search(r'^(\S+(?: \S+)*)\s*\|\s*(\S+)\s*$', txt, re.M)
    rows = re.findall(r'^([01](?: [01])*)\s*\|\s*([01x])\s*$', txt, re.M)
    if not hdr or not rows:
        return None, 'no truth table'
    if any(o == 'x' for _, o in rows):
        return None, 'undefined (x) - combinational loop'
    cols = hdr.group(1).split()
    if not SENSORS <= set(cols):
        return None, f'inputs are {cols}, not the 3 sensors'
    bits = 0
    for vals, out in rows:
        v = dict(zip(cols, [int(x) for x in vals.split()]))
        if out == '1':
            bits |= 1 << (4 * (1 - v['pBAD']) + 2 * (1 - v['pTet']) + 1 * (1 - v['pTac']))
    return f'0x{bits:02X}', None


def topology_problems(net):
    gates, wires = net['gates'], net['wires']
    produced = {g['output'] for g in gates if not g.get('is_output')}
    consumed = {i for g in gates for i in g['inputs']}
    primary = set(net.get('primary_inputs', [])) | SENSORS
    touched = {x for w in wires for x in w}
    p = []
    if not wires and len(gates) > 1:
        p.append('NO WIRES (discrete)')
    p += [f"gate output is not a promoter: {g['id']}:{g['cds']}"
          for g in gates if not g.get('is_output') and g['output'] == g['cds']]
    p += [f'input produced by nobody: {i}' for i in sorted(consumed - produced - primary)]
    p += [f'produced but unconsumed: {o}' for o in sorted(produced - consumed)]
    # An isolated *reporter/output* gate is legitimate: an MD5 quorum-sensing output
    # gene driven straight off a sensor promoter has neither an in- nor an out-edge.
    # An isolated logic gate is not.
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
    print(f"{'circuit':<34}{'netlist':<10}{'topology':<10}{'truth table':<26}result")
    print('-' * 96)
    failed = skipped_tt = 0
    for c in ran:
        got = json.load(open(os.path.join(out_dir, f'{c}_circuit_netlist.json')))
        exp_p = os.path.join(exp_dir, f'{c}_circuit_netlist.json')
        net_ok = os.path.exists(exp_p) and norm(got) == norm(json.load(open(exp_p)))
        net_txt = 'match' if net_ok else ('DIFFERS' if os.path.exists(exp_p) else 'no ref')
        probs = topology_problems(got)
        if HEX_NAME.match(c) or c.endswith('_fixed'):
            tt, why = hex_from_log(os.path.join(out_dir, f'{c}.log'))
            want = c.replace('_fixed', '')
            tt_ok = tt is not None and tt.upper() == want.upper()
            tt_txt = f'{tt} vs {want}' if tt else f'({why})'
        else:
            tt_ok, tt_txt = True, 'n/a (name is not a table)'
            skipped_tt += 1
        verdict = 'PASS' if (net_ok and not probs and tt_ok) else 'FAIL'
        failed += verdict == 'FAIL'
        print(f'{c:<34}{net_txt:<10}{"clean" if not probs else "BROKEN":<10}{tt_txt:<26}{verdict}')
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
