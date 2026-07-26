#!/usr/bin/env python3
"""Convert a gate netlist (gate_assembler output) to GraphML for Cytoscape.

Produces the CONNECTED gate-level graph: primary inputs -> gates -> primary
output, with each gate carrying its type (NOR/NOT/OUTPUT) and each edge carrying
the promoter/signal that flows along it.

Node attributes : label, kind (input|gate|output), gate_type, cds
Edge attributes : signal (the promoter/wire name), relation

Usage:
    python netlist_to_graphml.py netlist.json -o circuit.graphml
In Cytoscape: File > Import > Network from File... , then style by `kind` /
`gate_type` (node) and `signal` (edge).
"""
import argparse
import json
import xml.sax.saxutils as xml


def netlist_to_graphml(netlist):
    inputs = netlist.get('primary_inputs', [])
    gates = netlist.get('gates', [])
    # promoter -> gate id that drives it (produces that promoter)
    producer = {g['output']: g['id'] for g in gates if not g.get('is_output')}

    nodes = {}   # id -> dict(label, kind, gate_type, cds)
    edges = []   # (source, target, signal, relation)

    def add_node(nid, **attrs):
        if nid not in nodes:
            nodes[nid] = attrs

    for p in inputs:
        add_node(p, label=p, kind='input', gate_type='', cds='')

    for g in gates:
        label = f"{g['id']}: {g.get('cds','')} ({g['type']})"
        add_node(g['id'], label=label, kind='gate', gate_type=g['type'], cds=g.get('cds', ''))

    # wire gate inputs
    for g in gates:
        for inp in g.get('inputs', []):
            if inp in producer and producer[inp] != g['id']:
                edges.append((producer[inp], g['id'], inp, 'wire'))       # gate -> gate
            else:
                add_node(inp, label=inp, kind='input', gate_type='', cds='')
                edges.append((inp, g['id'], inp, 'input'))                # primary input -> gate

    # output nodes: reporter gates emit the circuit output signal
    consumed = {inp for g in gates for inp in g.get('inputs', [])}
    for g in gates:
        out = g['output']
        if g.get('is_output'):
            add_node(out, label=out, kind='output', gate_type='', cds='')
            edges.append((g['id'], out, out, 'output'))
        elif out not in consumed:
            # dangling promoter (incomplete netlist) -- still show it as an output
            add_node(out, label=out, kind='output', gate_type='', cds='')
            edges.append((g['id'], out, out, 'output'))

    return _render(nodes, edges)


def _render(nodes, edges):
    keys = [
        ('label', 'node', 'string'), ('kind', 'node', 'string'),
        ('gate_type', 'node', 'string'), ('cds', 'node', 'string'),
        ('signal', 'edge', 'string'), ('relation', 'edge', 'string'),
    ]
    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">']
    for name, dom, typ in keys:
        out.append(f'  <key id="{name}" for="{dom}" attr.name="{name}" attr.type="{typ}"/>')
    out.append('  <graph edgedefault="directed">')
    for nid, a in nodes.items():
        out.append(f'    <node id="{xml.quoteattr(nid)[1:-1]}">')
        for k in ('label', 'kind', 'gate_type', 'cds'):
            out.append(f'      <data key="{k}">{xml.escape(str(a.get(k, "")))}</data>')
        out.append('    </node>')
    for i, (s, t, sig, rel) in enumerate(edges):
        out.append(f'    <edge id="e{i}" source="{xml.quoteattr(s)[1:-1]}" target="{xml.quoteattr(t)[1:-1]}">')
        out.append(f'      <data key="signal">{xml.escape(str(sig))}</data>')
        out.append(f'      <data key="relation">{xml.escape(str(rel))}</data>')
        out.append('    </edge>')
    out.append('  </graph>')
    out.append('</graphml>')
    return '\n'.join(out)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='gate netlist JSON -> GraphML (Cytoscape)')
    ap.add_argument('netlist', help='netlist JSON from gate_assembler')
    ap.add_argument('-o', '--out', required=True, help='output .graphml path')
    args = ap.parse_args()
    with open(args.netlist) as f:
        nl = json.load(f)
    graphml = netlist_to_graphml(nl)
    with open(args.out, 'w') as f:
        f.write(graphml)
    print(f'wrote {args.out}  ({len(nl.get("gates", []))} gates)')
