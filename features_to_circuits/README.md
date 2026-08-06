# features_to_circuits — from annotated sequence to logic

Four scripts, run in this order. The first is upstream SYNBICT; the other three are the
gate layer added in 2.1.

```
annotated SBOL ──features_to_circuits.py──► circuit SBOL (molecular interaction graph)
                        │ with -gn
                        ▼
                  gate_assembler.py ──► *_circuit_netlist.json  (gates + wires)
                        │                        │
                        │                        ├──circuit_to_truth_table.py──► truth table
                        │                        └──netlist_to_graphml.py──────► *.graphml
```

| Script | Input | Output |
|---|---|---|
| `features_to_circuits.py` | annotated SBOL + sub-circuit library | circuit SBOL (`ModuleDefinition`s of Interactions) |
| `gate_assembler.py` | that circuit SBOL | netlist JSON: `{gates, wires, primary_inputs, primary_outputs, gaps}` |
| `circuit_to_truth_table.py` | netlist JSON | truth table (via Yosys), optionally the Verilog |
| `netlist_to_graphml.py` | netlist JSON | GraphML for Cytoscape |

---

## Why the gate layer exists

`features_to_circuits.py` emits a **molecular interaction graph**: production
(CDS → protein), repression / activation (protein → promoter) and transcription
(promoter → CDS). For `0xEA` that is 9 `ModuleDefinition`s and 14 nodes.

That graph does not say **which parts form a gate, what type it is, or how gates wire
together**, so the circuit's function cannot be read off it. `gate_assembler.py` closes
that gap for Cello-style repressor circuits.

## The gate model

Every Cello gate is one transcriptional unit:

```
[input promoter(s)] [RBS] [repressor CDS] [terminator]
repressor protein  --represses-->  its cognate output promoter
```

so:

| netlist field | comes from |
|---|---|
| `inputs` | the promoters in the **same transcriptional unit** as the CDS (positional) |
| `cds` | the repressor CDS |
| `output` | the promoter repressed by that CDS's protein (CDS → protein → promoter, from Interactions) |
| `type` | `NOR` (≥2 inputs) / `NOT` (1 input) / `OUTPUT` (CDS that represses nothing = reporter) |
| `wires` | gate A's `output` promoter appears in gate B's `inputs` |

Four steps in `GateAssembler.assemble()`:

1. `parse_regulation()` — Interactions → `cds_to_protein` (production) and
   `protein_to_repressed` (inhibition).
2. `_construct_parts()` — SequenceAnnotations → parts ordered by position, typed by SO role,
   collapsing duplicate annotations at the same locus.
3. `_segment()` — **split into transcriptional units at terminators**.
4. `assemble()` — one gate per CDS; wire by promoter-name identity.

Note that the promoter→CDS link is *not* taken from the transcription Interaction; it is
positional co-membership in a TU. This is what makes terminator annotation critical
(see Troubleshooting).

## Usage

```bash
# 1. circuit SBOL + gate netlist   (-gn is what triggers the gate layer)
python features_to_circuits/features_to_circuits.py \
    -n http://examples.org \
    -c example/jet_libs/cello_library.xml \
    -t 0xEA_annotated.xml \
    -o 0xEA_circuit.xml \
    -m 1000 -gn
#    -> 0xEA_circuit.xml  and  0xEA_circuit_netlist.json

# 2. truth table
python features_to_circuits/circuit_to_truth_table.py 0xEA_circuit_netlist.json \
    --yosys $(which yosys) [--verilog 0xEA.v]

# 3. graph for Cytoscape
python features_to_circuits/netlist_to_graphml.py 0xEA_circuit_netlist.json -o 0xEA_gates.graphml
```

Run it by **script path, not `python -m`** — the package has no `__main__`.
(`sequences_to_features` does support `-m`.)

### `-gn` / `--gate_netlist`

| Argument | Short | Type | Description |
|---|---|---|---|
| `--gate_netlist` | `-gn` | Boolean | Also assemble the circuit into a logic-gate netlist, written next to the output file as `<output_base>_netlist.json`. Default off. |

### Netlist format

```json
{
  "construct": "_0xEA_comp",
  "gates": [
    {"id": "g1", "cds": "SrpR", "inputs": ["pBAD", "pTet"], "output": "pSrpR",
     "type": "NOR", "is_output": false},
    {"id": "g3", "cds": "YFP",  "inputs": ["pSrpR", "pAmtR"], "output": "YFP",
     "type": "OUTPUT", "is_output": true}
  ],
  "wires": [["g1", "g3"], ["g2", "g3"]],
  "primary_inputs":  ["pBAD", "pTac", "pTet"],
  "primary_outputs": ["YFP"],
  "gaps": []
}
```

`gaps` lists what the assembler could not resolve (e.g. a promoter-only TU with no
annotated CDS). A non-empty `gaps` means the netlist is incomplete.

## circuit_to_truth_table.py

Netlist → structural Verilog → Yosys `eval` over all 2^n input assignments. Gate
semantics:

```
NOT / NOR gate : output_promoter = ~(OR of input promoters)
OUTPUT gate    : reporter        =  (OR of input promoters)
```

A pure-Python evaluation runs alongside and is cross-checked against Yosys (Yosys is the
authority; the Python pass catches a broken Yosys invocation). Output looks like:

```
pBAD pTac pTet  |  YFP
----------------------
0 0 0  |  0
...
output bitstring: 0x37 (00110111)
```

`x` in the table means Yosys could not resolve the value — almost always a **combinational
loop** caused by a mis-segmented netlist, not a Yosys problem.

To compare against a published Cello circuit name, the sensor inputs must be mapped to the
paper's convention (all inputs active-low, bit index
`j = 4*(1-pBAD) + 2*(1-pTet) + 1*(1-pTac)`). `test_bundle/check_results.py` implements this;
it reproduces the published name for 43 of the 44 hex circuits.

## Assumptions and limits

Documented in the `gate_assembler.py` docstring; the practical ones:

1. **Repressor logic only.** `NOR` / `NOT` / reporter `OUTPUT`. An *activator*-based gate is
   not typed as a gate — its CDS gets no output promoter and is emitted as an `OUTPUT`, so
   its downstream edge is missing. In Cello libraries activation appears only in the input
   sensors (AraC→pBAD, LuxR→pLuxStar), which is exactly how it is treated here.
2. **Tandem promoters OR-combine**, so the gate is `NOR(inputs)`.
3. **One output promoter per repressor** — fan-out to a second promoter is not modelled.
4. **Pure Boolean** — cooperativity/multimerisation are ignored (they change the analog
   transfer function, not the truth table).
5. **Forward strand only.**
6. **Single reporter.** Multi-output circuits (`demultiplexer`, `priority_detector`) get
   only the branch that reaches the annotated reporter; other repressors become dangling
   `OUTPUT` gates.

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Truth table is all `x`, netlist has a loop | a terminator was not annotated, so two TUs merged and the downstream CDS inherited the upstream promoters | check the annotation: BLAST must run `-task blastn` (megablast's word_size 28 misses the 47 bp `L3S3P11`) |
| One giant gate, `wires: []`, all promoters as inputs | same merge, extreme case — usually **`-nms` was used with a library containing composite cassettes** (`engineered_region` parts such as `S3_SrpR` swallow their own terminator) | drop `-nms` from the annotation step |
| Every gate is an isolated `OUTPUT` with `wires: []` | `cds_to_protein` is empty — the production Interaction's template role URI was not recognised | `parse_regulation()` accepts `SBO:0000645`, `SBO:0000010` and `http://sbols.org/v2#template`; check which one your converter emits |
| `SBOL_ERROR_URI_NOT_UNIQUE` | re-running on a file that already contains proteins/interactions | re-annotate from FASTA first |
| A repressor CDS shows up as `OUTPUT` | the promoter it represses is not on this sequence, or the library has no inhibition Interaction for it | expected for multi-output circuits; otherwise check the library |

## Tests

`test_bundle/` runs both features over 6 published Cello circuits and checks netlist,
topology and truth table against the circuit names. See `test_bundle/TESTING.md`.
