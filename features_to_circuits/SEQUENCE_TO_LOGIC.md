# Sequence → Circuit → Logic: design notes

How SYNBICT gets from a DNA sequence to a testable logic function, and where the
three complementary approaches — **deterministic assembly**, **SAT**, and a
**GNN** — each belong.

## Pipeline

```
sequence (FASTA)
  └─ sequences_to_features ──▶ annotated SBOL            (parts on the sequence)
       └─ features_to_circuits ──▶ circuit SBOL          (discrete interaction graph:
            │                                              repression / production / activation)
            └─ gate_assembler ──▶ gate netlist (JSON)     (gates + gate-to-gate wiring)
                 └─ circuit_to_truth_table ──▶ truth table (via Yosys)
```

The first two stages already existed. `gate_assembler.py` and
`circuit_to_truth_table.py` are the new pieces that turn the *discrete* molecular
graph into a *connected* gate netlist and then into a function.

### New modules

| Module | Role |
| --- | --- |
| `gate_assembler.py` | interaction graph + linear annotations → gate netlist. Groups a transcriptional unit `[input promoter(s)][RBS][repressor CDS][terminator]` into one gate; the gate's output is the promoter its product represses. Wires gate A → gate B when A's output promoter drives B's CDS. Reports `gaps` for anything it cannot resolve. |
| `circuit_to_truth_table.py` | netlist → structural Verilog → **Yosys `eval`** over all `2^n` inputs → truth table. A pure-Python evaluation cross-checks Yosys. |
| `features_to_circuits.py` | new `-gn/--gate_netlist` flag: after building the circuit, also emit `<output>_netlist.json`. |

### Gate model (Cello repressor logic)

* **repressor gate** (`NOT`/`NOR`): `output_promoter = NOR(input promoters)` — the
  repressor is expressed when any input promoter is on, turning the output
  promoter off.
* **reporter gate** (`OUTPUT`): a CDS that represses nothing (e.g. YFP); its logic
  is `OR(input promoters)` and it is a primary output.

## Why the graph was "discrete", and how it is closed

`features_to_circuits` emits per-motif interactions but never wires gates, so the
raw graph fragments into disconnected pieces with no readable function. The
assembler closes it using the **linear order** of the annotated parts
(transcriptional units) plus the **repression map** (which CDS represses which
promoter). Gaps appear only when a part is missing from the annotation — e.g. the
output reporter that is not in the parts library.

## The three approaches and when each applies

| Situation | Tool |
| --- | --- |
| Structure complete / unambiguous | **deterministic assembler** — pure parse, no search |
| Structure incomplete but recoverable (missing reporter that is identifiable from sequence, cognate promoter known) | assembler + fill — still deterministic |
| Structure genuinely ambiguous **and** the truth table is known | **SAT/SMT** — constrain "assembled logic == truth table" and solve for the wiring |
| Structure ambiguous **and** no ground truth (unknown parts, real-world) | **GNN** — predict edges / roles / function from sequence + graph context |

### Role of the two external inputs

The **UCF** (parts collection) is the *grammar*: part identities **and** the
regulation map (who represses whom). The **Verilog** is the *answer*: the truth
table.

|  | needs UCF (grammar) | needs Verilog (answer) |
| --- | --- | --- |
| deterministic assembler | **yes** | no |
| SAT | **yes** (to form variables) | **yes** (to constrain) |
| GNN | no | no |

Consequences:

* **No Verilog** → SAT loses its disambiguating constraint. If the structure is
  sufficient, use the assembler; if not, it is a prediction problem → GNN. SAT is
  rarely the answer.
* **No UCF** → there is no vocabulary of parts/regulation at all; SAT cannot even
  build its variables. Inferring "who represses whom" must come from the
  **sequence** (operator/binding-site motifs) → this is squarely GNN /
  sequence-inference territory.

**SAT has the strictest prerequisites** (UCF *and* Verilog) and only fills a
genuinely ambiguous wiring gap. For generalizing to unknown circuits — the
motivating goal — the GNN is the tool; SAT is a narrow, edge-case helper.

## Validation

* Synthetic complete circuit (two `NOT` → `OR` reporter): assembler yields 3 gates
  + 2 wires; Yosys truth table = **NAND (`0x7`)**; Python cross-check agrees.
* Real **0xEA** (`example` Cello design): assembler recovers `g1 = NOR(pBAD,pTet)→SrpR→pSrpR`,
  `g2 = NOT(pTac)→AmtR→pAmtR`, inputs `{pBAD,pTac,pTet}`, and flags the one gap
  (the reporter transcriptional unit `pSrpR,pAmtR → ?` has no annotated CDS). That
  region is the **YFP reporter** (identifiable directly from the sequence:
  `ATGGTGAGCAAGGGCGAGGAG…`). Adding the reporter gate and running Yosys reproduces
  `out = NOT(pTac AND (pBAD OR pTet))`, which **matches `0xEA.v` exactly** (input
  map `A=pTac, B=pBAD, C=pTet`). No SAT was needed.

## Known limitations

* Forward-strand annotations only (reverse-complement parts are skipped in v1).
* Cello repressor logic only (`NOR`/`NOT` + `OR` output); activation-based gates
  are not yet typed.
* A missing reporter is reported as a gap; fix by adding the reporter (YFP/GFP +
  RBS) to the feature library so it gets annotated, or by an ORF heuristic on the
  trailing unannotated region.
* The packed output bitstring uses this tool's sorted input order; to compare to a
  Cello design name (`{C,B,A}`) re-order inputs via the UCF sensor assignment (the
  logic function is invariant to this).

## Running it

```bash
conda activate synbict_conda

# annotate: sequence -> annotated SBOL   (sequences_to_features, see repo README)

# assemble circuit + emit gate netlist
python features_to_circuits/features_to_circuits.py -n http://examples.org \
    -c example/Eco1C1G1T1_collection.xml -t OXEA_annotated.xml -o oxea_circuit.xml -gn

# netlist -> truth table (Yosys lives in the base env)
python features_to_circuits/circuit_to_truth_table.py oxea_circuit_netlist.json \
    --yosys /home/sophia/miniconda3/bin/yosys
```
