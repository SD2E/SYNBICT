# Test instructions — gate netlist & truth table

Two features are under test. Everything they need is in this folder (see `MANIFEST.md`).

| # | Feature | Command | Turns | Into |
|---|---------|---------|-------|------|
| 1 | **Gate netlist** | `features_to_circuits.py … -gn` | annotated SBOL | `*_circuit_netlist.json` (gates + wires) |
| 2 | **Truth table** | `circuit_to_truth_table.py` | that netlist | a truth table, via Yosys |

Annotation (FASTA → annotated SBOL) runs first so the test starts from raw sequence, but it
is not what is being tested.

Two datasets:

| | circuits | why it is a good test |
|---|---|---|
| **cello** | 66 single-plasmid circuits, 3 sensor inputs | 44 of them are named after their own truth table — `0xEA` must compute `0xEA`, so the expected answer needs no interpretation |
| **md5** | 70 circular plasmids from the 2-bit MD5 paper | circular targets, 14 bp parts, quorum-sensing wiring — exercises paths the Cello set does not |

---

## 1. Prerequisites

```bash
# SYNBICT + deps (from the repo root; creates the synbict_conda env)
conda env create -f environment.yml
conda activate synbict_conda
pip install -e .

# Yosys — the truth table is computed by simulating the netlist
conda install -c conda-forge yosys      # or: apt install yosys
```

On `PATH`: `blastn`, `makeblastdb`, `yosys`. Python needs `sbol2`, `pysam`, `pandas`.

```bash
which blastn makeblastdb yosys
python -c "import sbol2, pysam, pandas; print('deps ok')"
```

## 2. Run it

```bash
cd test_bundle
./run_test.sh                 # quick: 6 Cello circuits          ~30 s
./run_test.sh --cello         # all 66 Cello circuits            ~5 min
./run_test.sh --md5           # all 70 MD5 plasmids              ~10 min
./run_test.sh --all           # both

python check_results.py out/cello --dataset cello
python check_results.py out/md5   --dataset md5
```

If your Python or Yosys is not the default:

```bash
SYNBICT_PY=/path/to/env/bin/python YOSYS=/path/to/yosys ./run_test.sh --all
```

Each dataset gets its own annotation parameters (the script handles this):

| | min target length | min feature length | extra |
|---|---|---|---|
| cello | `-m 1000` (shortest circuit is 1889 bp; the 2000 default drops it) | `-M 40` | — |
| md5 | `-m 2000` | `-M 14` (parts go down to 14 bp) | `-cir` (circular), `-pid 90` |

`-pid 90` matters: `PJR1` is a 62 bp part that matches over 58 bp (93.5% coverage-weighted
identity), so the 95% default drops it and that locus gets annotated as the spurious
constitutive `Plambda` instead — two plasmids then assemble a different gate.

## 3. Expected result

```
circuit                           netlist   topology  truth table               result
------------------------------------------------------------------------------------------------
0x01                              match     clean     0x01 vs 0x01              PASS
0x07                              match     clean     0x07 vs 0x07              PASS
...
------------------------------------------------------------------------------------------------
6/6 PASS
```

Three independent checks per circuit:

- **netlist** — gate types, CDS, input promoters, output promoter and wire count match
  `<dataset>/expected/netlists/`. Gate *ids* (`g1`, `g2`, …) are positional and not compared.
- **topology** — no broken edges. Four rules:
  - *input produced by nobody* — a gate's input promoter is neither produced by another gate
    nor a primary input, so the edge feeding it is missing;
  - *produced but unconsumed* — a gate makes an output promoter that no gate takes as input,
    so the edge leaving it is missing;
  - *isolated gate* — a gate that appears in no wire at all (the same defect seen from the
    node side). **OUTPUT/reporter gates are exempt**: an MD5 quorum-sensing output gene driven
    straight off a sensor promoter legitimately has no edges, as does the unused branch of a
    multi-output Cello circuit such as `demultiplexer`;
  - *gate output is not a promoter* — a logic gate whose repressor could not be resolved to a
    repressed promoter, so it degenerates to its own CDS name.

  For the **md5** dataset the middle two are switched off, plus the "no wires at all" rule:
  its plasmids are partitions of one circuit spread across cells, so an unwired fragment is
  the expected shape. Note there is no separate "the reporter is driven" rule — a reporter
  with a dangling input is caught by the first rule, but a reporter wired straight to a sensor
  promoter is accepted.
- **truth table** — our Yosys table is compared **row by row** against Cello's own table in
  `cello/reference/cello_logic/<C>_A000_*.txt`, matching rows by the input promoters' states.
  Nothing about the input convention is hard-coded: the reference file's `INPUT` rows give each
  sensor promoter's waveform and its `OUTPUT` rows give each reporter's, column by column, so a
  column *is* an input state. That matters because the column order is **not** constant — across
  the 56 reference files the three sensors appear in 6 different orders, and only 12 of 32 hex
  circuits have an `OUTPUT_OR` string that reads as the circuit's own name. Matching by state
  sidesteps all of it, and covers the 2-input circuits and the multi-output ones (whose reporter
  row is selected by the reporter gene our netlist annotates). Circuits with no reference file
  are marked `n/a` and checked on the first two only.

Known non-PASS cases, expected:

- **`0xB9`** — 7/8 rows. The published sequence is missing the 63 bp `pHlyIIR` promoter, so the
  AmtR gate reads as `NOT(pBAD)` instead of `NOR(pHlyIIR, pBAD)`. `0xB9_fixed` (repaired
  sequence, also in the bundle) matches 8/8. A source defect, not a pipeline bug.
- **`demultiplexer` (1/4 rows) and `priority_detector` (6/8 rows)** — multi-output circuits.
  Cello drives their four/three reporters from separate branches (its `output_YFP` comes from
  the `F1_AmeR` gate alone), but only the YFP branch is on the published sequence, so the other
  branches' promoters end up in the annotated YFP transcriptional unit and YFP reads as an OR of
  two branches. This is limitation 6 in `features_to_circuits/README.md` — now measured instead
  of assumed. Both were invisible before the checker started comparing against the reference.

Current baseline on this code: **cello 63/66 PASS** (the three above), **md5 63/70 PASS**.
Of the 66 Cello circuits, **56 are compared against a published table** (10 have no reference
file: `ANDv2`/`NANDv2`/`XORv2`/… whose Cello entries ship only a `.v`, plus `majority_alt` and
`multiplexer_alt`). The previous checker verified only the 32 hex circuits that had one. The seven
MD5 differences are all the same ambiguity — `PLuxB` and `PLux_u42_TA` are near-identical
quorum-sensing sensor promoters that both match the same locus, and since neither is a
*repressed* promoter the regulation-aware collapse in `gate_assembler._construct_parts()`
cannot tell them apart, so it keeps whichever sorts first. The recorded run picked the other
one. Affected: `seq_002`, `seq_012`, `seq_051`, `seq_055`, `seq_060`, `seq_066`, `seq_067`.
This changes the promoter *name* on a primary input, not the logic — but it does matter when
merging plasmids into a cell, because signals are matched across plasmids by promoter name.

Because the tie is unresolved, **which** plasmids land on which side is not stable: making the
annotated SBOL self-contained (so the aligner paths embed the library definitions they
reference) shuffled the set from six to seven without changing anything about the logic —
`seq_018` and `seq_056` started matching, `seq_012`, `seq_055` and `seq_067` stopped. Any
change that perturbs annotation order can move this boundary. Giving the collapse a
deterministic tie-break is the open issue.

## 4. Reading a result by hand

`out/cello/0xEA.log` ends with the gate table and the truth table:

```
gates:
  g1   NOR    inputs=['pBAD', 'pTet'] cds=SrpR -> output=pSrpR
  g2   NOT    inputs=['pTac'] cds=AmtR -> output=pAmtR
  g3   OUTPUT inputs=['pSrpR', 'pAmtR'] cds=YFP -> output=YFP
wires: [['g1', 'g3'], ['g2', 'g3']]
```

Compare with `cello/reference/cello_logic/0xEA_A000_logic_circuit.txt`:

```
OUTPUT_OR   11101010   output_YFP   0  (1,2)
NOT         10101010   A1_AmtR      1  (5)
NOR         11000000   S3_SrpR      2  (3,4)
INPUT       00001111   input_pTet   3
INPUT       00110011   input_pBAD   4
INPUT       01010101   input_pTac   5
```

Same three gates, same types, same wiring. `S3_SrpR` / `A1_AmtR` are Cello's names for the
whole gate cassette; SYNBICT names the gate after the repressor CDS inside it.

**Input convention.** The `INPUT` rows above are the convention, per circuit: column 3 of
`input_pTet` / `input_pBAD` / `input_pTac` together with column 3 of `OUTPUT_OR` is one row of
the truth table. `check_results.py` reads them straight from this file, so it never has to
assume a bit order or a polarity. (For the record, the order the circuit *names* are written in
is `j = 4*(1-pBAD) + 2*(1-pTet) + 1*(1-pTac)`; re-indexing every reference table that way
reproduces all 32 names, which is how the convention was confirmed — but the checker does not
depend on it.)

For MD5, compare a *cell* (not a plasmid) against `md5/reference/verilog/scNN.v`: look up the
cell's plasmids in `md5/expected/sc_to_plasmids.json`, merge their gates by repressor CDS
(a repressor is present if **any** of its promoters is active, so the promoter it represses is
`NOR(union of all promoters driving any TU that produces it)`), then compare. The script that
did this in the recorded run is **not in the repo** — see `MANIFEST.md`.

## 5. Known traps

**Leave NMS off.** It is off by default (`-nms` is the switch that turns it *on*; there is no
flag to turn it off), so this only matters if you add it yourself. The Cello library contains
21 `engineered_region` parts that are whole gate *cassettes* (`S3_SrpR` = RBS + ribozyme + CDS
+ terminator). NMS keeps the highest-scoring hit per locus — the cassette — and suppresses
everything nested inside it, **including the terminators**. `gate_assembler` splits
transcriptional units only at terminators, so every unit merges into one:

```
0xEA without -nms : 19 annotations -> 3 gates, 2 wires, truth table 0xEA         ✅
0xEA with    -nms :  9 annotations -> 1 gate,  0 wires, truth table 0xFFFFFFFE   ❌
```

`run_test.sh` deliberately omits it.

**Re-annotating in place.** Running `features_to_circuits.py` on a file that already contains
proteins/interactions raises `SBOL_ERROR_URI_NOT_UNIQUE`. Always start from a freshly
annotated file (`run_test.sh` does).

**The BLAST index is written to the current directory** as `test.*`. `run_test.sh` cd's into
the output dir first so it does not overwrite an index in the repo root.

**`sbol2.Document.write()` validates online by default.** With no network you may see
validator timeouts; set `sbol2.Config.setOption('validate', False)`.

**Run `features_to_circuits.py` by script path, not `python -m`** — the package has no
`__main__`. (`sequences_to_features` does support `-m`.)

## 6. Optional: see it fail

To confirm the checks actually catch a regression, break the annotation step and re-run — in
`sequences_to_features/BlastAligner.py`, delete the four lines

```python
'-task', 'blastn', '-word_size', '9',
```

so BLAST falls back to its megablast default, then

```bash
./run_test.sh --cello out_broken && python check_results.py out_broken/cello --dataset cello
```

Megablast cannot seed the short terminators (47 bp `L3S3P11`), transcriptional units merge,
and about half the circuits develop extra input edges; 9 of them close into a combinational
loop whose truth table is undefined. Expected: roughly 24 circuits still PASS, the rest FAIL
with `DIFFERS` / `undefined (x)`. Restore with
`git checkout -- sequences_to_features/BlastAligner.py`.

## 7. Reporting back

Please send:

1. the `check_results.py` summary table for each dataset you ran,
2. `out/<dataset>/*.log` for any circuit that is not PASS (and not on the known-cases list),
3. `conda list -n synbict_conda | grep -E "sbol2|pysam|blast|pandas"` and `yosys -V`,
4. `git rev-parse HEAD` and `git status --porcelain` from the repo.
