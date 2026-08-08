# test_bundle contents

Everything the two features need, for both validation datasets. 4.5 MB total.

```
test_bundle/
├── TESTING.md          how to run and what to expect   <- start here
├── MANIFEST.md         this file
├── run_test.sh         runs both features over a dataset
├── check_results.py    compares a run against expected/ and the references
│
├── cello/              66 Cello circuits, single plasmid each, 3 sensor inputs
│   ├── library/cello_library.xml            146 parts (Cello Eco1C1G1T1)
│   ├── sequences/*.fasta                    66 plasmid sequences  <- pipeline input
│   ├── reference/verilog/*.v                74 Cello reference Verilogs
│   ├── reference/cello_logic/*_A000_*.txt   56 Cello gate structures + input mapping
│   ├── reference/name_map.tsv               which reference belongs to which sequence
│   └── expected/
│       ├── netlists/*.json                  66 expected gate netlists
│       └── logs/*.log                       66 reference run logs (gate table + truth table)
│
└── md5/                70 MD5 plasmids, circular, wired across cells by quorum sensing
    ├── library/Eco2C1G5T1_public_library.xml    251 parts (public MD5 library)
    ├── sequences/*.fasta                        70 plasmid sequences  <- pipeline input
    ├── reference/verilog/*.v                    48 reference Verilogs (gNN_boolean.v = per
    │                                            subcircuit, scNN.v = per cell)
    └── expected/
        ├── netlists/*.json                      70 expected gate netlists
        ├── truth_tables.txt                     per-plasmid truth tables from the reference run
        └── sc_to_plasmids.json                  which plasmids make up each of the 41 cells
```

## Provenance

| What | Where it came from |
|---|---|
| `cello/library/cello_library.xml` | `example/jet_libs/` in the repo — **not tracked by git**, which is why it ships here |
| `cello/sequences/` | `oxea_fastas/`, parsed from `circuit_DNA_sequences_v2.csv` |
| `cello/reference/` | `~/git_repo/cello/resources/tested_circuits/<name>/` |
| `cello/expected/` | a full pipeline run on the current code (43/44 hex circuits reproduce their name) |
| `md5/library/` | `md5_public_run/Eco2C1G5T1_public_library.xml` (public synbiohub.org/public/Eco2C1G5T1) |
| `md5/sequences/` | `~/git_repo/annotate_plasmid/simple_split/` |
| `md5/reference/verilog/` | `~/git_repo/Cello-v2-1-Core/library/verilogs/` |
| `md5/expected/` | `md5_public_run/` (the run recorded in `MD5_validation_record.md`) |

## Known gaps

**`0xB9` is a defective source sequence, not a pipeline failure.** The published CSV row is
missing the 63 bp `pHlyIIR` promoter on the AmtR gate, so it computes `0xF9`. `0xB9_fixed`
is the repaired sequence (promoter re-inserted); it computes `0xB9`. Both ship here.

**Named Cello circuits have no per-name truth table.** `demultiplexer`, `priority_detector`,
`majority`, … are checked on netlist + topology only. `demultiplexer` additionally has one
legitimately isolated gate: it is a 4-output circuit, and the `pAmeR`/`pPhlF`-driven reporter
branches are not present on the published sequence.

**MD5 validation here is per plasmid, not per cell.** The 70 plasmids implement a 2-bit MD5
*partitioned across 41 cells*, wired by quorum-sensing signals, so a single plasmid's netlist
is legitimately a fragment. Validating a whole cell means merging its plasmids by repressor
CDS and comparing to `scNN.v` — the script that did this in the recorded run was kept in a
scratch directory and is **not in the repo**; `sc_to_plasmids.json` and the reference Verilogs
are shipped so it can be rewritten. Per-plasmid netlist regression is what `check_results.py`
covers today.

**Seven MD5 plasmids disagree with the recorded run** (`seq_002, 012, 051, 055, 060, 066,
067`), all on the same ambiguity: `PLuxB` vs `PLux_u42_TA` are near-identical Lux sensor
promoters matching one locus, and neither is a repressed promoter, so the regulation-aware
same-locus collapse cannot choose between them and keeps whichever sorts first. Logic is
unaffected; the promoter name on a primary input is not. Which plasmids fall on which side is
not stable — the self-contained-SBOL fix moved the set from six to seven with no change in
logic. Open issue: the collapse needs a deterministic tie-break.

**Name mapping.** `ANDv2`, `NANDv2`, `XORv2`, `Majority2..6` etc. are v2 builds whose Cello
references are named `AND`, `NAND`, `XOR`, `majority`. `reference/name_map.tsv` records each
mapping; the reference files are also copied under the sequence's own name so tooling can
find them directly.
