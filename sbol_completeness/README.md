# Self-contained SBOL output — evidence files

The annotated SBOL that the aligner paths (`-blastn`, `-bwa`, `-minimap2`) produced did not
contain the library ComponentDefinitions its own SequenceAnnotations referenced. Every
annotation pointed at a SynBioHub URI that was not in the file, so a reader that resolves
those references to get a part's SO role and glyph — SBOLCanvas — had nothing to resolve
and crashed. The FlashText path was never affected.

Cause: `sequences_to_features.py` called the aligner-path annotator with `in_place=True`
pinned, and `FeatureAnnotatorBase.annotate()` derives `copy_definitions` from that flag, so
the matched library parts were never copied into the document that gets written. Fixed by
threading the caller's `in_place` through, which also makes `-p/--in_place` actually work on
those paths.

All files here annotate `test_bundle/cello/sequences/0xEA.fasta` against
`test_bundle/cello/library/cello_library.xml` with `-m 1000 -M 40 -np -ni`.

| File | Produced by | CDs | Annotations | Dangling |
|---|---|---|---|---|
| `0xEA_annotated_blastn_BEFORE_fix.xml` | old code, `-blastn` | 4 | 19 | **16** |
| `0xEA_annotated_blastn_fixed.xml` | fixed code, `-blastn` | 22 | 19 | 0 |
| `0xEA_annotated_blastn_exact.xml` | fixed code, `-blastn -exact` | 19 | 16 | 0 |
| `0xEA_annotated_bwa_fixed.xml` | fixed code, `-bwa` | 22 | 19 | 0 |
| `0xEA_annotated_minimap2_fixed.xml` | fixed code, `-minimap2` | 22 | 19 | 0 |
| `0xEA_annotated_flashtext.xml` | `-flashText` (unaffected either way) | 19 | 9 | 0 |
| `0xEA_annotated_blastn_inplace.xml` | fixed code, `-blastn -p` | 4 | 19 | 16 |

The last row is not a bug: `-p/--in_place` asks for annotation in place, and honouring it is
the point of the fix. It is kept here so the difference is visible.

The three extra definitions in the similar-match aligner files are `S1_SrpR_v1`,
`S2_SrpR_v1` and `S4_SrpR_v1` -- variant definitions that similarity matching creates for
parts not identical to a library entry. All three sit on 410-1258, the same span as the
exactly-matched `S3_SrpR`: they are the other RBS variants of the same SrpR gate cassette.
`-exact` drops them, which is the whole difference between the `_fixed` and `_exact` files
(22 CDs / 19 annotations vs 19 / 16).

Comparing `_exact` with `_flashtext` is more interesting, since **both are exact matching**
and they disagree: 16 annotations against 9. FlashText misses `RiboJ10`, `SrpR`,
`ECK120019600`, `ECK120029600`, `BydvJ`, `AmtR` and `L3S2P55` -- every one of them a part
nested inside a composite cassette (`S3_SrpR` at 410-1258, `A1_AmtR` at 1334-2173). It
annotates the outer cassette but not the sub-parts within it; the aligner annotates both
levels. That is a separate limitation from exact-vs-similar matching.

## Re-checking

```bash
python sbol_completeness/check_selfcontained.py sbol_completeness/*.xml
```

Reports, per file, whether every `component -> definition`, `CD -> sequence` and
`SequenceAnnotation -> component` reference resolves inside the file. Exit code 1 if any
dangles. `0xEA_annotated_blastn_fixed.xml` also passes the official SBOL validator
(`validator.sbolstandard.org`) with `Valid.`

One thing the fix does **not** change: the target construct itself (`_0xEA_comp`) carries no
SO role, before or after. The parts all do. If a viewer still misbehaves on a file whose
references all resolve, that is the next thing to look at.
