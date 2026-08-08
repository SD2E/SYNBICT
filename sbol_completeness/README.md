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

| File | Produced by | ComponentDefinitions | Dangling |
|---|---|---|---|
| `0xEA_annotated_blastn_BEFORE_fix.xml` | old code, `-blastn` | 4 | **16** |
| `0xEA_annotated_blastn_fixed.xml` | fixed code, `-blastn` | 22 | 0 |
| `0xEA_annotated_bwa_fixed.xml` | fixed code, `-bwa` | 22 | 0 |
| `0xEA_annotated_minimap2_fixed.xml` | fixed code, `-minimap2` | 22 | 0 |
| `0xEA_annotated_flashtext.xml` | `-flashText` (unaffected either way) | 19 | 0 |
| `0xEA_annotated_blastn_inplace.xml` | fixed code, `-blastn -p` | 4 | 16 |

The last row is not a bug: `-p/--in_place` asks for annotation in place, and honouring it is
the point of the fix. It is kept here so the difference is visible.

The aligner files carry 3 more definitions than FlashText: `S1_SrpR_v1`, `S2_SrpR_v1` and
`S4_SrpR_v1`, the variant definitions similarity matching creates for parts that are not
identical to a library entry.

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
