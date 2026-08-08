#!/usr/bin/env python3
"""Report whether an SBOL file is self-contained.

A reader such as SBOLCanvas follows a SequenceAnnotation to its Component and then to
that Component's `definition` URI to get the part's SO role, name and sequence. If the
ComponentDefinition behind that URI is not in the file, the reference dangles and the
reader has nothing to resolve -- which is what made SBOLCanvas fail on annotated output
from the aligner paths (-blastn / -bwa / -minimap2) before the in_place fix.

Usage: python check_selfcontained.py <file.xml> [<file.xml> ...]
Exit code 1 if any file has a dangling reference.
"""
import sys, os
import sbol2

sbol2.Config.setOption('validate', False)


def audit(path):
    doc = sbol2.Document()
    doc.read(path)
    cds = {cd.identity for cd in doc.componentDefinitions}
    seqs = {s.identity for s in doc.sequences}
    dangling = {'component -> definition': [], 'CD -> sequence': [], 'SA -> component': []}
    for cd in doc.componentDefinitions:
        local = {c.identity for c in cd.components}
        for c in cd.components:
            if c.definition not in cds:
                dangling['component -> definition'].append(c.definition)
        for s in cd.sequences:
            if s not in seqs:
                dangling['CD -> sequence'].append(s)
        for sa in cd.sequenceAnnotations:
            if sa.component and sa.component not in local:
                dangling['SA -> component'].append(sa.identity)
    total = sum(len(v) for v in dangling.values())
    print(f'{os.path.basename(path)}')
    print(f'   ComponentDefinition {len(doc.componentDefinitions)}, Sequence {len(doc.sequences)}')
    for kind, refs in dangling.items():
        mark = 'ok' if not refs else f'{len(refs)} DANGLING'
        print(f'   {kind:<24} {mark}')
        for r in sorted(set(refs))[:5]:
            print(f'        {r}')
    print(f'   => {"self-contained" if total == 0 else "NOT self-contained"}\n')
    return total


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    sys.exit(1 if sum(audit(p) for p in sys.argv[1:]) else 0)
