#!/usr/bin/env python3
"""Assemble a SYNBICT circuit SBOL into a logic-gate netlist.

SYNBICT's features_to_circuits emits a *discrete* molecular interaction graph
(repression / activation / production motifs + promoter->CDS transcription).
By itself that graph does not say which parts form a gate, what type each gate
is, or how the gates wire together -- so the circuit function can't be read off.

This module closes that gap for Cello-style repressor circuits, where every
gate is stereotyped:

    [input promoter(s)] [RBS] [repressor CDS] [terminator]     (one transcriptional unit)
    repressor protein  --represses-->  its cognate output promoter

so a gate is:
    inputs  = the promoters that drive the repressor CDS (same transcriptional unit)
    cds     = the repressor CDS
    output  = the promoter repressed by the CDS's protein product
    type    = NOR (>=2 inputs) / NOT (1 input); a non-repressing CDS = reporter/OUTPUT

Gates wire when gate A's output promoter is an input promoter of gate B.

Input:  an sbol2.Document holding (a) the annotated construct ComponentDefinition
        (SequenceAnnotations = linear parts with positions + SO roles) and
        (b) the circuit's Interactions (repression/production/activation).
Output: a netlist dict {gates, wires, primary_inputs, primary_outputs}.

This is the *complete-input* assembler (deterministic parse). Missing / unknown
parts (e.g. an unannotated output CDS) are reported as gaps; resolving those is
the job of the later SAT-based inference, not this module.

Assumptions and limitations
---------------------------
The gate model above is correct for standard Cello-style *repressor* logic; it
holds under these assumptions, and simplifies outside them:

1. Repressor gates only. Gate typing covers NOR / NOT (repressor -> promoter) and
   an OR-style reporter OUTPUT. *Activation*-based gates (an activator driving its
   output promoter) are not typed as logic gates -- in Cello these appear only as
   input sensors (primary inputs), which is how they are handled here.
2. Input promoters OR-combine. Several promoters driving one CDS in the same
   transcriptional unit are treated as a logical OR (tandem promoters), so the
   output promoter is NOR(inputs). This is the standard Cello layout assumption.
3. One output promoter per repressor. If a repressor represses more than one
   promoter (fan-out), only the first is used as the gate output; multi-target
   fan-out is not fully modelled.
4. Pure Boolean abstraction. Multimerisation and cooperativity are ignored -- they
   change the analog transfer function (Hill steepness) but not the Boolean truth
   table, so they do not affect the netlist this module produces.
5. Forward strand only (v1). Reverse-complement annotations are skipped.
"""
import json
import logging
import re

import sbol2

logger = logging.getLogger('synbict')

SO_PROMOTER = sbol2.SO_PROMOTER
SO_CDS = sbol2.SO_CDS
SO_TERMINATOR = getattr(sbol2, 'SO_TERMINATOR', 'http://identifiers.org/so/SO:0000141')
SO_RBS = getattr(sbol2, 'SO_RBS', 'http://identifiers.org/so/SO:0000139')  # ribosome entry site

# SBO interaction types
SBO_INHIBITION = sbol2.SBO_INHIBITION          # repression
SBO_STIMULATION = sbol2.SBO_STIMULATION        # activation / transcription
SBO_GENETIC_PRODUCTION = 'http://identifiers.org/biomodels.sbo/SBO:0000589'
# SBO participant roles
SBO_INHIBITOR = sbol2.SBO_INHIBITOR
SBO_INHIBITED = sbol2.SBO_INHIBITED
SBO_STIMULATOR = sbol2.SBO_STIMULATOR
SBO_STIMULATED = sbol2.SBO_STIMULATED
SBO_PRODUCT = 'http://identifiers.org/biomodels.sbo/SBO:0000011'
SBO_TEMPLATE = 'http://identifiers.org/biomodels.sbo/SBO:0000645'
SBO_REACTANT = 'http://identifiers.org/biomodels.sbo/SBO:0000010'
# some converters (e.g. the MD5Collection SBOL) tag the production template with
# the SBOL role http://sbols.org/v2#template instead of the SBO template role.
SBOL_TEMPLATE = 'http://sbols.org/v2#template'


def _short(identity):
    """Short human name from an SBOL identity URI (…/<name>/<version>)."""
    if not identity:
        return identity
    parts = identity.rstrip('/').split('/')
    # identities look like http://ns/<name>/<version>; take the name segment
    return parts[-2] if len(parts) >= 2 and parts[-1].isdigit() else parts[-1]


# Name patterns for Cello parts whose ComponentDefinition is often dropped from
# the circuit SBOL (so their SO role can't be read) -- fall back to the name.
# B001x are terminator BioBricks (B0015 ...); B003x/B006x are RBSs -- do NOT match those.
_TERMINATOR_RE = re.compile(r'^(ECK\d|L\dS\d|DT\d|Ret\d|.*term|rrnB|B001\d|T\d)', re.IGNORECASE)
_RBS_RE = re.compile(r'^(B00[36]\d|Rj\d|RBS\d|.*rbs)', re.IGNORECASE)
_RIBOZYME_RE = re.compile(r'(RiboJ|BydvJ|PlmJ|ScmJ|SarJ|ElvJ|LtsvJ|.*J\d*$)')
# Selection/resistance markers and origins are not part of the logic; a marker CDS
# is non-repressing so it would otherwise surface as a spurious OUTPUT gate. Only
# applied to a CDS that represses nothing, so real repressors (TetR, ...) are safe.
_MARKER_RE = re.compile(r'^(amp[R]?\d*|bla|kan[R]?\d*|neoR|cam[R]?\d*|cat[R]?\d*|'
                        r'cml|chlor|spec[R]?\d*|gent[R]?\d*|ori\b|rop\b|repA)',
                        re.IGNORECASE)


def _norm_prom(name):
    """Canonicalize a promoter name so annotation variants of the same promoter
    compare equal: PJR1_variant_1 / PJR3b / PJR3_2b -> PJR1 / PJR3 / PJR3."""
    name = re.sub(r'(_variant_\d+|_v\d+)$', '', name or '')
    name = re.sub(r'(PJR\d+).*', r'\1', name)
    return name


class GateAssembler:
    def __init__(self, doc, circuit_doc=None, library_docs=None):
        self.doc = doc                       # construct + part definitions (annotated SBOL)
        self.circuit_doc = circuit_doc or doc # interactions (features_to_circuits output)
        # part libraries (the -c sub-circuit files). The annotated construct only
        # references a part's definition by URI without embedding it, so the SO role
        # (e.g. terminator) lives here -- consulted before the name-based fallback.
        self.library_docs = list(library_docs) if library_docs else []
        self.kind_cache = {}      # key -> 'promoter'|'cds'|'terminator'|'protein'|'unknown'
        # regulation maps keyed by definition identity
        self.cds_to_protein = {}  # repressor CDS -> its protein product
        self.protein_to_repressed = {}   # repressor protein -> [repressed promoters]
        self.protein_to_activated = {}   # activator protein -> [activated promoters]
        # every promoter that is repressed by some regulator, harvested from the
        # construct AND the part libraries. Used to disambiguate overlapping same-
        # locus promoter annotations: the regulated one is the real wired promoter.
        self.repressed_promoter_ids = set()
        self.gaps = []            # human-readable notes about missing info

    # ------------------------------------------------------------------ kinds
    def _kind(self, identity, name=None):
        key = (identity, name)
        if key in self.kind_cache:
            return self.kind_cache[key]
        kind = 'unknown'
        # 1) authoritative: SO role / biopax type from the ComponentDefinition.
        # search the construct/circuit docs first, then the part libraries (where a
        # by-URI-referenced definition and its SO role actually live).
        for src in (self.doc, self.circuit_doc, *self.library_docs):
            try:
                cd = src.getComponentDefinition(identity)
            except Exception:
                continue
            roles, types = set(cd.roles), set(cd.types)
            if SO_PROMOTER in roles:
                kind = 'promoter'
            elif SO_CDS in roles:
                kind = 'cds'
            elif SO_TERMINATOR in roles:
                kind = 'terminator'
            elif SO_RBS in roles:
                kind = 'ribozyme'   # an RBS is an insulator inside a TU, not a splitter
            elif sbol2.BIOPAX_PROTEIN in types:
                kind = 'protein'
            if kind != 'unknown':
                break
        # 2) fall back to name heuristics (definition dropped from the SBOL)
        if kind == 'unknown':
            nm = name or _short(identity) or ''
            if nm.endswith('_protein'):
                kind = 'protein'
            elif _RBS_RE.match(nm):
                kind = 'ribozyme'      # RBS: insulator inside a TU, checked before terminator
            elif _TERMINATOR_RE.match(nm):
                kind = 'terminator'
            elif _RIBOZYME_RE.search(nm):
                kind = 'ribozyme'      # insulator: sits inside a TU, does not split it
            elif nm[:1] == 'p' and nm[1:2].isupper():
                kind = 'promoter'
        self.kind_cache[key] = kind
        return kind

    # ------------------------------------------------------- parse regulation
    def parse_regulation(self):
        """Read every Interaction and record production / repression / activation
        by *definition identity*, distinguishing real regulation (protein acts on
        promoter) from transcription stimulation (promoter acts on CDS)."""
        for md in self.circuit_doc.moduleDefinitions:
            fc_def = {fc.identity: fc.definition for fc in md.functionalComponents}

            for intxn in md.interactions:
                types = set(intxn.types)
                # collect (role, definition) per participation
                parts = []
                for par in intxn.participations:
                    dfn = fc_def.get(par.participant)
                    if dfn:
                        parts.append((set(par.roles), dfn))

                if SBO_GENETIC_PRODUCTION in types:
                    template = next((d for r, d in parts
                                     if SBO_TEMPLATE in r or SBO_REACTANT in r
                                     or SBOL_TEMPLATE in r), None)
                    product = next((d for r, d in parts if SBO_PRODUCT in r), None)
                    if template and product:
                        self.cds_to_protein[template] = product

                elif SBO_INHIBITION in types:
                    inhibitor = next((d for r, d in parts if SBO_INHIBITOR in r), None)
                    inhibited = next((d for r, d in parts if SBO_INHIBITED in r), None)
                    if inhibitor and inhibited:
                        self.protein_to_repressed.setdefault(inhibitor, []).append(inhibited)
                        self.repressed_promoter_ids.add(_norm_prom(_short(inhibited)))

                elif SBO_STIMULATION in types:
                    stim = next((d for r, d in parts if SBO_STIMULATOR in r), None)
                    stimd = next((d for r, d in parts if SBO_STIMULATED in r), None)
                    # keep only protein->promoter activation; drop promoter->CDS
                    # transcription stimulation (that is layout, not regulation)
                    if stim and stimd and self._kind(stim) == 'protein' \
                            and self._kind(stimd) == 'promoter':
                        self.protein_to_activated.setdefault(stim, []).append(stimd)

        # also learn which promoters are repressed in the part libraries, so an
        # overlapping same-locus promoter that IS wired (e.g. PJR1) can be told
        # apart from a spurious constitutive match (e.g. Plambda) at that locus.
        for lib in self.library_docs:
            for md in lib.moduleDefinitions:
                fc_def = {fc.identity: fc.definition for fc in md.functionalComponents}
                for intxn in md.interactions:
                    if SBO_INHIBITION not in set(intxn.types):
                        continue
                    for par in intxn.participations:
                        if SBO_INHIBITED in set(par.roles):
                            dfn = fc_def.get(par.participant)
                            if dfn:
                                self.repressed_promoter_ids.add(_norm_prom(_short(dfn)))

    def _cds_output_promoter(self, cds_identity):
        """The promoter a repressor CDS ultimately represses (CDS -> protein ->
        repressed promoter). None if the CDS represses nothing (reporter)."""
        protein = self.cds_to_protein.get(cds_identity)
        if not protein:
            return None
        repressed = self.protein_to_repressed.get(protein)
        return repressed[0] if repressed else None

    # -------------------------------------------------------- parse construct
    def _construct_parts(self):
        """Return the construct's inline parts as sorted dicts by position,
        collapsing overlapping same-locus CDS annotations (variant duplicates)."""
        construct = next((cd for cd in self.doc.componentDefinitions
                          if len(cd.sequenceAnnotations) > 0), None)
        if construct is None:
            raise ValueError('no annotated construct (ComponentDefinition with '
                             'SequenceAnnotations) found in document')

        raw = []
        for sa in construct.sequenceAnnotations:
            if len(sa.locations) != 1 or sa.locations[0].getTypeURI() != sbol2.SBOL_RANGE:
                continue
            rng = sa.locations.getRange()
            if rng.orientation != sbol2.SBOL_ORIENTATION_INLINE:
                continue  # v1: forward strand only
            # resolve the part definition behind this annotation
            definition = None
            if sa.component:
                try:
                    definition = construct.components.get(sa.component).definition
                except Exception:
                    definition = None
            nm = _short(definition) if definition else (sa.name or sa.displayId)
            raw.append({'start': rng.start, 'end': rng.end,
                        'definition': definition,
                        'name': nm,
                        'kind': self._kind(definition, nm)})

        raw.sort(key=lambda p: (p['start'], p['end']))

        # collapse overlapping same-locus variant annotations into one node
        parts = []
        for p in raw:
            if parts and p['kind'] == parts[-1]['kind'] and p['start'] < parts[-1]['end']:
                if p['kind'] == 'cds':
                    # prefer the CDS that has a known output promoter
                    if self._cds_output_promoter(p['definition']) and \
                            not self._cds_output_promoter(parts[-1]['definition']):
                        parts[-1] = p
                    continue
                if p['kind'] == 'promoter':
                    # prefer the promoter that is actually regulated (wired), so a
                    # spurious constitutive match at the same locus is dropped
                    p_reg = _norm_prom(p['name']) in self.repressed_promoter_ids
                    prev_reg = _norm_prom(parts[-1]['name']) in self.repressed_promoter_ids
                    if p_reg and not prev_reg:
                        parts[-1] = p
                    continue
            parts.append(p)
        return construct, parts

    # ------------------------------------------------------------- segment TUs
    @staticmethod
    def _segment(parts):
        """Split the linear part list into transcriptional units at terminators."""
        tus, cur = [], []
        for p in parts:
            if p['kind'] == 'terminator':
                if cur:
                    tus.append(cur)
                    cur = []
                continue
            if p['kind'] in ('promoter', 'cds'):
                cur.append(p)
        if cur:
            tus.append(cur)
        return tus

    # --------------------------------------------------------------- assemble
    def assemble(self):
        self.parse_regulation()
        construct, parts = self._construct_parts()
        tus = self._segment(parts)

        gates = []
        for idx, tu in enumerate(tus):
            promoters = [p['name'] for p in tu if p['kind'] == 'promoter']
            cds_list = [p for p in tu if p['kind'] == 'cds']
            if not cds_list:
                # promoter-only TU with no CDS: its promoters just drive whatever
                # comes next; nothing to assemble, but note it as a gap
                if promoters:
                    self.gaps.append(f'transcriptional unit with promoters '
                                     f'{promoters} has no annotated CDS')
                continue
            for cds in cds_list:
                out_prom = self._cds_output_promoter(cds['definition'])
                is_output = out_prom is None
                # a non-repressing selection/resistance marker (ampR, kanR, ori...)
                # is not circuit logic -- drop it instead of emitting a stray OUTPUT
                if is_output and _MARKER_RE.match(cds['name'] or ''):
                    continue
                gate = {
                    'id': f'g{len(gates) + 1}',
                    'cds': cds['name'],
                    'inputs': promoters,
                    'output': (_short(out_prom) if out_prom else cds['name']),
                    'type': ('OUTPUT' if is_output
                             else ('NOR' if len(promoters) >= 2 else 'NOT')),
                    'is_output': is_output,
                }
                gates.append(gate)

        # wiring: promoter produced by a gate -> feeds any gate that lists it as input
        producer = {g['output']: g['id'] for g in gates if not g['is_output']}
        wires = []
        for g in gates:
            for inp in g['inputs']:
                if inp in producer and producer[inp] != g['id']:
                    wires.append([producer[inp], g['id']])

        consumed = {inp for g in gates for inp in g['inputs']}
        primary_inputs = sorted({inp for g in gates for inp in g['inputs']
                                 if inp not in producer})
        primary_outputs = [g['output'] for g in gates
                           if g['is_output'] or g['output'] not in consumed]

        return {
            'construct': _short(construct.identity),
            'gates': gates,
            'wires': wires,
            'primary_inputs': primary_inputs,
            'primary_outputs': primary_outputs,
            'gaps': self.gaps,
        }


def assemble_from_document(doc):
    return GateAssembler(doc).assemble()


def assemble_from_file(path):
    doc = sbol2.Document()
    doc.read(path)
    return GateAssembler(doc).assemble()


def write_netlist(netlist, path):
    with open(path, 'w') as f:
        json.dump(netlist, f, indent=2)
    logger.info('Wrote gate netlist %s', path)


def summarize(netlist):
    print(f"construct: {netlist['construct']}")
    print(f"primary inputs : {netlist['primary_inputs']}")
    print(f"primary outputs: {netlist['primary_outputs']}")
    print("gates:")
    for g in netlist['gates']:
        print(f"  {g['id']:4s} {g['type']:6s} inputs={g['inputs']} "
              f"cds={g['cds']} -> output={g['output']}")
    print(f"wires: {netlist['wires']}")
    if netlist['gaps']:
        print("gaps (need SAT / re-annotation):")
        for gap in netlist['gaps']:
            print(f"  - {gap}")


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='SYNBICT circuit SBOL -> gate netlist')
    ap.add_argument('input', help='circuit SBOL2 XML (features_to_circuits output)')
    ap.add_argument('-o', '--out', help='write netlist JSON to this path')
    args = ap.parse_args()

    nl = assemble_from_file(args.input)
    summarize(nl)
    if args.out:
        write_netlist(nl, args.out)
