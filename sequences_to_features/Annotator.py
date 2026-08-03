from .FeatureAnnotatorBase import FeatureAnnotatorSimple
from .Feature import Feature
import json, pysam, math
import pandas as pd
def canonical_intervals(start, end, target_length=None):
    """Map a half-open query interval onto the target, as a list of intervals.

    For a linear target this is just ``[(start, end)]``. For a circular target
    the query was extended to ``seq + seq[:overlap]``, so a hit can sit in the
    appended copy (``start >= target_length``) or straddle the origin
    (``end > target_length``); both are folded back onto ``[0, target_length)``,
    the second one as two intervals. A hit longer than the whole target covers
    the entire circle.
    """
    if not target_length:
        return [(start, end)]

    span = end - start
    if span >= target_length:
        return [(0, target_length)]

    s = start % target_length
    e = s + span
    if e <= target_length:
        return [(s, e)]
    return [(s, target_length), (0, e - target_length)]


def intervals_overlap(a_intervals, b_intervals):
    """Total number of positions covered by both interval lists."""
    total = 0
    for a_start, a_end in a_intervals:
        for b_start, b_end in b_intervals:
            ov = min(a_end, b_end) - max(a_start, b_start)
            if ov > 0:
                total += ov
    return total


def apply_non_maximum_suppression(candidates, overlap_frac=0.5, target_length=None):
    """Non-maximum suppression over query coordinates.

    `candidates` is a sequence of tuples whose first three elements are
    ``(score, start, end)``; anything after that is carried through untouched, so
    every mapper can use its own payload. Sorts by score (desc) and keeps a hit
    only if it does not substantially overlap a higher-scoring one. Adjacent
    parts (promoter/RBS/CDS) barely overlap and all survive; competing
    annotations for one locus collapse to the single best-scoring reference.

    Suppression is only against a STRICTLY higher-scoring overlap. Tied hits
    (e.g. PJR1 62bp vs Plambda 58bp, identical bitscore) are all kept so the
    regulation-aware promoter collapse downstream picks the right one instead of
    a coin-flip tie-break.

    `target_length` (the length of the ORIGINAL, un-extended target) switches
    overlap testing to circular coordinates. Without it, a backbone feature that
    spans the origin -- aligned as one block ending past `target_length` -- lies
    entirely to the right of the short parts near the 5' end and so never
    suppresses them, while the same feature does suppress their counterparts at
    the 3' end. Folding both onto the circle first makes the two ends behave
    identically. Duplicates from the appended origin overlap fold onto their
    first-copy twins and score identically, so the tie rule keeps both;
    normalize_circular_matches drops the appended copy afterwards.

    Off by default at every call site: NMS discards nested parts, which is wanted
    for circuit reconstruction but not for exhaustive annotation.
    """
    ordered = sorted(candidates, key=lambda c: -c[0])
    kept = []
    kept_spans = []  # (score, intervals, length) parallel to `kept`
    for cand in ordered:
        cscore, s, e = cand[0], cand[1], cand[2]
        c_intervals = canonical_intervals(s, e, target_length)
        c_length = e - s
        conflict = False
        for k_score, k_intervals, k_length in kept_spans:
            ov = intervals_overlap(c_intervals, k_intervals)
            if ov > 0 and ov >= overlap_frac * min(c_length, k_length) and k_score > cscore:
                conflict = True
                break
        if not conflict:
            kept.append(cand)
            kept_spans.append((cscore, c_intervals, c_length))
    return kept


def suppress_short_matches(kept_matches, short_inline, short_rc, overlap_frac=0.5,
                           target_length=None):
    """Extend non-maximum suppression to the exhaustive short-feature matches.

    The seed-based mappers run NMS over their own candidates before the
    ShortFeatureMatcher (and Prokka) hits are merged in, so a short part nested
    inside a longer aligner hit used to survive unconditionally -- NMS collapsed
    the >=14 bp parts at a locus but left the <14 bp ones on top of them.

    `kept_matches` (the aligner/Prokka matches, both strands, already resolved
    among themselves) are fixed and act only as suppressors; they are never
    dropped here. A short match is discarded when a strictly longer match --
    aligner hit or another short hit -- covers at least `overlap_frac` of it.
    Ties are kept, matching apply_non_maximum_suppression, so synonymous parts of
    equal length all survive for the downstream collapse to arbitrate.

    Matches use the shared ``(feature_list, start, end)`` tuple format; strands
    are pooled for the overlap test (as in the mappers, where inline and rc
    candidates compete in one list) and split again on return.

    Returns (short_inline, short_rc) filtered, in their original order.
    """
    fixed_spans = [(canonical_intervals(m[1], m[2], target_length), m[2] - m[1])
                   for m in kept_matches]

    tagged = [(m, True) for m in short_inline] + [(m, False) for m in short_rc]
    # Longest first: within the short set the longer part wins, as bitscore does
    # for the aligner candidates.
    tagged.sort(key=lambda t: -(t[0][2] - t[0][1]))

    kept_inline_ids = set()
    kept_rc_ids = set()
    kept_spans = list(fixed_spans)
    for match, is_inline in tagged:
        c_intervals = canonical_intervals(match[1], match[2], target_length)
        c_length = match[2] - match[1]
        conflict = False
        for k_intervals, k_length in kept_spans:
            if k_length <= c_length:
                continue  # only a STRICTLY longer match suppresses
            ov = intervals_overlap(c_intervals, k_intervals)
            if ov > 0 and ov >= overlap_frac * min(c_length, k_length):
                conflict = True
                break
        if not conflict:
            (kept_inline_ids if is_inline else kept_rc_ids).add(id(match))
            kept_spans.append((c_intervals, c_length))

    return ([m for m in short_inline if id(m) in kept_inline_ids],
            [m for m in short_rc if id(m) in kept_rc_ids])


class ProkkaTableFeatureMapper:
    def __init__(self):
        self.inline_matches = []
        self.rc_matches = []
        
    def extend_list(self, dna_matches, protein_matches):
        """
        Extend dna_matches with protein_matches whose identities are not already present. Prioritizes protein matches.
        """
        existing_identities = {
            match[0][0].identity for match in protein_matches
        }

        for match in dna_matches:
            identity = match[0][0].identity
            if identity not in existing_identities or identity == " ":
                protein_matches.append(match)
                existing_identities.add(identity)  # prevent future duplicates

        return protein_matches


    def extract_matches(self, cds_df: pd.DataFrame, min_feature_length=40, mode='exact'):
        """
        Build inline_matches / rc_matches from Prokka-parsed CDS dataframe.

        mode='exact'   -> keep only identity_pct == 100.0 (exact protein match)
        mode='similar' -> keep all identity_pct < 100%, exclude "hypothetical protein"
        mode='all'     -> keep every row regardless of identity or product name
        """

        # reset containers each call
        self.inline_matches = []
        self.rc_matches = []

        if cds_df is None or len(cds_df) == 0:
            return self.inline_matches, self.rc_matches

        def _is_missing(x) -> bool:
            # True for None, NaN, and empty string
            if x is None:
                return True
            if isinstance(x, float) and math.isnan(x):
                return True
            if isinstance(x, str) and x.strip() == "":
                return True
            return False

        def _is_hypothetical(row) -> bool:
            product = row.get("product", None)
            if _is_missing(product):
                return True
            return "hypothetical" in str(product).lower()

        for _, row in cds_df.iterrows():
            start = int(row["start"])
            end = int(row["end"])
            strand = str(row.get("strand", "+")).strip()
            type = row["type"]

            feature_len = abs(end - start) + 1
            if feature_len < min_feature_length:
                continue

            identity_pct = row.get("identity_pct", None)

            if mode == 'exact':
                # keep only 100% protein identity
                if _is_missing(identity_pct):
                    continue
                try:
                    if not math.isclose(float(identity_pct), 100.0, rel_tol=0.0, abs_tol=1e-6):
                        continue
                except (TypeError, ValueError):
                    continue
            elif mode == 'similar':
                # keep any identity_pct (including 100%), exclude hypothetical proteins
                if _is_hypothetical(row):
                    continue
            # mode == 'all': keep everything, no filtering

            # (3) ref_name is ids_sequence, or None if missing
            ids_sequence = row.get("ids_sequence", None)
            ref_name = " " if _is_missing(ids_sequence) else str(ids_sequence)

            feature = Feature(
                nucleotides="",
                identity=ref_name,   # NOTE: can be " ""
                roles="",
                sub_identities="",
                parent_identities="",
            )
            # prokka has one more results
            match = ([feature], start, end, identity_pct, type)

            # Use strand to decide rc vs inline
            if strand == "-":
                self.rc_matches.append(match)
            else:
                self.inline_matches.append(match)

        return self.inline_matches, self.rc_matches

class TableFeatureMapper:
    def __init__(self, tab_path, min_mapq=20):
        self.tab_path = tab_path
        self.min_mapq = min_mapq
        #self.metadata_dict = self._load_metadata(metadata_path)
        self.inline_matches = []
        self.rc_matches = []

    def _load_metadata(self, metadata_path):
        with open(metadata_path, "r") as f:
            return json.load(f)

    def extract_matches(self, min_feature_length=40, exact_match=True,
                        pid_threshold=95.0, overlap_frac=0.5, apply_nms=False,
                        target_length=None):
        """Parse the blast (or vsearch) tabular output and return the best-scoring,
        non-overlapping set of feature matches.

        Every hit that clears the identity / length filters is collected as a
        candidate; overlapping candidates for the same locus are then resolved by
        keeping the highest-scoring one (bitscore for blast, identities for
        vsearch). This replaces the previous "emit every hit >= 95%" behaviour,
        which annotated a shorter partial part over the correct longer one when both
        passed threshold. Lower `min_feature_length` to recover short parts (RBS,
        short terminators); `pid_threshold` is coverage-weighted identity
        (identical bases / reference length). Pass `target_length` (the
        un-extended target length) for circular targets so suppression happens in
        circular coordinates."""
        blast_output = self.tab_path
        candidates = []  # (score, start, end, ref_name, sstart, send)
        with open(blast_output) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue  # skip headers or blank lines
                segs = line.strip().split('\t')
                sstart = send = 0  # subject coords; set in the blast branch below
                if len(segs) < 14:
                    # vsearch output
                    thi = int(segs[9])  # sstart
                    tlo = int(segs[8])
                    ids = int(segs[12])
                    ref_name = segs[1]
                    ref_length = abs(thi - tlo) + 1
                    pid_ref = 100.0 * ids / ref_length
                    start = int(segs[6]) - 1
                    end = int(segs[7])
                    score = ids
                    if(exact_match):
                        if not (int(segs[3]) == ref_length and math.isclose(pid_ref, 100.0, rel_tol=0.0, abs_tol=1e-6)):
                            continue
                    else:
                        if(pid_ref < pid_threshold):
                            continue

                else:
                    # change to query start and end
                    ref_name = segs[1]  # ref_name, # sseqid
                    align_len = int(segs[3]) # alignment length (matches+mismatches+gaps)
                    start = int(segs[6]) - 1 # qstart (convert from 1-based to 0-based)
                    end = int(segs[7]) # qend
                    # subject (feature) coords: blastn keeps the query forward and
                    # reports sstart>send for minus-strand hits, so strand must be
                    # read from the subject, not the (always-ascending) query.
                    sstart = int(segs[8]) # sstart
                    send = int(segs[9]) # send
                    pident = float(segs[2]) # pident
                    ref_length = int(segs[13]) # slen
                    score = float(segs[11]) # bitscore -- used to rank overlapping hits
                    if exact_match:
                        # Check for exact match, e.g., if the alignment length matches the reference length
                        if not ((end - start) == ref_length
                                and math.isclose(pident, 100.0, rel_tol=0.0, abs_tol=1e-6)):
                            continue
                    else:
                        if(len(segs) > 14 and segs[14] is not None): # improve
                            nident = int(segs[14]) # nident
                            pid_ref = 100.0 * nident / ref_length
                        else:
                            pid_ref = 100.0 * align_len / ref_length # estimate
                        if(pid_ref < pid_threshold):
                            continue

                if ref_length < min_feature_length:
                    continue
                candidates.append((score, start, end, ref_name, sstart, send))

        kept = (apply_non_maximum_suppression(candidates, overlap_frac, target_length)
                if apply_nms else candidates)

        for score, start, end, ref_name, sstart, send in kept:
            feature = Feature(
                        nucleotides='',
                        identity=ref_name,# will be replaced later
                        roles='',
                        sub_identities='',
                        parent_identities=''
                        )
            match = ([feature], start, end)
            # route by subject strand (sstart>send == minus). The previous
            # `start > end` test used the query, which blastn always reports
            # ascending, so rc_matches was never populated and every hit was
            # annotated inline regardless of its true orientation.
            if sstart > send:
                self.rc_matches.append(match)
            else:
                self.inline_matches.append(match)
        return self.inline_matches, self.rc_matches
    
class SAMFeatureMapper:
    def __init__(self, sam_path, min_mapq=20):
        self.sam_path = sam_path
        self.min_mapq = min_mapq
        #self.metadata_dict = self._load_metadata(metadata_path)
        self.inline_matches = []
        self.rc_matches = []

    def _load_metadata(self, metadata_path):
        with open(metadata_path, "r") as f:
            return json.load(f)

    def parse_cigar_for_query_coords(self, read):
        cigar_tuples = read.cigartuples
        query_len = 0
        query_consuming_ops = {0, 1, 7, 8}  # M, I, =, X

        for op, length in cigar_tuples:
            if op in query_consuming_ops:
                query_len += length

        hard_clip_front = cigar_tuples[0][1] if cigar_tuples[0][0] in {4, 5} else 0
        hard_clip_end = cigar_tuples[-1][1] if cigar_tuples[-1][0] in {4, 5} else 0

        full_query_len = hard_clip_front + query_len + hard_clip_end
        if read.is_reverse:
            query_end = full_query_len - hard_clip_front
            query_start = query_end - query_len
        else:
            query_start = hard_clip_front
            query_end = query_start + query_len

        return read.reference_name, query_start, query_end

    def extract_matches(self, min_feature_length=40, exact_match=True, is_bowtie2=False,
                        pid_threshold=95.0, overlap_frac=0.5, apply_nms=False,
                        target_length=None):
        """Parse the SAM output and return the feature matches.

        `pid_threshold` was previously hardcoded to 95.0 here, so the CLI/API knob
        had no effect on the BWA and Minimap2 paths. `apply_nms` uses the same
        suppression as the tabular path, ranked by number of identical bases
        (by reference length for exact matches, where every hit is 100%).
        `target_length` (the un-extended target length) puts that suppression in
        circular coordinates for circular targets.
        """
        candidates = []  # (score, start, end, ref_name, is_reverse)
        try:
            samfile = pysam.AlignmentFile(self.sam_path, "r")
            for aln in samfile.fetch(until_eof=True):
                ref_name = aln.reference_name
                if aln.is_unmapped or ref_name is None:
                    continue
                length = samfile.get_reference_length(ref_name)
                # enforce the minimum feature length (TableFeatureMapper does this
                # too); without it short RBS/terminator/promoter refs leak through.
                if length < min_feature_length:
                    continue
                # Ranking for optional NMS. Every exact hit is 100% identical,
                # so those rank by reference length -- the longer part wins.
                # Similar hits override this with their identical-base count.
                score = length
                if exact_match:
                    # Check for exact match AS = len(ref)
                    if is_bowtie2:
                        print("yes bowtie2")
                        if not (aln.get_tag("XO") == 0 and aln.get_tag("XM") == 0 and aln.cigartuples[1][1] == length):
                            continue
                    else:
                        if not (aln.has_tag("NM") and aln.get_tag("NM") == 0 and aln.cigartuples[1][1] == length):
                            continue
                else:
                    # if identity >= threshold 
                    M = I = D = EQ = X = 0
                    for op, ln in (aln.cigartuples or []):
                        if op == 0: M += ln      # M (match+mismatch)
                        elif op == 1: I += ln;   # insertion (run)
                        elif op == 2: D += ln;   # deletion (run)
                        elif op == 7: EQ += ln   # '=' exact match
                        elif op == 8: X  += ln   # 'X' mismatch

                    block = (EQ + X) if (EQ + X) > 0 else M
                    
                    # mismatches:
                    if (EQ + X) > 0:
                        mismatches = X
                        matches = EQ
                    else:
                        NM = aln.get_tag("NM")
                        mismatches = max(NM - I - D, 0)
                        matches = max(block - mismatches, 0)

                    # BLAST alignment length includes gaps (I + D)
                    aln_len = block + I + D
                    if aln_len == 0:
                        continue
                    pident = 100.0 * matches / (length + I + D)
                    if pident < pid_threshold:
                        continue
                    # Rank similar hits by number of identical bases.
                    score = matches

                reference_name, start, end = self.parse_cigar_for_query_coords(aln)

                #print("annotation: ", reference_name, start, end)
                #ref_name = samfile.get_reference_name(read.reference_id)
                #feature_pre = self.metadata_dict.get(ref_name)
                #print("feature_pre: ", feature_pre)
                # construct feature from metadata dictionary only for the mapped parts
                feature = Feature(
                    nucleotides='',
                    identity=reference_name,# will be replaced later
                    roles='',
                    sub_identities='',
                    parent_identities=''
                    #identity=feature_pre['original_identity'],
                    #roles=feature_pre['roles'],
                    #sub_identities=feature_pre.get('sub_identities', []),
                    #parent_identities=feature_pre.get('parent_identities', [])
                    )
                    
                candidates.append((score, start, end, feature, aln.is_reverse))

            kept = (apply_non_maximum_suppression(candidates, overlap_frac, target_length)
                    if apply_nms else candidates)

            for _score, start, end, feature, is_reverse in kept:
                match = ([feature], start, end)
                if is_reverse:
                    self.rc_matches.append(match)
                else:
                    self.inline_matches.append(match)
        except Exception as e:
            print("Failed to process alignment","Error:", e)

        return self.inline_matches, self.rc_matches

    def insert_into_sbol(self, target_library, min_target_length=20, in_place=False,
                         output_library=None, complete_matches=False, strip_prefixes=[], output_matches=False):
        feature_annotater = FeatureAnnotatorSimple(self.inline_matches, self.rc_matches)
        return feature_annotater.annotate(self.inline_matches, self.rc_matches,
            target_library, min_target_length, in_place, output_library, complete_matches,
                 strip_prefixes, output_matches)
        
    def write_sbol_to_file(target_library):
        target_library.write('output_file.xml')
