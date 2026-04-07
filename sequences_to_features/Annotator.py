from .FeatureAnnotatorBase import FeatureAnnotatorSimple
from .Feature import Feature
import json, pysam, math
import pandas as pd
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

    def extract_matches(self, min_feature_length=40, exact_match=True):
        blast_output = self.tab_path # error, filename error
        with open(blast_output) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue  # skip headers or blank lines
                segs = line.strip().split('\t')
                if len(segs) < 14:
                    # vsearch output

                    thi = int(segs[9])  # sstart
                    tlo = int(segs[8])
                    ids = int(segs[12])
                    ref_name = segs[1]
                    ref_length = abs(thi - tlo) + 1
                    pid_ref = 100.0 * ids / ref_length
                    if(exact_match):
                        if not (segs[3] == ref_length and math.isclose(pid_ref, 100.0, rel_tol=0.0, abs_tol=1e-6)):
                            continue
                    else:
                        if(pid_ref < 95.0):
                            continue
                    
                else:
                    # change to query start and end
                    ref_name = segs[1]  # ref_name, # sseqid
                    align_len = int(segs[3]) # alignment length (matches+mismatches+gaps)
                    start = int(segs[6]) - 1 # qstart (convert from 1-based to 0-based)
                    end = int(segs[7]) # qend
                    pident = float(segs[2]) # pident
                    ref_length = int(segs[13]) # slen
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
                        if(pid_ref < 95.0):
                            continue

                if ref_length < min_feature_length:
                    continue
                feature = Feature(
                            nucleotides='',
                            identity=ref_name,# will be replaced later
                            roles='',
                            sub_identities='',
                            parent_identities=''
                            #identity=feature_pre['original_identity'],
                            #roles=feature_pre['roles'],
                            #sub_identities=feature_pre.get('sub_identities', []),
                            #parent_identities=feature_pre.get('parent_identities', [])
                            )
                    
                match = ([feature], start, end)
                if start > end:
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

    def extract_matches(self, min_feature_length=40, exact_match=True, is_bowtie2=False): # need change back
        try:
            samfile = pysam.AlignmentFile(self.sam_path, "r")
            for aln in samfile.fetch(until_eof=True):    
                ref_name = aln.reference_name
                length = samfile.get_reference_length(ref_name)  
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
                    if pident < 95.0:
                        continue

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
                    
                match = ([feature], start, end)
                if aln.is_reverse:
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
