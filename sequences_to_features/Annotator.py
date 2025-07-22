# write an Annotator python class, that read sam file and convert the alignment result and find the coresponding metadata dictionary from FeatureExtractor class or from the saved json file to find the metadata and wrtie all this insert into a SBOL file.
from FeatureAnnotatorBase import FeatureAnnotatorSimple
from sequences_to_features import Feature
import json, pysam

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

    def extract_matches(self, exact_match=True):
        try:
            samfile = pysam.AlignmentFile(self.sam_path, "r")
            for read in samfile.fetch(until_eof=True):
                if(exact_match):
                    if not (read.has_tag("NM") and read.get_tag("NM") == 0):
                        continue
                reference_name, start, end = self.parse_cigar_for_query_coords(read)

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
                if read.is_reverse:
                    self.rc_matches.append(match)
                else:
                    self.inline_matches.append(match)
        except Exception as e:
            print("Failed to process read:", read.query_name, "Error:", e)

        return self.inline_matches, self.rc_matches

    def insert_into_sbol(self, target_library, min_target_length=20, in_place=False,
                         output_library=None, complete_matches=False, strip_prefixes=[], output_matches=False):
        feature_annotater = FeatureAnnotatorSimple(self.inline_matches, self.rc_matches)
        return feature_annotater.annotate(self.inline_matches, self.rc_matches,
            target_library, min_target_length, in_place, output_library, complete_matches,
                 strip_prefixes, output_matches)
        
    def write_sbol_to_file(target_library):
        target_library.write('output_file.xml')
