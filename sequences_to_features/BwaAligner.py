import subprocess
import tempfile
from .sbol_utils import sbol_sequence
from .Aligner import Aligner


class BwaAligner(Aligner):
    def __init__(self, index_prefix):
        super().__init__(index_prefix)

    def align(self, query_sbol, output_sam_path, exact_match=False, query_seq=None):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            # query_seq lets the caller inject a circular-extended query (origin wrap);
            # fall back to the SBOL sequence for linear targets.
            seq = query_seq if query_seq is not None else sbol_sequence(query_sbol)
            fasta_file.write(f">query_sequence\n{seq}\n")
        with open(output_sam_path, 'w') as out_sam, open(output_sam_path + '.log', 'w') as err_log:
            if(exact_match):#'bwa', 'mem', '-a', '-B', '100', '-O', '100', '-E', '100', self.index_prefix, fasta_path
                # -k 9 short seed (match blastn word_size); -D 0 keep chains
                # shorter than the best overlapping chain (nested features);
                # -W 10 is the lever that recovers the remaining short/RC parts.
                # (-c/-r/-y were tested and had no effect on this reference.)
                subprocess.run(
                    ['bwa', 'mem', '-a', '-T', '0', '-k', '9', '-D', '0', '-W', '10',
                     self.index_prefix, fasta_path],
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
            else:
                # similar mode: same maximum-recall seeding as exact (-k 9 short
                # seed, -D 0 keep nested chains, -W 10 recover short/RC parts). The
                # aligner casts the same wide net in both modes; only the downstream
                # filter differs (exact = 100% full length; similar = >=90% + NMS),
                # so near-match variants (RiboJ*, BydvJ ...) are reported and kept.
                subprocess.run(
                   ['bwa', 'mem', '-a', '-T', '0', '-k', '9', '-D', '0', '-W', '10',
                    self.index_prefix, fasta_path],
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
                
