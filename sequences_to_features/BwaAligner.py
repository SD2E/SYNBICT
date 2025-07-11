import subprocess
import tempfile
from sbol_utils import sbol_sequence
from Aligner import Aligner


class BwaAligner(Aligner):
    def __init__(self, index_prefix):
        super().__init__(index_prefix)

    def align(self, query_sbol, output_sam_path, exact_match=False):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            seq = sbol_sequence(query_sbol)
            fasta_file.write(f">query_sequence\n{seq}\n") #'bwa', 'mem', '-B', '100', '-O', '100', '-E', '100',
        with open(output_sam_path, 'w') as out_sam, open(output_sam_path + '.log', 'w') as err_log:
            if(exact_match):
                subprocess.run(
                    ['bwa', 'mem', '-a', '-B', '100', '-O', '100', '-E', '100', self.index_prefix, fasta_path],
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
            else:
                subprocess.run(
                    ['bwa', 'mem', '-a', '-T', '0', self.index_prefix, fasta_path],
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
                
