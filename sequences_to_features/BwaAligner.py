import subprocess
import tempfile
from sbol_utils import load_sbol, sbol_sequence
from Aligner import Aligner


class BwaAligner(Aligner):
    def __init__(self, index_prefix):
        super().__init__(index_prefix)

    def align(self, query_sbol, output_sam_path):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            seq = sbol_sequence(query_sbol)
            fasta_file.write(f">query_sequence\n{seq}\n")
        with open(output_sam_path, 'w') as out_sam, open(output_sam_path + '.log', 'w') as err_log:
            subprocess.run(
                ['bwa', 'mem', self.index_prefix, fasta_path],
                stdout=out_sam,
                stderr=err_log,
                check=True
            )
