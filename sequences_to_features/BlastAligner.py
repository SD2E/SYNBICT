import subprocess
import tempfile
from sbol_utils import sbol_sequence
from Aligner import Aligner


class BlastAligner(Aligner):
    def __init__(self, index_prefix):
        super().__init__(index_prefix)

    def align(self, query_sbol, output_sam_path, exact_match=False):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            seq = sbol_sequence(query_sbol)
            fasta_file.write(f">query_sequence\n{seq}\n")
        output_path = output_sam_path
        blast_command_simi = [
            'blastn', 
            '-query', 
            fasta_path,
            '-db',
            self.index_prefix,
            '-out',
            output_path,
            '-outfmt',
            '6 std qlen slen nident'
        ]
        with open(output_path, 'w') as out, open(output_path + '.log', 'w') as err_log:
            subprocess.run(
                blast_command_simi,
                stdout=out,
                stderr=err_log,
                check=True
            )
