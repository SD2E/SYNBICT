import subprocess
import tempfile
from .sbol_utils import sbol_sequence
from .Aligner import Aligner


class BlastAligner(Aligner):
    def __init__(self, index_prefix):
        super().__init__(index_prefix)

    def align(self, query_sbol, output_sam_path, exact_match=False, query_seq=None):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            # query_seq lets the caller inject a circular-extended query (origin wrap);
            # fall back to the SBOL sequence for linear targets.
            seq = query_seq if query_seq is not None else sbol_sequence(query_sbol)
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
