import subprocess
import tempfile
from sbol_utils import sbol_sequence
from Aligner import Aligner
class BowtieAligner(Aligner):
    def __init__(self, index_prefix):
        super().__init__(index_prefix)

    def align(self, query_sbol, output_sam_path, exact_match=False):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            print(f"Writing query sequence to {fasta_path}")
            seq = sbol_sequence(query_sbol)
            fasta_file.write(f">query_sequence\n{seq}\n")
        bowtie2_command = [
            "bowtie2",
            "-x", "test",
            "-f",
            "-U", fasta_path,
            "--local",
            "-a",
            "-S", output_sam_path
        ]
        bowtie2_command_simi = [
            "bowtie2",
            "-x", "test",
            "-f",
            "-U", fasta_path,
            "-a",
            "--very-sensitive-local",
            "-S", output_sam_path
        ]
        with open(output_sam_path, 'w') as out_sam, open(output_sam_path + '.log', 'w') as err_log:
            if(exact_match): # For long reads (ONT or PacBio)
                subprocess.run(
                    bowtie2_command,
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
            else:
                subprocess.run(
                    bowtie2_command_simi,
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
                
