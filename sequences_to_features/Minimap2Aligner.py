import subprocess
import tempfile
from .sbol_utils import sbol_sequence
from .Aligner import Aligner
class Minimap2Aligner(Aligner):
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
            if(exact_match): # For long reads (ONT or PacBio)
                subprocess.run(
                    ['minimap2', '-ax', 'map-ont', '-N', '1', '--secondary=no', '-p', '0.9', f'{self.index_prefix}.mmi', fasta_path],
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
            else:
                subprocess.run(
                    ['minimap2', '-ax', 'map-ont', f'{self.index_prefix}.mmi', fasta_path],
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
                
