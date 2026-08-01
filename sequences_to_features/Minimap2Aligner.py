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
                # Maximum-recall settings. Emit many secondary alignments (-N 200,
                # -p 0.01) so nested/duplicated short features are reported, and
                # index the reference on the fly with a short seed (-k 9 -w 3) for
                # short features. Uses the fasta (not the .mmi) so -k/-w take effect.
                # FP stays 0 under the downstream exact full-length filter.
                subprocess.run(
                    ['minimap2', '-ax', 'map-ont', '-N', '200', '-p', '0.01',
                     '-k', '9', '-w', '3', f'{self.index_prefix}.fasta', fasta_path],
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
            else:
                # similar mode: same maximum-recall settings as exact — emit many
                # secondary alignments (-N 200, -p 0.01) and index on the fly with a
                # short seed (-k 9 -w 3, from the fasta so -k/-w take effect) so
                # short and near-match variant parts are reported. The downstream
                # >=90% + NMS filter (not a full-length exact check) keeps the good
                # ones. Using the default .mmi here would drop short/variant parts.
                subprocess.run(
                    ['minimap2', '-ax', 'map-ont', '-N', '200', '-p', '0.01',
                     '-k', '9', '-w', '3', f'{self.index_prefix}.fasta', fasta_path],
                    stdout=out_sam,
                    stderr=err_log,
                    check=True
                )
                
