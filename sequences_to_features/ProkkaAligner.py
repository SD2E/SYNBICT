import os
import shutil
import subprocess
import tempfile
from .sbol_utils import sbol_sequence

DATABASE_PROTEIN_PATH = "./database_protein.fasta"
PROKKA_BIN = "./prokka-1.14.6/bin/prokka"


class ProkkaAligner():
    def __init__(self, query_sbol):
        self.query_sbol = query_sbol
        self.output_sam_path = "PROKKA_SYNBICT"
        self.prefix = "PROKKA_SYNBICT"

    # query_sbol: SBOLDocument, same as BwaAligner
    # output_sam_path: str, prokka output path
    # exact_match: bool, don't need this parameter
    def align(self): 
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            seq = sbol_sequence(self.query_sbol)
            if not seq or not str(seq).strip():
                raise ValueError("sbol_sequence(query_sbol) returned an empty sequence")
            fasta_file.write(f">query_sequence\n{seq}\n")
        output_path = self.output_sam_path
        prefix = self.prefix
        if os.path.isdir(output_path):
            shutil.rmtree(output_path)
        prokka_cmd = [
            PROKKA_BIN, 
            '--debug', 
            '--quiet',
            '--rfam', 
            '--proteins', 
            DATABASE_PROTEIN_PATH,
            '--outdir', output_path, 
            '--prefix', prefix,
            '--locustag', 'LOCUS',
            fasta_path
        ]
        try:
            subprocess.run(
                prokka_cmd,
                check=False,
            )
        finally:
            # Clean up temp fasta
            try:
                os.remove(fasta_path)
            except OSError:
                pass