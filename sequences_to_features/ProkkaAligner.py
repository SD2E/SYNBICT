import os
import shutil
import subprocess
import tempfile
from .sbol_utils import sbol_sequence

DATABASE_PROTEIN_PATH = "./database_protein.fasta"
PROKKA_BIN = shutil.which("prokka") or "./prokka-1.14.6/bin/prokka"


class ProkkaAligner():
    # output_dir and database_path default to the historical fixed paths in the
    # working directory. A caller that runs Prokka concurrently (e.g. a
    # multi-threaded server) must pass a distinct output_dir per call, and the
    # protein database for that call, or the runs overwrite each other.
    def __init__(self, query_sbol, output_dir="PROKKA_SYNBICT",
                 database_path=DATABASE_PROTEIN_PATH):
        self.query_sbol = query_sbol
        self.output_sam_path = output_dir
        self.prefix = "PROKKA_SYNBICT"
        self.database_path = database_path

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
            self.database_path,
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
