import subprocess
import tempfile, os, glob
from sbol_utils import sbol_sequence
from Aligner import Aligner


class MmseqAligner(Aligner):
    def __init__(self, index_prefix):
        super().__init__(index_prefix)

    def align(self, query_sbol, output_sam_path, exact_match=False):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            seq = sbol_sequence(query_sbol)
            fasta_file.write(f">query_sequence\n{seq}\n")
        #output_path = 'output.txt'
        # command = [
        #     'mmseqs', 'createdb', fasta_path, 'query_db', '&&',
        #     'mmseqs', 'search', 'query_db', f'{self.index_prefix}.mmseqs', 'result_db', 'mmseqs_tmp',
        #     '--search-type', '3', '--min-seq-id', '1.0', '--alignment-mode', '3', '-s', '7.5', '&&',
        #     'mmseqs', 'convertalis', 'query_db', f'{self.index_prefix}.mmseqs', 'result_db', output_sam_path
        # ]
        # reverse query and db can get more exact matches, but it is slower
        exact_command = [
            'vsearch', '--usearch_global', fasta_path, '--db', f'{self.index_prefix}.udb',
            '--id', '1.0', '--mincols', '40', '--strand', 'both',
            '--maxaccepts', '0', '--maxhits', '0', '--blast6out', output_sam_path
        ]
        # Clean up previous result files
        for f in glob.glob("result_db*"):
            if os.path.isfile(f):
                os.remove(f)
                
        with open(output_sam_path, 'w') as out:
            if(exact_match):
                exact_command
            # Step 1: createdb
            subprocess.run(['mmseqs', 'createdb', fasta_path, 'query_db'], check=True)

            # Step 2: search
            subprocess.run([
                'mmseqs', 'search', 'query_db', f'{self.index_prefix}.mmseqs', 'result_db', 'mmseqs_tmp',
                '--search-type', '3', '--min-seq-id', '1.0', '--alignment-mode', '3', '-s', '7.5'
            ], check=True)

            # Step 3: convertalis
            subprocess.run([
                'mmseqs', 'convertalis', 'query_db', f'{self.index_prefix}.mmseqs', 'result_db', output_sam_path
            ], stdout=out, check=True)

        
        # filter aligned.txt
        # Now filter results in Python
        # if exact_match:
        #     # Filter for exact matches (100% identity)
        #     with open(output_path) as infile, open(output_sam_path, 'w') as outfile:
        #         for line in infile:
        #             fields = line.strip().split('\t')
        #             if float(fields[2]) == 100.0:  # identity column is the 3rd field (0-based index 2)
        #                 outfile.write(line)
        # else:
        #     # Filter for matches with at least 90% identity
        #     with open(output_path) as infile, open(output_sam_path, 'w') as outfile:
        #         for line in infile:
        #             fields = line.strip().split('\t')
        #             if float(fields[2]) >= 90.0:  # identity column is the 3rd field (0-based index 2)
        #                 outfile.write(line)
