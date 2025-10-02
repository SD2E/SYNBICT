import subprocess
import tempfile
import pandas as pd
from sbol_utils import sbol_sequence
from Aligner import Aligner

class VSearchAligner(Aligner):
    def __init__(self, index_prefix):
        super().__init__(index_prefix)

    def align(self, query_sbol, output_path, exact_match=False):
        with tempfile.NamedTemporaryFile(mode='w+', suffix='.fasta', delete=False) as fasta_file:
            fasta_path = fasta_file.name
            seq = sbol_sequence(query_sbol)
            fasta_file.write(f">query_sequence\n{seq}\n")
        command = [
            'vsearch', '--usearch_global', fasta_path,
            '--db', f'{self.index_prefix}.udb',
            '--id', '0.95',
            '--mincols', '40',
            '--maxaccepts', '16', '--maxhits', '1',
            '--strand', 'both',
            '--userout', output_path,
            '--userfields', "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+ql+tl+ids"
        ]
        
        command_sam = [
            'vsearch', '--usearch_global', fasta_path,
            '--db', f'{self.index_prefix}.udb',
            '--id', '0.95',
            '--mincols', '40',
            '--maxaccepts', '16', '--maxhits', '1',
            '--strand', 'both',
            '--samout', f'{output_path}.sam'
        ]
        
        with open(output_path, 'w') as out, open(output_path + '.log', 'w') as err_log:
            subprocess.run(
                command,
                stderr=err_log,
                stdout=out,
                check=True
            )
        with open(f'{output_path}.sam', 'w') as out, open(f'{output_path}.sam.log', 'w') as err_log:
            subprocess.run(
                command_sam,
                stderr=err_log,
                stdout=out,
                check=True
            )
        # write code to merge output_path and {output_path}.sam table together, for the sam table only keep the column 6 and 17
        df1 = pd.read_csv(output_path, sep="\t", header=None)
        df2 = pd.read_csv(f'{output_path}.sam', sep="\t", header=None)
        
        merged = pd.merge(df1, df2[[2, 5, 16]], left_on=1, right_on=2, how='left')
        # save this as new table, named merged.txt
        merged.to_csv(f'{output_path}.merged.txt', sep="\t", header=None, index=None)

        