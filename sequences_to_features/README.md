cd /home/sophia/git_repo/SYNBICT

## First step, building index
python -m sequences_to_features.sequences_to_features -n http://mynamespace.org -f /home/sophia/git_repo/SYNBICT/example/jet_libs/CIDAR_MoClo_*.xml -bi

## Second step, run alignment.
### Similar Match 
python -m sequences_to_features.sequences_to_features -n http://mynamespace.org -f /home/sophia/git_repo/SYNBICT/example/jet_libs/CIDAR_MoClo_*.xml -t /home/sophia/git_repo/SYNBICT/11508_addgene_out.xml -bwa -np -o 11508_out.xml
### Exact Match
python -m sequences_to_features.sequences_to_features -n http://mynamespace.org -f /home/sophia/git_repo/SYNBICT/example/jet_libs/CIDAR_MoClo_*.xml -t /home/sophia/git_repo/SYNBICT/11508_addgene_out.xml -bwa -np -exact 
### Protein Match
python -m sequences_to_features.sequences_to_features -n http://mynamespace.org -f /home/sophia/git_repo/SYNBICT/example/jet_libs/CIDAR_MoClo_*.xml -t /home/sophia/git_repo/SYNBICT/11508_addgene_out.xml -bwa -np -prokka -o 11508_out_protein.xml
