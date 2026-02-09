import sbol2
import subprocess
import json
import logging
from collections.abc import Mapping, Iterable
from Bio import SeqIO
from Bio.Seq import Seq

# SBOL ➜ FASTA + metadata + indexing (one-time) 
class FeatureExtractor():
    def __init__(self, docs, require_sequence=True):
        #self.metadata_dict = {}
        self.fasta_records = []
        self.protein_fasta_records = []
        self.cds_id_map = {}
        self.id_cds_map = {}
        self.__extract_features(docs, require_sequence)
        self.__extract_protein_sequences(docs, require_sequence)

    def __extract_features(self, docs, require_sequence):
        """
        Accepts:
        - dict-like: {label: doc}
        - list/tuple of docs
        - single doc object
        Appends (id, seq) tuples to self.fasta_records.
        """

        # --- normalize to an iterator of (label, doc) ---
        if isinstance(docs, Mapping):
            pairs = docs.items()                               # (key, doc)
        elif isinstance(docs, Iterable) and not isinstance(docs, (str, bytes)):
            pairs = enumerate(docs)                            # (index, doc)
        else:
            pairs = [(0, docs)]                                # single doc

        for label, doc in pairs:
            comp_defs = getattr(doc, "componentDefinitions", None)
            if not comp_defs:
                logging.warning(f"No component definitions found in document {label}. Skipping.")
                continue

            for comp_def in comp_defs:
                if sbol2.BIOPAX_DNA not in getattr(comp_def, "types", []):
                    continue

                dna_seqs = self.get_DNA_sequences(comp_def, doc)
                if require_sequence and not dna_seqs:
                    continue

                seq = dna_seqs[0].elements if dna_seqs else ""
                new_id = getattr(comp_def, "identity", None)

                # Store plain (ID, sequence) tuple instead of SeqRecord
                self.fasta_records.append((new_id, seq))
                # Store metadata separately
                # self.metadata_dict[new_id] = {
                #     'original_identity': comp_def.identity,
                #     'name': comp_def.name,
                #     'displayId': comp_def.displayId,
                #     'roles': comp_def.roles,
                #     'wasDerivedFrom': comp_def.wasDerivedFrom,
                #     'doc_index': doc_index
                # }

    def write_fasta(self, fasta_path):
        with open(fasta_path, "w") as fasta_file:
            for record_id, sequence in self.fasta_records:
                fasta_file.write(f">{record_id}\n{sequence}\n")

    def write_metadata(self, metadata_path):
        with open(metadata_path, 'w') as f:
            json.dump(self.metadata_dict, f, indent=2)


    def build_index(self, fasta_path, index_prefix, tool='bwa'):
        if tool == 'bwa':
            subprocess.run(['bwa', 'index', '-p', index_prefix, fasta_path], check=True)
        elif tool == 'minimap2':
            subprocess.run(['minimap2', '-d', f'{index_prefix}.mmi', fasta_path], check=True)
        elif tool == 'bowtie2':
            subprocess.run(['bowtie2-build', fasta_path, index_prefix], check=True)
        elif tool == 'blast':
            subprocess.run(['makeblastdb', '-in', fasta_path, '-dbtype', 'nucl', '-out', index_prefix], check=True)
        elif tool == 'vsearch':
            subprocess.run(['vsearch', '--makeudb_usearch', fasta_path, '-output', f'{index_prefix}.udb'], check=True)
        elif tool == 'mmseqs2':
            subprocess.run(['mmseqs', 'createdb', fasta_path, f'{index_prefix}.mmseqs'], check=True)
        else:
            raise ValueError("Unsupported tool for indexing")

    @classmethod
    def get_DNA_sequences(cls, comp_definition, doc):
        dna_seqs = []
        for seq_URI in comp_definition.sequences:
            try:
                seq = doc.getSequence(seq_URI)
            except (RuntimeError, sbol2.Document.NotFoundError):
                seq = None

            if seq and seq.encoding == sbol2.SBOL_ENCODING_IUPAC:
                dna_seqs.append(seq)

        return dna_seqs
    
    def __extract_protein_sequences(self, docs, require_sequence):

        if isinstance(docs, Mapping):
            pairs = docs.items()  
        elif isinstance(docs, Iterable) and not isinstance(docs, (str, bytes)):
            pairs = enumerate(docs)     
        else:
            pairs = [(0, docs)] 
        counter = 1
        for label, doc in pairs:
            comp_defs = getattr(doc, "componentDefinitions", None)
            
            if not comp_defs:
                print(f"No component definitions found in document {label}. Skipping.")
                continue

            for comp_def in comp_defs:
                if sbol2.BIOPAX_DNA not in getattr(comp_def, "types", []):
                    continue
                roles = comp_def.roles[0]
                
                if(roles == "http://identifiers.org/so/SO:0000316"):
                    dna_seqs = self.get_DNA_sequences(comp_def, doc)
                    if require_sequence and not dna_seqs:
                        continue

                    seq = dna_seqs[0].elements if dna_seqs else ""
                    protein_seq = Seq(seq).translate(to_stop=True)
                    
                    id = getattr(comp_def, "identity")
                    cds_id = f"CDS_{counter:06d}"
                    self.cds_id_map[cds_id] = id
                    self.id_cds_map[id] = cds_id
                    counter += 1
                    old_des = getattr(comp_def, "description", None)
                    if old_des is None:
                        old_des = cds_id
                    total_id = f"{cds_id} {old_des}"

                    # Store plain (ID, sequence) tuple instead of SeqRecord
                    self.protein_fasta_records.append((total_id, protein_seq)) # this is DNA seq
                    
    def write_protein_fasta(self, fasta_path):
        with open(fasta_path, "w") as fasta_file:
            for record_id, sequence in self.protein_fasta_records:
                fasta_file.write(f">{record_id} \n{sequence}\n")