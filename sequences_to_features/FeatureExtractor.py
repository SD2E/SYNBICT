import sbol2
import subprocess
import json
import logging

# SBOL ➜ FASTA + metadata + indexing (one-time) 
class FeatureExtractor():
    def __init__(self, docs, require_sequence=True):
        #self.metadata_dict = {}
        self.fasta_records = []
        self.__extract_features(docs, require_sequence)

    def __extract_features(self, docs, require_sequence):
        for doc_index, doc in enumerate(docs):
            if(len(doc.componentDefinitions) == 0):
                logging.warning(f"No component definitions found in document {doc_index}. Skipping.")
                continue
            for comp_def in doc.componentDefinitions:
                if sbol2.BIOPAX_DNA not in comp_def.types:
                    continue
                dna_seqs = self.get_DNA_sequences(comp_def, doc)
                if require_sequence and not dna_seqs:
                    continue

                seq = dna_seqs[0].elements if dna_seqs else ''
                #new_id = str(uuid.uuid4())
                new_id = comp_def.identity

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