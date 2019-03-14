import argparse
import sys

from Bio.Seq import Seq
from sbol import *
from flashtext import KeywordProcessor

class SequenceAnnotater():
    def __init__(self, feature_files=[], sbh_URL=None, user=None, password=None, feature_URLs=[]):
        self.feature_matcher = KeywordProcessor()
        self.seq_anno_dict = {}

        # Config.setOption('sbol_typed_uris', False)

        self.load_file_features(feature_files)

        if sbh_URL is not None and user is not None and password is not None and len(feature_URLs) > 0:
            self.load_remote_features(sbh_URL, user, password, feature_URLs, part_shop)

    def load_file_features(self, feature_files):
        for feature_file in feature_files:
            doc = Document()

            doc.read(feature_file)

            self.__load_document_components(doc)
            self.__load_document_annotations(doc)
        
    def load_remote_features(self, sbh_URL, user, password, feature_URLs, part_shop):
        part_shop = PartShop(sbh_URL + '/')
        part_shop.login(user, password)

        doc = Document()

        part_shop.pull(feature_URLs, doc)

        self.__load_document_components(doc)
        self.__load_document_annotations(doc)

    def __load_document_components(self, doc):
        for root_comp in doc.componentDefinitions:
            for leaf_comp in self.__get_leaf_components(root_comp, doc):
                if leaf_comp.sequence is not None:
                    seq = leaf_comp.sequence.elements
                    rcomplement_seq = str(Seq(seq).reverse_complement())

                    self.feature_matcher.add_keyword(' '.join(seq), (leaf_comp.identity, SBOL_ORIENTATION_INLINE))
                    self.feature_matcher.add_keyword(' '.join(rcomplement_seq), (leaf_comp.identity,
                                                     SBOL_ORIENTATION_REVERSE_COMPLEMENT))

    def __get_leaf_components(self, root_comp, doc):
        leaf_comps = []

        if len(root_comp.components) == 0:
            leaf_comps.append(root_comp)
        else:
            for sub_comp in root_comp.components:
                try:
                    leaf_comp = doc.getComponentDefinition(sub_comp.definition)

                    leaf_comps.extend(self.__get_leaf_components(leaf_comp, doc))
                except RuntimeError:
                    pass

        return leaf_comps

    def __load_document_annotations(self, doc):
        for comp in doc.componentDefinitions:
            for seq_anno in comp.sequenceAnnotations:
                if (seq_anno.component is None and len(seq_anno.locations) == 1
                    and seq_anno.locations[0].getTypeURI() == 'http://sbols.org/v2#Range'
                    and comp.sequence is not None):
                    seq_range = seq_anno.locations.getRange()

                    if seq_range.orientation == SBOL_ORIENTATION_INLINE:
                        seq = comp.sequence.elements[seq_range.start - 1:seq_range.end]
                        rcomplement_seq = str(Seq(seq).reverse_complement())
                    else:
                        rcomplement_seq = comp.sequence.elements[seq_range.start - 1:seq_range.end]
                        seq = str(Seq(seq).reverse_complement())
                        
                    self.feature_matcher.add_keyword(' '.join(seq), (comp.identity, SBOL_ORIENTATION_INLINE,
                                                     seq_anno.identity, seq_anno.name, seq_anno.description,
                                                     seq_anno.roles))
                    self.feature_matcher.add_keyword(' '.join(rcomplement_seq), (comp.identity,
                                                     SBOL_ORIENTATION_REVERSE_COMPLEMENT, seq_anno.identity,
                                                     seq_anno.name, seq_anno.description, seq_anno.roles))

                    self.seq_anno_dict[seq_anno.identity] = seq_anno

    def annotate_sequences(self, sequence_files=[], sbh_URL=None, user=None, password=None, sequence_URLs=[], 
                           sbh_output=None, validate=False):
        Config.setOption('validate', validate)

        docs = []

        for sequence_file in sequence_files:
            docs.append(Document())

            docs[-1].read(sequence_file)

        if sbh_URL is not None and user is not None and password is not None and len(sequence_URLs) > 0:
            part_shop = PartShop(sbh_URL + '/')
            part_shop.login(user, password)

            docs.append(Document())

            part_shop.pull(sequence_URLs, docs[-1])

        for doc in docs:
            for comp in doc.componentDefinitions:
                if comp.sequence is not None:
                    unannotated_seq = ' '.join(comp.sequence.elements)

                    feature_matches = self.feature_matcher.extract_keywords(unannotated_seq, span_info=True)

                    for i in range(0, len(feature_matches)):
                        annotated = False
                        j = i

                        if len(feature_matches[i][0]) == 2:
                            while not annotated:
                                try:
                                    sub_comp = comp.components.create('comp' + str(j))

                                    sub_comp.definition = feature_matches[i][0][0]

                                    annotated = True
                                except RuntimeError:
                                    j = j + 1

                            annotated = False
                            j = i

                            while not annotated:
                                try:
                                    seq_anno = comp.sequenceAnnotations.create('anno' + str(j))

                                    seq_anno.component = sub_comp.identity
                                    
                                    location = seq_anno.locations.createRange('loc1')

                                    location.orientation = feature_matches[i][0][1]
                                    location.start = feature_matches[i][1]//2 + 1
                                    location.end = (feature_matches[i][2] + 1)//2

                                    annotated = True
                                except RuntimeError:
                                    j = j + 1
                        elif len(feature_matches[i][0]) == 6:
                            while not annotated:
                                try:
                                    seq_anno = comp.sequenceAnnotations.create('anno' + str(j))

                                    seq_anno.name = feature_matches[i][0][3]
                                    seq_anno.description = feature_matches[i][0][4]
                                    seq_anno.roles = seq_anno.roles + feature_matches[i][0][5]
                                    seq_anno.wasDerivedFrom = seq_anno.wasDerivedFrom + [feature_matches[i][0][2]]
                                    
                                    location = seq_anno.locations.createRange('loc1')

                                    location.orientation = feature_matches[i][0][1]
                                    location.start = feature_matches[i][1]//2 + 1
                                    location.end = (feature_matches[i][2] + 1)//2

                                    annotated = True
                                except RuntimeError:
                                    j = j + 1

        for i in range(0, len(sequence_files)):
            sequence_filename = sequence_files[i][:sequence_files[i].index('.xml')]

            docs[i].write(sequence_filename + '_annotated.xml')

        if len(docs) == len(sequence_files) + 1 and sbh_output is not None:
            docs[-1].write(sbh_output + '.xml')


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sequence_files', nargs='*', default=[])
    parser.add_argument('-S', '--sequence_URLs', nargs='*', default=[])
    parser.add_argument('-f', '--feature_files', nargs='*', default=[])
    parser.add_argument('-F', '--feature_URLs', nargs='*', default=[])
    parser.add_argument('-u', '--sbh_URL', nargs='?', default=None)
    parser.add_argument('-n', '--username', nargs='?', default=None)
    parser.add_argument('-p', '--password', nargs='?', default=None)
    parser.add_argument('-o', '--sbh_output', nargs='?', default='sbh_annotated')
    
    args = parser.parse_args(args)

    sequence_annotater = SequenceAnnotater(args.feature_files, args.sbh_URL, args.username, args.password,
                                           args.feature_URLs)

    sequence_annotater.annotate_sequences(args.sequence_files, args.sbh_URL, args.username, args.password,
                                          args.sequence_URLs, args.sbh_output)

    # def extract_sequence(sub_def, seq_anno, comp_def, doc):
    #     if len(comp_def.sequences) == 1 and len(seq_anno.locations) == 1:
    #         comp_seq = doc.sequences.get(comp_def.sequences[0])

    #         rang = seq_anno.locations.getRange()
    #         if rang.start > len(comp_seq.elements):
    #             return False
    #         else:
    #             elements = comp_seq.elements[rang.start - 1:rang.end]

    #         if rang.orientation == 'http://sbols.org/v2#reverseComplement':
    #             elements = str(Seq(elements).reverse_complement())

    #         sub_seq = Sequence(sub_def.displayId + '_seq', elements, SBOL_ENCODING_IUPAC, '1')
            
    #         doc.addSequence(sub_seq)

    #         sub_def.sequences = sub_def.sequences + [sub_seq.identity]
            
    #         return True
    #     else:
    #         return False
        
    # def reconcile_roles(canonical_def, comp_def, generic_roles):
    #     same_roles = True

    #     for role in comp_def.roles:
    #         if role not in canonical_def.roles:
    #             same_roles = False

    #     if same_roles:
    #         for role in canonical_def.roles:
    #             if role not in comp_def.roles:
    #                 same_roles = False

    #     if not same_roles:
    #         canonical_roles = []

    #         for role in comp_def.roles:
    #             if role not in generic_roles:
    #                 canonical_roles.append(role)

    #         if len(canonical_roles) > 0:
    #             canonical_def.roles = canonical_roles

    # setHomespace('http://hub.sd2e.org/')
    # Config.setOption('sbol_typed_uris', False)
    # Config.setOption('validate', True)
    # # Config.setOption('serialization_format', 'rdfxml')
    # Config.setOption('version', '1')

    # doc = Document()

    # part_shop = PartShop('https://hub.sd2e.org/')
    # part_shop.login('sd2_service@sd2e.org', 'jWJ1yztJl2f7RaePHMtXmxBBHwNt')

    # # part_shop.pull('https://hub.sd2e.org/user/sd2e/design/yeast_gates_plasmids/1', doc)
    # for i in range(1, 10):
    #     part_shop.pull('https://hub.sd2e.org/user/sd2e/design/YG_plasmid_00' + str(i) + '/1', doc)
    #     print('load ' + str(i))
    # for i in range(10, 26):
    #     part_shop.pull('https://hub.sd2e.org/user/sd2e/design/YG_plasmid_0' + str(i) + '/1', doc)
    #     print('load ' + str(i))
        
    # setHomespace('http://hub.sd2e.org/user/sd2e/design')

    # non_canonical_ids = set()
    # non_canonical_names = set()
    # name_to_def = {}
    # generic_roles = {
    #     SO_GENE, 
    #     'http://identifiers.org/so/SO:0000001', 
    #     'http://identifiers.org/so/SO:0000110', 
    #     'http://identifiers.org/so/SO:0000804'
    # }

    # for comp_def in doc.componentDefinitions:
    #     for seq_anno in comp_def.sequenceAnnotations:
    #         if seq_anno.component is None:
    #             print('SequenceAnnotation ' + seq_anno.identity + ' has no associated Component.')
    #         else:
    #             comp = comp_def.components.get(seq_anno.component)
    #             sub_def = doc.componentDefinitions.get(comp.definition)
                
    #             if len(sub_def.sequences) > 0 or extract_sequence(sub_def, seq_anno, comp_def, doc):
    #                 sub_seq = doc.sequences.get(sub_def.sequences[0])
                    
    #                 if len(sub_seq.elements) > 1:
    #                     if sub_def.name in name_to_def:
    #                         canonical_def = name_to_def[sub_def.name]
    #                         canonical_seq = doc.sequences.get(canonical_def.sequences[0])

    #                         if canonical_seq.elements.lower() == sub_seq.elements.lower():
    #                             reconcile_roles(canonical_def, sub_def, generic_roles)
    #                         else:
    #                             non_canonical_ids.add(sub_def.identity)
    #                             non_canonical_names.add(sub_def.name)
    #                     else:
    #                         name_to_def[sub_def.name] = sub_def
    #             else:
    #                 non_canonical_ids.add(sub_def.identity)
    #                 non_canonical_names.add(sub_def.name)
                    
    # for comp_def in doc.componentDefinitions:  
    #     deleted_comp_ids = set()
        
    #     for comp in comp_def.components:
    #         if comp.definition not in non_canonical_ids:
    #             sub_def = doc.componentDefinitions.get(comp.definition)
                
    #             if sub_def.name in name_to_def:
    #                 canonical_def = name_to_def[sub_def.name]
                    
    #                 if sub_def.identity != canonical_def.identity:
    #                     comp.definition = canonical_def.identity
                        
    #                     for sub_seq_id in sub_def.sequences:
    #                         doc.sequences.remove(sub_seq_id)

    #                     doc.componentDefinitions.remove(sub_def.identity)
    #             else:
    #                 deleted_comp_ids.add(comp.identity)
                    
    #                 for sub_seq_id in sub_def.sequences:
    #                     doc.sequences.remove(sub_seq_id)

    #                 doc.componentDefinitions.remove(sub_def.identity)
                
    #     for deleted_comp_id in deleted_comp_ids:
    #         comp_def.components.remove(deleted_comp_id)
            
    #     for anno in comp_def.sequenceAnnotations:
    #         if anno.component in deleted_comp_ids:
    #             anno.component = None
              
    # print('writing')
    # doc.write('bograth.xml')
    # print(non_canonical_names)

if __name__ == '__main__':
    main()