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

        self.encode_file_features(feature_files)

        if sbh_URL is not None and user is not None and password is not None and len(feature_URLs) > 0:
            self.encode_remote_features(sbh_URL, user, password, feature_URLs, part_shop)

    @classmethod
    def __load_sbol_files(cls, sbol_files):
        docs = []

        for sbol_file in sbol_files:
            docs.append(Document())

            docs[-1].read(sbol_file)

        return docs

    @classmethod
    def __load_remote_sbol(cls, sbh_URL, user, password, sbol_URLs, part_shop):
        doc = Document()

        part_shop = PartShop(sbh_URL)
        part_shop.login(user, password)

        part_shop.pull(sbol_URLs, doc)

        return doc

    def encode_file_features(self, feature_files):
        docs = self.__load_sbol_files(feature_files)

        for doc in docs:
            self.__encode_document_components(doc)
            self.__encode_document_annotations(doc)

    def encode_remote_features(self, sbh_URL, user, password, feature_URLs, part_shop):
        doc = self.__load_remote(sbh_URL, user, password, feature_URLs, part_shop)

        self.__encode_document_components(doc)
        self.__encode_document_annotations(doc)

    def __get_leaf_components(self, root_comp, doc):
        leaf_comps = []

        if len(root_comp.components) == 0:
            leaf_comps.append(root_comp)
        else:
            for sub_comp in root_comp.components:
                try:
                    leaf_comp = doc.getComponentDefinition(sub_comp.definition)

                except RuntimeError:
                    leaf_comp = None

                if leaf_comp is not None:
                    leaf_comps.extend(self.__get_leaf_components(leaf_comp, doc))

        return leaf_comps

    def __encode_document_components(self, doc):
        for root_comp in doc.componentDefinitions:
            for leaf_comp in self.__get_leaf_components(root_comp, doc):
                if leaf_comp.sequence is not None:
                    seq = leaf_comp.sequence.elements
                    rcomplement_seq = str(Seq(seq).reverse_complement())

                    self.feature_matcher.add_keyword(' '.join(seq), (leaf_comp.identity, leaf_comp.displayId,
                                                     SBOL_ORIENTATION_INLINE))
                    self.feature_matcher.add_keyword(' '.join(rcomplement_seq), (leaf_comp.identity,
                                                     leaf_comp.displayId, SBOL_ORIENTATION_REVERSE_COMPLEMENT))

    def __encode_document_annotations(self, doc):
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
                        
                    self.feature_matcher.add_keyword(' '.join(seq), (seq_anno.displayId, seq_anno.name,
                                                     seq_anno.description, seq_anno.roles, seq_anno.identity,
                                                     SBOL_ORIENTATION_INLINE))
                    self.feature_matcher.add_keyword(' '.join(rcomplement_seq), (seq_anno.displayId, seq_anno.name,
                                                     seq_anno.description, seq_anno.roles, seq_anno.identity,
                                                     SBOL_ORIENTATION_REVERSE_COMPLEMENT))

                    self.seq_anno_dict[seq_anno.identity] = seq_anno

    @classmethod
    def __create_sub_component(cls, feature_URI, feature_ID):
        i = 1

        while not i < 0:
            try:
                sub_comp = comp.components.create('_'.join([feature_ID, str(i)]))
            except RuntimeError:
                sub_comp = None

            if sub_comp is None:
                i = i + 1
            else:
                sub_comp.definition = feature_URI

                i = -1

        return sub_comp

    @classmethod
    def __create_sub_component_annotation(cls, sub_comp, feature_URI, orientation, start, end):
        i = 1

        while not i < 0:
            try:
                seq_anno = comp.sequenceAnnotations.create('_'.join([sub_comp.displayId, 'anno', str(i)]))
            except RuntimeError:
                seq_anno = None

            if seq_anno is None:
                i = i + 1
            else:
                seq_anno.component = sub_comp.identity
                seq_anno.wasDerivedFrom = seq_anno.wasDerivedFrom + [feature_URI]
                
                location = seq_anno.locations.createRange('loc_1')

                location.orientation = orientation
                location.start = start
                location.end = end

                i = -1

        return seq_anno

    @classmethod
    def __create_sequence_annotation(cls, anno_ID, name, description, roles, parent_URI, orientation, start, end):
        i = 1

        while not i < 0:
            try:
                seq_anno = comp.sequenceAnnotations.create('_'.join([anno_ID, str(i)]))
            except RuntimeError:
                seq_anno = None

            if seq_anno is None:
                i = i + 1
            else:
                seq_anno.name = name
                seq_anno.description = description
                seq_anno.roles = seq_anno.roles + roles
                seq_anno.wasDerivedFrom = seq_anno.wasDerivedFrom + [parent_URI]
                
                location = seq_anno.locations.createRange('loc1')

                location.orientation = orientation
                location.start = start
                location.end = end

                i = -1

        return seq_anno

    def annotate_sequences(self, sequence_files=[], sbh_URL=None, user=None, password=None, sequence_URLs=[], 
                           output_file=None, validate=False):
        Config.setOption('validate', validate)

        docs = self.__load_sbol_files(sequence_files)

        if sbh_URL is not None and user is not None and password is not None and len(sequence_URLs) > 0:
            docs.append(self.__load_remote_sbol(sequence_URLs))

        for doc in docs:
            for comp in doc.componentDefinitions:
                if comp.sequence is not None:
                    unannotated_seq = ' '.join(comp.sequence.elements)

                    feature_matches = self.feature_matcher.extract_keywords(unannotated_seq, span_info=True)

                    for feature_match in feature_matches:
                        if len(feature_match[0]) == 2:
                            sub_comp = self.__create_sub_component(feature_match[0][0], feature_match[0][1])

                            self.__create_sub_component_annotation(sub_comp, feature_match[0][0], feature_match[0][2],
                                                                   feature_match[1]//2 + 1, (feature_match[2] + 1)//2)
                        elif len(feature_match[0]) == 6:
                            self.__create_sequence_annotation(feature_match[0][0], feature_match[0][1],
                                                              feature_match[0][2], feature_match[0][3],
                                                              feature_match[0][4], feature_match[0][5],
                                                              feature_match[1]//2 + 1, (feature_match[2] + 1)//2)

        for i in range(0, len(sequence_files)):
            sequence_filename = sequence_files[i][:sequence_files[i].index('.xml')]

            docs[i].write(sequence_filename + '_annotated.xml')

        if len(docs) == len(sequence_files) + 1 and output_file is not None:
            docs[-1].write(output_file + '.xml')


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
    parser.add_argument('-o', '--output_file', nargs='?', default='sbh_annotated')
    
    args = parser.parse_args(args)

    sequence_annotater = SequenceAnnotater(args.feature_files, args.sbh_URL, args.username, args.password,
                                           args.feature_URLs)

    sequence_annotater.annotate_sequences(args.sequence_files, args.sbh_URL, args.username, args.password,
                                          args.sequence_URLs, args.output_file)

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