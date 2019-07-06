import logging
import argparse
import os
import sys

from Bio.Seq import Seq
from sbol import *
from flashtext import KeywordProcessor

def load_sbol(sbol_file):
    print('Loading ' + sbol_file + '.')

    doc = Document()
    doc.read(sbol_file)

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')

    return doc

class Feature():

    def __init__(self, nucleotides, identity, roles):
        self.nucleotides = nucleotides
        self.identity = identity
        self.roles = roles

    def reverse_complement_nucleotides(self):
        return str(Seq(self.nucleotides).reverse_complement())

class FeatureAnnotater():
    SO_REGION = 'http://identifiers.org/so/SO:0000001'
    SO_SEQUENCE_FEATURE = 'http://identifiers.org/so/SO:0000110'

    GENERIC_ROLES = {
        SO_REGION,
        SO_SEQUENCE_FEATURE
    }

    def __init__(self, feature_library, min_feature_length=40):
        self.feature_library = feature_library
        self.feature_matcher = KeywordProcessor()

        for feature in feature_library.features:
            inline_elements = ' '.join(feature.nucleotides)

            if self.__has_min_length(feature, min_feature_length):
                if inline_elements in self.feature_matcher:
                    if self.__has_non_generic_role(feature):
                        canonical_features = [cf for cf in self.feature_matcher.get_keyword(inline_elements) if self.__has_non_generic_role(cf)]

                        canonical_features.append(feature)
                else:
                    canonical_features = [feature]

                self.feature_matcher.add_keyword(inline_elements, canonical_features)

    @classmethod
    def __has_non_generic_role(cls, feature):
        return len(feature.roles.difference(cls.GENERIC_ROLES)) > 0

    # @classmethod
    # def __has_canonical_role(cls, feature, canonical_features):
    #     for canonical_feature in canonical_features:
    #         if len(feature.roles.difference(canonical_feature.roles)) == 0:
    #             return False

    #     return True

    @classmethod
    def __has_min_length(cls, feature, min_feature_length):
        return len(feature.nucleotides) >= min_feature_length

    @classmethod
    def __create_sub_component(cls, parent_definition, child_definition):
        i = 1

        while i > 0:
            try:
                sub_comp = parent_definition.components.create('_'.join([child_definition.displayId, str(i)]))
            except RuntimeError:
                sub_comp = None

            if sub_comp is None:
                i = i + 1
            else:
                sub_comp.name = child_definition.name
                sub_comp.definition = child_definition.identity

                i = -1

        return sub_comp

    @classmethod
    def __create_sequence_annotation(cls, parent_definition, child_definition, orientation, start, end,
                                     sub_comp_URI=None, parent_URI=None):
        i = 1

        while i > 0:
            try:
                seq_anno = parent_definition.sequenceAnnotations.create('_'.join([child_definition.displayId,
                                                                                  'anno',
                                                                                  str(i)]))
            except RuntimeError:
                seq_anno = None

            if seq_anno is None:
                i = i + 1
            else:
                seq_anno.name = child_definition.name
                seq_anno.description = child_definition.description
                if sub_comp_URI is not None:
                    seq_anno.component = sub_comp_URI
                if parent_URI is not None:
                    seq_anno.roles = seq_anno.roles + child_definition.roles
                    seq_anno.wasDerivedFrom = seq_anno.wasDerivedFrom + [parent_URI]
                
                location = seq_anno.locations.createRange('loc_1')

                location.orientation = orientation
                location.start = start
                location.end = end

                i = -1

        return seq_anno

    @classmethod
    def __move_component_definition(cls, source_doc, sink_doc, component_definition):
        try:
            sink_doc.addComponentDefinition(component_definition)
        except:
            pass

        for sub_comp in component_definition.components:
            try:
                sub_definition = source_doc.getComponentDefinition(sub_comp.definition)

                cls.__move_component_definition(source_doc, sink_doc, sub_definition)
            except RuntimeError:
                pass

    def __process_feature_matches(self, target_doc, target_definition, feature_matches, orientation, rc_factor=0):
        for feature_match in feature_matches:
            temp_start = feature_match[1]//2 + 1
            temp_end = (feature_match[2] + 1)//2

            if rc_factor > 0:
                start = rc_factor - temp_end
                end = rc_factor - temp_start
            else:
                start = temp_start
                end = temp_end

            for feature in feature_match[0]:
                feature_definition = self.feature_library.get_definition(feature.identity)

                sub_comp = self.__create_sub_component(target_definition, feature_definition)
                self.__create_sequence_annotation(target_definition, feature_definition, orientation, start, end,
                                                 sub_comp.identity)

                feature_doc = self.feature_library.get_document(feature.identity)

                self.__move_component_definition(feature_doc, target_doc, feature_definition)

                logging.info('Annotated %s at [%s, %s] in %s.', feature_definition.identity, start, end, target_definition.identity)

    def annotate(self, target_library, min_target_length=1000):
        for target in target_library.features:
            print('Annotating ' + target.identity)

            if self.__has_min_length(target, min_target_length):
                inline_elements = ' '.join(target.nucleotides)
                rc_elements = ' '.join(target.reverse_complement_nucleotides())

                inline_matches = self.feature_matcher.extract_keywords(inline_elements, span_info=True)
                rc_matches = self.feature_matcher.extract_keywords(rc_elements, span_info=True)

                target_doc = target_library.get_document(target.identity)
                target_definition = target_library.get_definition(target.identity)
                
                self.__process_feature_matches(target_doc, target_definition, inline_matches, SBOL_ORIENTATION_INLINE)
                self.__process_feature_matches(target_doc, target_definition, rc_matches, SBOL_ORIENTATION_REVERSE_COMPLEMENT,
                                       len(target.nucleotides) + 1)

                logging.info('Finished annotating %s.\n', target.identity)

class FeaturePruner():

    def __init__(self):
        pass

    @classmethod
    def prune(cls, target_library, cover_offset=10):
        for target in target_library.features:
            parent_definition = target_library.get_definition(target.identity)

            seq_annos = [(sa.locations.getRange().start, sa.locations.getRange().end, sa.identity, sa.component) for sa in parent_definition.sequenceAnnotations if len(sa.locations) > 0 and sa.locations[0].getTypeURI() == SBOL_RANGE]
            seq_annos.sort()

            i = 0

            while i < len(seq_annos) - 1:
                if abs(seq_annos[i + 1][0] - seq_annos[i][0]) <= cover_offset and abs(seq_annos[i + 1][1] - seq_annos[i][1]) <= cover_offset:
                    selected_index = cls.__select_feature(parent_definition, seq_annos[i], seq_annos[i + 1])

                    if seq_annos[i + 1][3] is None:
                        feature_identity2 = seq_annos[i + 1][2]
                    else:
                        feature_identity2 = parent_definition.components.get(seq_annos[i + 1][3]).definition

                    if selected_index == 0 or selected_index == 3:
                        parent_definition.sequenceAnnotations.remove(seq_annos[i + 1][2])

                        if seq_annos[i + 1][3] is not None:
                            parent_definition.components.remove(seq_annos[i + 1][3])

                        logging.info('Removed %s at [%s, %s] in %s.', feature_identity2, seq_annos[i + 1][0], seq_annos[i + 1][1], parent_definition.identity)

                        del seq_annos[i + 1]
                    else:
                        logging.info('Kept %s at [%s, %s] in %s.', feature_identity2, seq_annos[i + 1][0], seq_annos[i + 1][1], parent_definition.identity)

                    if seq_annos[i][3] is None:
                        feature_identity1 = seq_annos[i][2]
                    else:
                        feature_identity1 = parent_definition.components.get(seq_annos[i][3]).definition

                    if selected_index == 1 or selected_index == 3:
                        parent_definition.sequenceAnnotations.remove(seq_annos[i][2])

                        if seq_annos[i][3] is not None:
                            parent_definition.components.remove(seq_annos[i][3])

                        logging.info('Removed %s at [%s, %s] in %s.', feature_identity1, seq_annos[i][0], seq_annos[i][1], parent_definition.identity)

                        del seq_annos[i]
                    else:
                        logging.info('Kept %s at [%s, %s] in %s.', feature_identity1, seq_annos[i][0], seq_annos[i][1], parent_definition.identity)

                    if selected_index == 0 or selected_index == 1 or selected_index == 3:
                        i = i - 1

                i = i + 1

            logging.info('Finished pruning %s.\n', target.identity)

        kept_identities = set()

        for target in target_library.features:
            target_doc = target_library.get_document(target.identity)

            kept_identities = kept_identities.union(cls.__get_kept_identities(target_doc, target.identity))

        for doc in target_library.docs:
            identities = set()

            for component_definition in doc.componentDefinitions:
                identities.add(component_definition.identity)

            pruned_identities = identities.difference(kept_identities)

            for pruned_identity in pruned_identities:
                doc.componentDefinitions.remove(pruned_identity)

    @classmethod
    def __get_kept_identities(cls, doc, identity):
        kept_identities = {identity}

        try:
            component_definition = doc.getComponentDefinition(identity)
        except RuntimeError:
            pass

        for parent_identity in component_definition.wasDerivedFrom:
            kept_identities.add(parent_identity)

        for sub_comp in component_definition.components:
            kept_identities = kept_identities.union(cls.__get_kept_identities(doc, sub_comp.definition))

        return kept_identities

    @classmethod
    def __select_feature(cls, parent_definition, seq_anno1, seq_anno2):
        if seq_anno1[3] is None:
            feature1 = seq_anno1[2]
        else:
            feature1 = parent_definition.components.get(seq_anno1[3]).definition
            
        if seq_anno2[3] is None:
            feature2 = seq_anno2[2]
        else:
            feature2 = parent_definition.components.get(seq_anno2[3]).definition
            
        select_message = 'There appear to be redundant features {f1} at [{s1}, {e1}] and {f2} at [{s2}, {e2}] in {pd}.\nPlease select which ones to keep: 0 = {f1}, 1 = {f2}, 2 = both, 3 = neither.\n'.format(f1=feature1, s1=seq_anno1[0], e1=seq_anno1[1], f2=feature2, s2=seq_anno2[0], e2=seq_anno2[1], pd=parent_definition.identity)

        selected_index = input(select_message)

        return int(selected_index)

class FeatureLibrary():

    def __init__(self, docs, is_copy=False):
        self.features = []
        self.docs = docs
        self.__feature_map = {}

        if is_copy:
            print('Loading and copying features.')
        else:
            print('Loading features.')

        for i in range(0, len(self.docs)):
            for comp_definition in self.docs[i].componentDefinitions:
                if BIOPAX_DNA in comp_definition.types:
                    dna_seqs = []

                    for seq_URI in comp_definition.sequences:
                        try:
                            seq = self.docs[i].getSequence(seq_URI)
                        except RuntimeError:
                            seq = None

                        if seq is not None and seq.encoding == SBOL_ENCODING_IUPAC:
                            dna_seqs.append(seq)

                    if len(dna_seqs) > 0:
                        self.features.append(Feature(dna_seqs[0].elements, comp_definition.identity, set(comp_definition.roles)))

                        self.__feature_map[comp_definition.identity] = i
                    else:
                        logging.warning('DNA sequence for %s not found.', comp_definition.identity)

        if is_copy:
            for i in range(0, len(self.features)):
                comp_definition = self.get_definition(self.features[i].identity)
                feature_doc = self.get_document(self.features[i].identity)

                definition_copy = self.__copy_component_definition(comp_definition, feature_doc)

                self.__feature_map[definition_copy.identity] = self.__feature_map.pop(self.features[i].identity)
                self.features[i].identity = definition_copy.identity

    def get_document(self, identity):
        return self.docs[self.__feature_map[identity]]

    def get_definition(self, identity):
        return self.get_document(identity).getComponentDefinition(identity)

    # def remove_feature(self, identity):
    #     if identity in self.__feature_map:
    #         feature_doc = self.get_document(identity)

    #         feature_doc.componentDefinitions.remove(identity)

    #         del self.__feature_map[identity]

    @classmethod
    def __copy_component_definition(cls, comp_definition, doc):
        definition_copy = comp_definition.copy(doc)

        for seq_anno in comp_definition.sequenceAnnotations:
            if seq_anno.component is not None:
                sub_comp = comp_definition.components.get(seq_anno.component)

                sub_copy = definition_copy.components.get(sub_comp.displayId)

                anno_copy = definition_copy.sequenceAnnotations.get(seq_anno.displayId)

                anno_copy.component = sub_copy.identity

        return definition_copy

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--namespace')
    parser.add_argument('-t', '--target_files', nargs='+')
    parser.add_argument('-f', '--feature_files', nargs='*', default=[])
    parser.add_argument('-l', '--curation_log', nargs='?', default='')
    parser.add_argument('-m', '--min_target_length', nargs='?', default=1000)
    parser.add_argument('-M', '--min_feature_length', nargs='?', default=40)
    parser.add_argument('-c', '--cover_offset', nargs='?', default=10)
    parser.add_argument('-v', '--validate', action='store_true')
    # parser.add_argument('-s', '--sbh_URL', nargs='?', default=None)
    # parser.add_argument('-u', '--username', nargs='?', default=None)
    # parser.add_argument('-p', '--password', nargs='?', default=None)
    # parser.add_argument('-F', '--feature_URLs', nargs='*', default=[])
    # parser.add_argument('-T', '--target_URLs', nargs='*', default=[])
    # parser.add_argument('-o', '--sbh_output_file', nargs='?', default=None)
    
    args = parser.parse_args(args)

    if len(args.curation_log) > 0:
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(level=logging.INFO, filename=args.curation_log, filemode='w',
                            format='%(levelname)s : %(message)s')

    setHomespace(args.namespace)
    Config.setOption('validate', args.validate)
    Config.setOption('sbol_typed_uris', False)

    feature_docs = []
    for feature_file in args.feature_files:
        feature_docs.append(load_sbol(feature_file))
    feature_library = FeatureLibrary(feature_docs)

    target_documents = []
    for target_file in args.target_files:
        target_documents.append(load_sbol(target_file))
    target_library = FeatureLibrary(target_documents, True)

    feature_annotater = FeatureAnnotater(feature_library, int(args.min_feature_length))
    feature_annotater.annotate(target_library, int(args.min_target_length))

    FeaturePruner.prune(target_library, int(args.cover_offset))

    for i in range (0, len(args.target_files)):
        (target_file, file_extension) = os.path.splitext(args.target_files[i])
        target_documents[i].write('_'.join([target_file, 'annotated', str(i)]) + file_extension)

    print('Finished curating.')

if __name__ == '__main__':
    main()
