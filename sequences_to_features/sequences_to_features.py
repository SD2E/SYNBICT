import logging
import argparse
import os
import sys

from Bio.Seq import Seq
from sbol import *
from flashtext import KeywordProcessor

def load_sbol(sbol_file):
    print('Loading ' + sbol_file + '...')

    doc = Document()
    doc.read(sbol_file)

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')

    return doc

class Feature():

    SO_REGION = 'http://identifiers.org/so/SO:0000001'
    SO_SEQUENCE_FEATURE = 'http://identifiers.org/so/SO:0000110'

    GENERIC_ROLES = {
        SO_REGION,
        SO_SEQUENCE_FEATURE
    }

    def __init__(self, nucleotides, identity, roles, parent_identities=[]):
        self.nucleotides = nucleotides
        self.identity = identity
        self.parent_identities = parent_identities
        self.roles = roles

    def reverse_complement_nucleotides(self):
        return str(Seq(self.nucleotides).reverse_complement())

    @classmethod
    def has_non_generic_role(cls, roles):
        return len(roles.difference(cls.GENERIC_ROLES)) > 0

    @classmethod
    def has_same_roles(cls, roles1, roles2):
        return len(roles1.difference(roles2)) == 0

    def is_non_generic(self):
        return self.has_non_generic_role(self.roles)

    def is_same_as(self, other_feature):
        return self.has_same_roles(self.roles, other_feature.roles)

class FeatureLibrary():

    def __init__(self, docs, is_copy=False):
        self.features = []
        self.docs = docs
        self.__feature_map = {}

        if is_copy:
            print('Loading and copying features...\n')
        else:
            print('Loading features...\n')

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
                        self.features.append(Feature(dna_seqs[0].elements, comp_definition.identity, set(comp_definition.roles),
                            comp_definition.wasDerivedFrom))

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

    def get_features(self, min_feature_length=0):
        features = []

        parent_identities = set()

        for feature in self.features:
            for parent_identity in feature.parent_identities:
                parent_identities.add(parent_identity)

        for feature in self.features:
            if len(feature.nucleotides) > min_feature_length and feature.identity not in parent_identities:
                features.append(feature)

        return features

    def get_document(self, identity):
        return self.docs[self.__feature_map[identity]]

    def get_definition(self, identity):
        return self.get_document(identity).getComponentDefinition(identity)

    def has_feature(self, identity):
        return identity in self.__feature_map

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

class FeatureAnnotater():

    def __init__(self, feature_library, min_feature_length):
        self.feature_library = feature_library
        self.feature_matcher = KeywordProcessor()

        for feature in feature_library.features:
            inline_elements = ' '.join(feature.nucleotides)

            if self.__has_min_length(feature, min_feature_length):
                if inline_elements in self.feature_matcher:
                    if feature.is_non_generic():
                        canonical_features = [cf for cf in self.feature_matcher.get_keyword(inline_elements) if cf.is_non_generic()]

                        canonical_features.append(feature)
                else:
                    canonical_features = [feature]

                self.feature_matcher.add_keyword(inline_elements, canonical_features)

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

    def annotate(self, target_doc, targets, min_target_length):
        for target in targets:
            print('Annotating ' + target.identity + '...')

            if self.__has_min_length(target, min_target_length):
                inline_elements = ' '.join(target.nucleotides)
                rc_elements = ' '.join(target.reverse_complement_nucleotides())

                inline_matches = self.feature_matcher.extract_keywords(inline_elements, span_info=True)
                rc_matches = self.feature_matcher.extract_keywords(rc_elements, span_info=True)

                target_definition = target_doc.getComponentDefinition(target.identity)
                
                self.__process_feature_matches(target_doc, target_definition, inline_matches, SBOL_ORIENTATION_INLINE)
                self.__process_feature_matches(target_doc, target_definition, rc_matches, SBOL_ORIENTATION_REVERSE_COMPLEMENT,
                                       len(target.nucleotides) + 1)

                logging.info('Finished annotating %s.\n', target.identity)

class FeaturePruner():

    def __init__(self, feature_library, roles=set()):
        self.feature_library = feature_library
        self.roles = roles

    @classmethod
    def __is_covered(cls, anno, cover_annos, cover_offset):
        for cover_anno in cover_annos:
            if not abs(cover_anno[0] - anno[0]) <= cover_offset or not abs(cover_anno[1] - anno[1]) <= cover_offset:
                return False

        return True

    @classmethod
    def __remove_annotations(cls, indices, annos, target_definition):
        for i in range(len(annos) - 1, -1, -1):
            if annos[i][5] is None:
                feature_identity = annos[i][2]
            else:
                feature_identity = target_definition.components.get(annos[i][5]).definition
            
            if i in indices:
                target_definition.sequenceAnnotations.remove(annos[i][2])

                if annos[i][5] is not None:
                    target_definition.components.remove(annos[i][5])

                logging.info('Removed %s at [%s, %s] in %s.', feature_identity, annos[i][0], annos[i][1], target_definition.identity)

                del annos[i]

    @classmethod
    def __select_annotations(cls, doc, target_definition, annos):
        feature_messages = []

        for i in range(0, len(annos)):
            if annos[i][5] is None:
                if annos[i][4] is None:
                    feature_ID = annos[i][3]
                else:
                    feature_ID = annos[i][4]
            else:
                feature_identity = target_definition.components.get(annos[i][5]).definition
                feature_definition = doc.getComponentDefinition(feature_identity)
                
                if feature_definition.name is None:
                    feature_ID = feature_definition.displayId
                else:
                    feature_ID = feature_definition.name

            feature_messages.append('{nx}: {fi} at [{st}, {en}]'.format(nx=str(i), fi=feature_ID, st=annos[i][0], en=annos[i][1]))

        if target_definition.name is None:
            target_ID = target_definition.displayId
        else:
            target_ID = target_definition.name

        select_message = '\nThere appear to be redundant features in {pi}:\n{fm}\nPlease select which ones to remove if any (comma-separated list of indices, for example 0,2,5):\n'.format(fm='\n'.join(feature_messages), pi=target_ID)

        selected_message = input(select_message)

        try:
            selected_indices = [int(si.strip()) for si in selected_message.split(',')]
        except ValueError:
            selected_indices = []

        return set(selected_indices)

    @classmethod
    def __get_linked_definitions_identities(cls, doc, identity):
        linked_identities = {identity}

        try:
            component_definition = doc.getComponentDefinition(identity)
        except RuntimeError:
            pass

        for parent_identity in component_definition.wasDerivedFrom:
            linked_identities.add(parent_identity)

        for sub_comp in component_definition.components:
            linked_identities = linked_identities.union(cls.__get_linked_definitions_identities(doc, sub_comp.definition))

        return linked_identities

    def __filter_annotations(self, annos, target_definition):
        for i in range(len(annos) - 1, -1, -1):
            if annos[i][5] is None:
                feature_identity = annos[i][2]

                feature_roles = annos[i][6]
            else:
                feature_identity = target_definition.components.get(annos[i][5]).definition

                feature_roles = set(self.feature_library.get_definition(feature_identity).roles)

            if len(feature_roles.intersection(self.roles)) == 0:
                target_definition.sequenceAnnotations.remove(annos[i][2])

                if annos[i][5] is not None:
                    target_definition.components.remove(annos[i][5])

                logging.info('Removed %s at [%s, %s] in %s.', feature_identity, annos[i][0], annos[i][1], target_definition.identity)

                del annos[i]

    def __merge_annotations(self, anno, sub_anno, target_definition):
        feature_identity = target_definition.components.get(sub_anno[5]).definition

        feature_definition = self.feature_library.get_definition(feature_identity)

        if not Feature.has_non_generic_role(anno[6]) or Feature.has_same_roles(anno[6], set(feature_definition.roles)):
            if anno[4] == None:
                anno_ID = anno[3]
            else:
                anno_ID = anno[4]

            if feature_definition.name == None:
                feature_ID = feature_definition.displayId
            else:
                feature_ID = feature_definition.name

            print('\nMerging {ai} and {fi}...'.format(ai=anno_ID, fi=feature_ID))

            seq_anno = target_definition.sequenceAnnotations.get(anno[2])

            seq_anno.roles = []
            seq_anno.component = sub_anno[5]

            target_definition.sequenceAnnotations.remove(sub_anno[2])

            logging.info('Merged %s at [%s, %s] and %s at [%s, %s] in %s.', anno[2], anno[0], anno[1], feature_identity, sub_anno[0], sub_anno[1], target_definition.identity)

    def prune(self, target_doc, targets, cover_offset):
        for target in targets:
            target_definition = target_doc.getComponentDefinition(target.identity)

            cut_annos = [(sa.locations.getCut().at, sa.locations.getCut().at, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles)) for sa in target_definition.sequenceAnnotations if len(sa.locations) > 0 and sa.locations[0].getTypeURI() == SBOL_CUT]
            annos = [(sa.locations.getRange().start, sa.locations.getRange().end, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles)) for sa in target_definition.sequenceAnnotations if len(sa.locations) > 0 and sa.locations[0].getTypeURI() == SBOL_RANGE]
            
            annos.extend(cut_annos)

            if len(self.roles) > 0:
                self.__filter_annotations(annos, target_definition)

            annos.sort()

            grouped_annos = [[]]

            for anno in annos:
                if len(grouped_annos) > 1 and self.__is_covered(anno, grouped_annos[-2], cover_offset):
                    grouped_annos[-2].append(anno)
                elif self.__is_covered(anno, grouped_annos[-1], cover_offset):
                    grouped_annos[-1].append(anno)
                else:
                    grouped_annos.append([anno])

            for anno_group in grouped_annos:
                if len(anno_group) > 1:
                    selected_indices = self.__select_annotations(target_doc, target_definition, anno_group)

                    self.__remove_annotations(selected_indices, anno_group, target_definition)

            for anno_group in grouped_annos:
                if len(anno_group) == 2:
                    if anno_group[0][5] is None and anno_group[1][5] is not None:
                        self.__merge_annotations(anno_group[0], anno_group[1], target_definition)
                    elif anno_group[0][5] is not None and anno_group[1][5] is None:
                        self.__merge_annotations(anno_group[1], anno_group[0], target_definition)
             
            logging.info('Finished pruning %s.\n', target.identity)

        kept_identities = set()

        for target in targets:
            kept_identities = kept_identities.union(self.__get_linked_definitions_identities(target_doc, target.identity))

        doc_identities = set()

        for component_definition in target_doc.componentDefinitions:
            doc_identities.add(component_definition.identity)

        pruned_identities = doc_identities.difference(kept_identities)

        for pruned_identity in pruned_identities:
            target_doc.componentDefinitions.remove(pruned_identity)

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
    parser.add_argument('-c', '--cover_offset', nargs='?', default=14)
    parser.add_argument('-r', '--roles', nargs='*', default=[])
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

    feature_annotater = FeatureAnnotater(feature_library, int(args.min_feature_length))

    feature_pruner = FeaturePruner(feature_library, set(args.roles))

    for i in range (0, len(args.target_files)):
        target_doc = load_sbol(args.target_files[i])
        
        target_library = FeatureLibrary([target_doc], True)
        feature_annotater.annotate(target_doc, target_library.features, int(args.min_target_length))

        feature_pruner.prune(target_doc, target_library.features, int(args.cover_offset))

        (target_file_base, file_extension) = os.path.splitext(args.target_files[i])
        target_doc.write('_'.join([target_file_base, 'annotated']) + file_extension)

    print('Finished curating.')

if __name__ == '__main__':
    main()
