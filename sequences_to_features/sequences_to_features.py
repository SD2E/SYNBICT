import logging
import argparse
import os
import sys
import requests
import json

from Bio.Seq import Seq
from Bio import Align
import sbol2
from flashtext import KeywordProcessor

import time

def load_target_file(target_file):
    logger = logging.getLogger('synbict')

    if target_file.endswith('.xml') or target_file.endswith('.sbol'):
        return load_sbol(target_file)
    elif (target_file.endswith('.gb')
            or target_file.endswith('.genbank')
            or target_file.endswith('.fasta')
            or target_file.endswith('.faa')
            or target_file.endswith('.fa')
            or target_file.endswith('.fas')
            or target_file.endswith('.fsa')):
        return load_non_sbol(target_file)
    else:
        logger.error('Extension of target file %s is unrecognized.', target_file)

        return None

# Set up the not found error for catching
try:
    # SBOLError is in the native python module
    NotFoundError = sbol2.SBOLError
except NameError:
    # The swig wrapper raises RuntimeError on not found
    NotFoundError = RuntimeError

# Set up the not unique error for catching
try:
    # SBOLError is in the native python module
    NotUniqueError = sbol2.SBOLError
except NameError:
    # The swig wrapper raises RuntimeError on not unique
    NotUniqueError = RuntimeError

def is_sbol_not_found(exc):
    return (exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_NOT_FOUND
        or exc.error_code() == sbol2.SBOLErrorCode.NOT_FOUND_ERROR)

def load_sbol(sbol_file):
    logger = logging.getLogger('synbict')

    logger.info('Loading %s', sbol_file)

    doc = sbol2.Document()
    doc.read(sbol_file)

    doc.name = sbol_file

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')
    doc.addNamespace('http://sbolstandard.org/gff3#', 'gff3')
    doc.addNamespace('http://cellocad.org/Terms/cello#', 'cello')

    return doc

def load_non_sbol(non_sbol_file):
    logger = logging.getLogger('synbict')

    logger.info('Loading %s', non_sbol_file)

    conversion_request = {
        'options': {
            'language' : 'SBOL2',
            'test_equality': False,
            'check_uri_compliance': False,
            'check_completeness': False,
            'check_best_practices': False,
            'fail_on_first_error': False,
            'provide_detailed_stack_trace': False,
            'subset_uri': '',
            'uri_prefix': sbol2.getHomespace(),
            'version': '1',
            'insert_type': False,
            'main_file_name': 'main file',
            'diff_file_name': 'comparison file'
        },
        'return_file': True,
        'main_file': open(non_sbol_file).read()
    }

    conversion_response = requests.post("https://validator.sbolstandard.org/validate/", json=conversion_request)

    response_dict = json.loads(conversion_response.content.decode('utf-8'))

    doc = sbol2.Document()
    doc.readString(response_dict['result'])

    (file_base, file_extension) = os.path.splitext(non_sbol_file)
    doc.name = file_base + '.xml'

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')

    return doc

class FeatureCurator():

    def __init__(self, target_library, output_library=None):
        self.target_library = target_library
        self.output_library = output_library

        self.logger = logging.getLogger('synbict')

    def annotate_features(self, feature_annotater, min_target_length, in_place=False):
        start_time = time.clock()

        annotated_identities = feature_annotater.annotate(self.target_library, min_target_length, in_place,
                                                          self.output_library)

        self.logger.info('Annotation Time: ' + str(time.clock() - start_time))

        if self.output_library and len(self.output_library.docs) > 0:
            added_features = self.output_library.update(False)
        else:
            added_features = self.target_library.update()

        annotated_features = []
        annotating_features = []

        for added_feature in added_features:
            if added_feature.identity in annotated_identities:
                annotated_features.append(added_feature)
            else:
                annotating_features.append(added_feature)

        return (annotated_features, annotating_features)

    def prune_features(self, feature_pruner, cover_offset, min_target_length, target_features=[],
            target_sub_features=[], delete_flat=False, auto_swap=False, ask_user=True):
        if self.output_library and len(self.output_library.docs) > 0:
            feature_pruner.prune(self.output_library, cover_offset, min_target_length,
                                 ask_user=ask_user, delete_flat=delete_flat, target_features=target_features,
                                 auto_swap=auto_swap, require_sequence=False)
        else:
            feature_pruner.prune(self.target_library, cover_offset, min_target_length,
                                 ask_user=ask_user, delete_flat=delete_flat, target_features=target_features,
                                 auto_swap=auto_swap)

            feature_pruner.clean(self.target_library, target_features, target_sub_features)

    def extend_features(self, feature_annotater, min_target_length, extension_threshold):
        start_time = time.clock()

        feature_annotater.extend_features_by_name(self.target_library,
                                                  min_target_length,
                                                  extension_threshold)

        self.logger.info('Extension Time: ' + str(time.clock() - start_time))

class Feature():

    SO_REGION = 'http://identifiers.org/so/SO:0000001'
    SO_SEQUENCE_FEATURE = 'http://identifiers.org/so/SO:0000110'

    GENERIC_ROLES = {
        SO_REGION,
        SO_SEQUENCE_FEATURE
    }

    def __init__(self, nucleotides, identity, roles, sub_identities=[], parent_identities=[]):
        self.nucleotides = nucleotides
        self.identity = identity
        self.sub_identities = sub_identities
        self.parent_identities = parent_identities
        self.roles = set(roles)

        self.logger = logging.getLogger('synbict')

    def reverse_complement_nucleotides(self):
        return str(Seq(self.nucleotides).reverse_complement())

    @classmethod
    def has_non_generic_role(cls, roles):
        return len(roles.difference(cls.GENERIC_ROLES)) > 0

    def is_non_generic(self):
        return self.has_non_generic_role(self.roles)

class FeatureLibrary():

    def __init__(self, docs, require_sequence=True):
        self.features = []
        self.docs = docs
        self.__updated_indices = set()

        self.__feature_map = {}
        self.__feature_dict = {}
        self.__name_to_idents = {}

        self.logger = logging.getLogger('synbict')

        self.logger.info('Loading features')

        for i in range(0, len(self.docs)):
            self.__load_features(self.docs[i], i, require_sequence)

    def update(self, require_sequence=True):
        added_features = []

        for i in range(0, len(self.docs)):
            added_features.extend(self.__load_features(self.docs[i], i, require_sequence))

        for added_feature in added_features:
            self.__updated_indices.add(self.get_document_index(added_feature.identity))

        return added_features

    def get_updated_documents(self):
        updated_docs = []

        for updated_index in self.__updated_indices:
            updated_docs.append(self.docs[updated_index])

        return updated_docs

    def get_non_updated_indices(self):
        non_updated_indices = []

        for i in range(0, len(self.docs)):
            if i not in self.__updated_indices:
                non_updated_indices.append(i)

        return non_updated_indices

    def __load_features(self, doc, doc_index, require_sequence=True):
        loaded_features = []

        comp_seq_identities = set()

        for comp_definition in doc.componentDefinitions:
            if sbol2.BIOPAX_DNA in comp_definition.types:
                dna_seqs = self.get_DNA_sequences(comp_definition, doc)

                for dna_seq in dna_seqs:
                    comp_seq_identities.add(dna_seq.identity)

                if comp_definition.identity not in self.__feature_map:
                    sub_identities = []

                    for sub_comp in comp_definition.components:
                        sub_identities.append(sub_comp.definition)

                    if len(dna_seqs) > 0:
                        feature = Feature(dna_seqs[0].elements,
                                          comp_definition.identity,
                                          comp_definition.roles,
                                          sub_identities,
                                          comp_definition.wasDerivedFrom)

                        loaded_features.append(feature)
                        self.features.append(feature)

                        self.__feature_map[comp_definition.identity] = doc_index
                        self.__feature_dict[comp_definition.identity] = feature

                        if comp_definition.name:
                            if comp_definition.name not in self.__name_to_idents:
                                self.__name_to_idents[comp_definition.name] = []

                            self.__name_to_idents[comp_definition.name].append(comp_definition.identity)
                    elif not require_sequence:
                        feature = Feature('',
                                          comp_definition.identity,
                                          comp_definition.roles,
                                          sub_identities,
                                          comp_definition.wasDerivedFrom)

                        loaded_features.append(feature)
                        self.features.append(feature)

                        self.__feature_map[comp_definition.identity] = doc_index
                        self.__feature_dict[comp_definition.identity] = feature

                        if comp_definition.name:
                            if comp_definition.name not in self.__name_to_idents:
                                self.__name_to_idents[comp_definition.name] = []

                            self.__name_to_idents[comp_definition.name].append(comp_definition.identity)
                    else:
                        self.logger.warning('%s not loaded since its DNA sequence was not found', comp_definition.identity)

        for seq in doc.sequences:
            if seq.identity not in comp_seq_identities and seq.encoding == sbol2.SBOL_ENCODING_IUPAC:
                seq_comp_definition = sbol2.ComponentDefinition(seq.displayId + '_comp', sbol2.BIOPAX_DNA, '1')
                seq_comp_definition.sequences = [seq.identity]

                try:
                    doc.addComponentDefinition(seq_comp_definition)

                    feature = Feature(seq.elements,
                                      seq_comp_definition.identity,
                                      [],
                                      [],
                                      [])

                    loaded_features.append(feature)
                    self.features.append(feature)

                    self.__feature_map[seq_comp_definition.identity] = doc_index
                    self.__feature_dict[seq_comp_definition.identity] = feature
                except RuntimeError:
                    self.logger.warning('Component could not be automatically generated for DNA sequence %s', seq.identity)
                except NotUniqueError as exc:
                    if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                        self.logger.warning('Component could not be automatically generated for DNA sequence %s', seq.identity)
                    else:
                        raise

        return loaded_features

    def get_features(self, min_feature_length=0, children_only=False):
        features = []

        if children_only:
            parent_identities = set()

            for feature in self.features:
                for parent_identity in feature.parent_identities:
                    parent_identities.add(parent_identity)

            for feature in self.features:
                if (min_feature_length == 0 or len(feature.nucleotides) > min_feature_length) and feature.identity not in parent_identities:
                    features.append(feature)
        else:
            for feature in self.features:
                if min_feature_length == 0 or len(feature.nucleotides) > min_feature_length:
                    features.append(feature)

        return features

    def get_added_feature_identities(self):
        added_feature_identities = set()

        for doc in self.docs:
            for comp_definition in doc.componentDefinitions:
                if comp_definition.identity not in self.__feature_map:
                    added_feature_identities.append(comp_definition.identity)

        return added_feature_identities

    def get_document(self, identity):
        return self.docs[self.get_document_index(identity)]

    def get_document_index(self, identity):
        if identity in self.__feature_map:
            return self.__feature_map[identity]
        else:
            return -1

    def get_definition(self, identity):
        return self.get_document(identity).getComponentDefinition(identity)

    def get_definitions_by_name(self, name):
        name_keys = []

        if name in self.__name_to_idents:
            name_keys.append(name)
        else:
            for other_name in self.__name_to_idents:
                if name in other_name or other_name in name:
                    name_keys.append(other_name)

        definitions = []

        for name_key in name_keys:
            identities = self.__name_to_idents[name_key]

            for identity in identities:
                definitions.append(self.get_definition(identity))

        return definitions

    def has_feature(self, identity):
        return identity in self.__feature_map

    def get_feature(self, identity):
        return self.__feature_dict[identity]

    @classmethod
    def get_DNA_sequences(cls, comp_definition, doc):
        dna_seqs = []

        for seq_URI in comp_definition.sequences:
            try:
                seq = doc.getSequence(seq_URI)
            except RuntimeError:
                seq = None
            except NotFoundError as exc:
                if is_sbol_not_found(exc):
                    seq = None
                else:
                    raise

            if seq and seq.encoding == sbol2.SBOL_ENCODING_IUPAC:
                dna_seqs.append(seq)

        return dna_seqs

    @classmethod
    def get_sequences(cls, comp_definition, doc):
        seqs = []

        for seq_URI in comp_definition.sequences:
            try:
                seq = doc.getSequence(seq_URI)
            except RuntimeError:
                seq = None
            except NotFoundError as exc:
                if is_sbol_not_found(exc):
                    seq = None
                else:
                    raise

            if seq:
                seqs.append(seq)

        return seqs

    @classmethod
    def copy_sequence(cls, seq, source_doc, sink_doc, import_namespace=False):
        if import_namespace:
            namespace = '/'.join(seq.identity.split('/')[:-2])

            if namespace == sbol2.getHomespace():
                try:
                    version = int(seq.version)
                except (TypeError, ValueError):
                    return None

                try:
                    seq_copy = seq.copy(sink_doc, namespace, str(version + 1))

                except RuntimeError:
                    return sink_doc.getSequence('/'.join([sbol2.getHomespace(), seq.displayId,
                                                          str(version + 1)]))
                except NotUniqueError as exc:
                    if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                        return sink_doc.getSequence('/'.join([sbol2.getHomespace(), seq.displayId,
                                                              str(version + 1)]))
                    else:
                        raise
                    
            else:
                try:
                    seq_copy = seq.copy(sink_doc, namespace, '1')
                except RuntimeError:
                    return sink_doc.getSequence('/'.join([sbol2.getHomespace(), seq.displayId, '1']))
                except NotUniqueError as exc:
                    if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                        return sink_doc.getSequence('/'.join([sbol2.getHomespace(), seq.displayId, '1']))
                    else:
                        raise

            cls.strip_origin_properties(seq_copy)
        else:
            try:
                sink_doc.getSequence(seq.identity)
                
                return None
            except RuntimeError:
                seq_copy = seq.copy(sink_doc)
            except NotFoundError as exc:
                if is_sbol_not_found(exc):
                    seq_copy = seq.copy(sink_doc)
                else:
                    raise

        return seq_copy

    @classmethod
    def make_variant_definition(cls, doc, definition_copy):
        doc.componentDefinitions.remove(definition_copy.identity)

        variant_index = 1
        unique_flag = False
        
        while not unique_flag:
            variant_ID = '_'.join([definition_copy.displayId, 'v' + str(variant_index)])

            split_identity = definition_copy.identity.split('/')
            variant_identity = '/'.join(split_identity[:-2] + [variant_ID, split_identity[-1]])
            variant_p_identity = '/'.join(split_identity[:-2] + [variant_ID])

            original_identity = definition_copy.identity
            original_ID = definition_copy.displayId
            original_p_identity = definition_copy.persistentIdentity

            definition_copy.identity = variant_identity
            definition_copy.displayId = variant_ID
            definition_copy.persistentIdentity = variant_p_identity

            try:
                doc.componentDefinitions.add(definition_copy)

                unique_flag = True
            except RuntimeError:
                definition_copy.identity = original_identity
                definition_copy.displayId = original_ID
                definition_copy.persistentIdentity = original_p_identity

                variant_index = variant_index + 1

                unique_flag = False
            except NotUniqueError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    definition_copy.identity = original_identity
                    definition_copy.displayId = original_ID
                    definition_copy.persistentIdentity = original_p_identity

                    variant_index = variant_index + 1

                    unique_flag = False
                else:
                    raise

    @classmethod
    def make_variant_sequence(cls, doc, sequence_copy):
        doc.sequences.remove(sequence_copy.identity)

        variant_index = 1
        unique_flag = False
        
        while not unique_flag:
            variant_ID = '_'.join([sequence_copy.displayId, 'v' + str(variant_index)])

            split_identity = sequence_copy.identity.split('/')
            variant_identity = '/'.join(split_identity[:-2] + [variant_ID, split_identity[-1]])
            variant_p_identity = '/'.join(split_identity[:-2] + [variant_ID])

            original_identity = sequence_copy.identity
            original_ID = sequence_copy.displayId
            original_p_identity = sequence_copy.persistentIdentity

            sequence_copy.identity = variant_identity
            sequence_copy.displayId = variant_ID
            sequence_copy.persistentIdentity = variant_p_identity

            try:
                doc.sequences.add(sequence_copy)

                unique_flag = True
            except RuntimeError:
                sequence_copy.identity = original_identity
                sequence_copy.displayId = original_ID

                variant_index = variant_index + 1

                unique_flag = False
            except NotUniqueError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    sequence_copy.identity = original_identity
                    sequence_copy.displayId = original_ID

                    variant_index = variant_index + 1

                    unique_flag = False
                else:
                    raise

    @classmethod
    def strip_origin_properties(cls, sbol_obj):
        origin_props = []

        for prop in sbol_obj.properties:
            if (prop.startswith('http://wiki.synbiohub.org/wiki/Terms/synbiohub#')
                    or prop.startswith('http://purl.org/dc/terms/created')
                    or prop.startswith('http://purl.org/dc/terms/modified')
                    or prop.startswith('http://www.ncbi.nlm.nih.gov/genbank#')
                    or prop.startswith('http://sbols.org/genBankConversion#')):
                origin_props.append(prop)

        for origin_prop in origin_props:
            del sbol_obj.properties[origin_prop]

        sbol_obj.wasGeneratedBy = []

    @classmethod
    def copy_component_definition(cls, comp_definition, source_doc, sink_doc, import_namespace=False,
                                  min_seq_length=0, import_sequences=False, seq_elements=None,
                                  parent_definitions=[], parent_doc=None, make_variant=False,
                                  shallow_copy=False):
        if sbol2.BIOPAX_DNA in comp_definition.types:
            seqs = cls.get_DNA_sequences(comp_definition, source_doc)
        else:
            seqs = cls.get_sequences(comp_definition, source_doc)

        if min_seq_length == 0 or (len(seqs) > 0 and len(seqs[0].elements) >= min_seq_length):
            namespace = '/'.join(comp_definition.identity.split('/')[:-2])

            if import_namespace:
                if namespace == sbol2.getHomespace():
                    try:
                        version = int(comp_definition.version)
                    except (TypeError, ValueError):
                        return None

                    try:
                        definition_copy = comp_definition.copy(sink_doc, namespace, str(version + 1))
                    except RuntimeError:
                        return sink_doc.getComponentDefinition('/'.join([sbol2.getHomespace(),
                                                                         comp_definition.displayId,
                                                                         str(version + 1)]))
                    except NotUniqueError as exc:
                        if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                            return sink_doc.getComponentDefinition('/'.join([sbol2.getHomespace(),
                                                                             comp_definition.displayId,
                                                                             str(version + 1)]))
                        else:
                            raise
                        
                else:
                    try:
                        definition_copy = comp_definition.copy(sink_doc, namespace, '1')
                    except RuntimeError:
                        return sink_doc.getComponentDefinition('/'.join([sbol2.getHomespace(),
                                                                         comp_definition.displayId, '1']))
                    except NotUniqueError as exc:
                        if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                            return sink_doc.getComponentDefinition('/'.join([sbol2.getHomespace(),
                                                                             comp_definition.displayId, '1']))
                        else:
                            raise

                cls.strip_origin_properties(definition_copy)
                for sub_comp_copy in definition_copy.components:
                    cls.strip_origin_properties(sub_comp_copy)
                for anno_copy in definition_copy.sequenceAnnotations:
                    cls.strip_origin_properties(anno_copy)

                    for loc_copy in anno_copy.locations:
                        cls.strip_origin_properties(loc_copy)

                if make_variant:
                    cls.make_variant_definition(sink_doc, definition_copy)
            else:
                try:
                    sink_doc.getComponentDefinition(comp_definition.identity)

                    return None
                except RuntimeError:
                    definition_copy = comp_definition.copy(sink_doc)
                except NotFoundError as exc:
                    if is_sbol_not_found(exc):
                        definition_copy = comp_definition.copy(sink_doc)
                    else:
                        raise

            if shallow_copy:
                definition_copy.sequences = list(comp_definition.sequences)
            elif import_sequences:
                if len(seqs) > 0:
                    seq_copy = cls.copy_sequence(seqs[0], source_doc, sink_doc, True)

                    if make_variant:
                        cls.make_variant_sequence(sink_doc, seq_copy)

                    if seq_elements:
                        seq_copy.elements = seq_elements

                    if parent_doc:
                        for parent_definition in parent_definitions:
                            if sbol2.BIOPAX_DNA in parent_definition.types:
                                parent_seqs = cls.get_DNA_sequences(parent_definition, parent_doc)

                                if len(parent_seqs) > 0:
                                    seq_copy.wasDerivedFrom = seq_copy.wasDerivedFrom + [parent_seqs[0].identity]
                            elif len(parent_definition.sequences) > 0:
                                seq_copy.wasDerivedFrom = seq_copy.wasDerivedFrom + [parent_definition.sequences[0].identity]

                    definition_copy.sequences = [seq_copy.identity]
                else:
                    return None
            else:
                for seq_URI in comp_definition.sequences:
                    seq = source_doc.getSequence(seq_URI)

                    cls.copy_sequence(seq, source_doc, sink_doc)

                definition_copy.sequences = list(comp_definition.sequences)

            if make_variant:
                definition_copy.sequenceAnnotations = []
                definition_copy.components = []
            else:
                for seq_anno in comp_definition.sequenceAnnotations:
                    if seq_anno.component:
                        sub_comp = comp_definition.components.get(seq_anno.component)

                        sub_copy = definition_copy.components.get(sub_comp.displayId)

                        anno_copy = definition_copy.sequenceAnnotations.get(seq_anno.displayId)
                        anno_copy.component = sub_copy.identity

                for sub_comp in comp_definition.components:
                    try:
                        sub_definition = source_doc.getComponentDefinition(sub_comp.definition)
                    except RuntimeError:
                        sub_definition = None
                    except NotFoundError as exc:
                        if is_sbol_not_found(exc):
                            sub_definition = None
                        else:
                            raise

                    sub_copy = definition_copy.components.get(sub_comp.displayId)

                    if sub_definition:
                        if shallow_copy:
                            sub_copy.definition = sub_definition.identity
                        else:
                            sub_definition_copy = cls.copy_component_definition(sub_definition, source_doc, sink_doc,
                                                                                import_namespace, min_seq_length,
                                                                                shallow_copy=shallow_copy)

                            if sub_definition_copy:
                                sub_copy.definition = sub_definition_copy.identity
                            else:
                                sub_copy.definition = sub_definition.identity
                    else:
                        sub_copy.definition = sub_comp.definition

            for parent_definition in parent_definitions:
                definition_copy.wasDerivedFrom = definition_copy.wasDerivedFrom + [parent_definition.identity]

            return definition_copy
        else:
            return None

class FeatureAnnotater():

    def __init__(self, feature_library, min_feature_length):
        self.feature_library = feature_library
        self.feature_matcher = KeywordProcessor()

        self.logger = logging.getLogger('synbict')

        for feature in feature_library.features:
            inline_elements = ' '.join(feature.nucleotides)

            if self.__has_min_length(feature, min_feature_length):
                if inline_elements in self.feature_matcher:
                    if feature.is_non_generic():
                        canonical_features = [cf for cf in self.feature_matcher.get_keyword(inline_elements) if
                                              cf.is_non_generic()]

                        canonical_features.append(feature)
                else:
                    canonical_features = [feature]

                self.feature_matcher.add_keyword(inline_elements, canonical_features)

    def get_updated_documents(self):
        return self.feature_library.get_updated_documents()

    @classmethod
    def __has_min_length(cls, feature, min_feature_length):
        return min_feature_length == 0 or len(feature.nucleotides) >= min_feature_length

    @classmethod
    def __create_sub_component(cls, parent_definition, child_definition):
        i = 1

        while i > 0:
            try:
                sub_comp = parent_definition.components.create('_'.join([child_definition.displayId,
                                                                         'comp',
                                                                         str(i)]))
            except RuntimeError:
                sub_comp = None
            except NotUniqueError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    sub_comp = None
                else:
                    raise

            if sub_comp is None:
                i = i + 1
            else:
                sub_comp.name = child_definition.name
                sub_comp.definition = child_definition.identity

                sub_comp.roleIntegration = None

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
            except NotUniqueError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    seq_anno = None
                else:
                    raise

            if seq_anno is None:
                i = i + 1
            else:
                seq_anno.name = child_definition.name
                seq_anno.description = child_definition.description
                if sub_comp_URI:
                    seq_anno.component = sub_comp_URI
                if parent_URI:
                    seq_anno.roles = seq_anno.roles + child_definition.roles
                    seq_anno.wasDerivedFrom = seq_anno.wasDerivedFrom + [parent_URI]
                
                location = seq_anno.locations.createRange('_'.join([seq_anno.displayId,
                                                                    'loc']))

                location.orientation = orientation
                location.start = start
                location.end = end

                i = -1

        return seq_anno

    def __process_feature_matches(self, target_doc, target_definition, feature_matches, orientation, target_length,
                                  rc_factor=0, copy_definitions=True):
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
                if len(feature.nucleotides) < target_length:
                    feature_definition = self.feature_library.get_definition(feature.identity)

                    if feature_definition.name is None:
                        feature_ID = feature_definition.displayId
                    else:
                        feature_ID = feature_definition.name

                    feature_role = FeaturePruner.get_common_role(feature_definition.roles)

                    sub_comp = self.__create_sub_component(target_definition, feature_definition)
                    self.__create_sequence_annotation(target_definition, feature_definition, orientation, start, end,
                                                      sub_comp.identity)

                    if copy_definitions:
                        feature_doc = self.feature_library.get_document(feature.identity)

                        FeatureLibrary.copy_component_definition(feature_definition, feature_doc, target_doc)

                    self.logger.debug('Annotated %s (%s, %s) at [%s, %s] in %s', feature_definition.identity, feature_ID,
                                      feature_role, start, end, target_definition.identity)

    def extend_features_by_name(self, target_library, min_target_length, mismatch_threshold):
        self.logger.info('Extending feature library')

        aligner = Align.PairwiseAligner()
        aligner.match_score = 1
        aligner.mismatch_score = -2
        aligner.internal_gap_score = -2.5

        for target in target_library.features:
            if self.__has_min_length(target, min_target_length):
                target_doc = target_library.get_document(target.identity)

                target_definition = target_doc.getComponentDefinition(target.identity)

                for seq_anno in target_definition.sequenceAnnotations:
                    if (seq_anno.name and not seq_anno.component and len(seq_anno.locations) == 1
                            and seq_anno.locations[0].getTypeURI() == sbol2.SBOL_RANGE):
                        anno_start = seq_anno.locations.getRange().start
                        anno_end = seq_anno.locations.getRange().end

                        target_nucleotides = target.nucleotides[anno_start - 1:anno_end]
                        rc_target_nucleotides = str(Seq(target_nucleotides).reverse_complement())

                        inline_elements = ' '.join(target_nucleotides)
                        rc_elements = ' '.join(rc_target_nucleotides)

                        if not inline_elements in self.feature_matcher and not rc_elements in self.feature_matcher:
                            feature_definitions = self.feature_library.get_definitions_by_name(seq_anno.name)

                            for feature_definition in feature_definitions:
                                if set(seq_anno.roles) == set(feature_definition.roles):
                                    feature_doc = self.feature_library.get_document(feature_definition.identity)

                                    feature_seqs = FeatureLibrary.get_DNA_sequences(feature_definition, feature_doc)
                                    
                                    if len(feature_seqs) > 0:
                                        feature_nucleotides = feature_seqs[0].elements

                                        score = aligner.score(target_nucleotides, feature_nucleotides)
                                        rc_score = aligner.score(rc_target_nucleotides, feature_nucleotides)

                                        if rc_score > score:
                                            best_score = rc_score
                                            best_nucleotides = rc_target_nucleotides
                                        else:
                                            best_score = score
                                            best_nucleotides = target_nucleotides

                                        if len(target_nucleotides) < len(feature_nucleotides):
                                            max_score = len(target_nucleotides)
                                        else:
                                            max_score = len(feature_nucleotides)

                                        if max_score - best_score < mismatch_threshold*max_score:
                                            variant_definition = FeatureLibrary.copy_component_definition(feature_definition,
                                                feature_doc, feature_doc, import_namespace=True, import_sequences=True,
                                                seq_elements=target_nucleotides, parent_definitions=[target_definition],
                                                parent_doc=target_doc, make_variant=True)

                                            if variant_definition:
                                                sub_identities = []
                                                for sub_comp in variant_definition.components:
                                                    sub_identities.append(sub_comp.definition)

                                                feature = Feature(feature_nucleotides,
                                                                  variant_definition.identity,
                                                                  variant_definition.roles,
                                                                  sub_identities,
                                                                  variant_definition.wasDerivedFrom)

                                                self.feature_matcher.add_keyword(inline_elements, [feature])

                                                self.logger.debug('Extended feature library with %s', variant_definition.identity)

        self.feature_library.update()

        self.logger.info('Finished extending feature library')

    def annotate_raw_sequences(self, raw_seqs, comp_IDs=[], min_target_length=0):
        annotated_comps = []

        if not isinstance(raw_seqs, list):
            raw_seqs = [raw_seqs]

        if not isinstance(comp_IDs, list):
            comp_IDs = [comp_IDs]

        for i in range(0, len(raw_seqs)):
            target_doc = sbol2.Document()

            if i < len(comp_IDs):
                comp_ID = comp_IDs[i]
            else:
                comp_ID = 'construct_' + str(i + 1)

            target_comp = sbol2.ComponentDefinition(comp_ID, sbol2.BIOPAX_DNA, '1')
            target_comp.sequence = sbol2.Sequence(comp_ID + '_seq', raw_seqs[i], sbol2.SBOL_ENCODING_IUPAC, '1')

            annotated_comps.append(target_comp)

            target_doc.addComponentDefinition(target_comp)

            target_library = FeatureLibrary([target_doc])

            self.annotate(target_library, min_target_length, True)

        if len(annotated_comps) == 1:
            return annotated_comps[0]
        else:
            return annotated_comps

    def annotate(self, target_library, min_target_length, in_place=False, output_library=None):
        annotated_identities = []

        for target in target_library.features:
            if self.__has_min_length(target, min_target_length):
                self.logger.info('Annotating %s', target.identity)

                inline_elements = ' '.join(target.nucleotides)
                rc_elements = ' '.join(target.reverse_complement_nucleotides())

                inline_matches = self.feature_matcher.extract_keywords(inline_elements, span_info=True)
                rc_matches = self.feature_matcher.extract_keywords(rc_elements, span_info=True)

                if len(inline_matches) > 0 or len(rc_matches) > 0:
                    target_doc = target_library.get_document(target.identity)

                    target_definition = target_doc.getComponentDefinition(target.identity)

                    doc_index = target_library.get_document_index(target.identity)
                    
                    if output_library and doc_index < len(output_library.docs):
                        output_doc = output_library.docs[doc_index]

                        if in_place:
                            definition_copy = FeatureLibrary.copy_component_definition(target_definition,
                                                                                       target_doc,
                                                                                       output_doc,
                                                                                       min_seq_length=min_target_length,
                                                                                       shallow_copy=True)
                        else:
                            definition_copy = FeatureLibrary.copy_component_definition(target_definition,
                                                                                       target_doc,
                                                                                       output_doc, True,
                                                                                       min_target_length,
                                                                                       shallow_copy=True)
                    elif in_place:
                        definition_copy = target_definition
                    else:
                        definition_copy = FeatureLibrary.copy_component_definition(target_definition, target_doc,
                                                                                   target_doc, True,
                                                                                   min_target_length)

                    if definition_copy:
                        self.__process_feature_matches(target_doc, definition_copy, inline_matches,
                            sbol2.SBOL_ORIENTATION_INLINE, len(target.nucleotides),
                            copy_definitions=(not output_library or doc_index >= len(output_library.docs)))
                        self.__process_feature_matches(target_doc, definition_copy, rc_matches,
                            sbol2.SBOL_ORIENTATION_REVERSE_COMPLEMENT, len(target.nucleotides), len(target.nucleotides) + 1,
                            (not output_library or doc_index >= len(output_library.docs)))

                        annotated_identities.append(definition_copy.identity)
                    else:
                        self.logger.warning('%s was not annotated because its version could not be incremented.',
                                        target.identity)

                self.logger.info('Finished annotating %s', target.identity)

        return annotated_identities

class FeaturePruner():

    COMMON_ROLE_DICT = {
        sbol2.SO_PROMOTER: 'promoter',
        sbol2.SO_CDS: 'CDS',
        sbol2.SO_TERMINATOR: 'terminator',
        'http://identifiers.org/so/SO:0001977': 'ribonuclease_site'
    }

    def __init__(self, feature_library, roles=set()):
        self.feature_library = feature_library
        self.roles = roles

        self.logger = logging.getLogger('synbict')

    @classmethod
    def __has_min_length(cls, feature, min_feature_length):
        return min_feature_length == 0 or len(feature.nucleotides) >= min_feature_length

    @classmethod
    def __is_covered(cls, anno, cover_annos, cover_offset):
        for cover_anno in cover_annos:
            if not abs(cover_anno[0] - anno[0]) <= cover_offset or not abs(cover_anno[1] - anno[1]) <= cover_offset:
                return False

        return True

    def __remove_annotations(self, indices, annos, target_definition):
        for i in range(len(annos) - 1, -1, -1):
            if annos[i][5] is None:
                feature_identity = annos[i][2]
            else:
                feature_identity = target_definition.components.get(annos[i][5]).definition
            
            if i in indices:
                target_definition.sequenceAnnotations.remove(annos[i][2])

                if annos[i][5]:
                    target_definition.components.remove(annos[i][5])

                self.logger.debug('Removed %s at [%s, %s] in %s', feature_identity, annos[i][0], annos[i][1],
                              target_definition.identity)

                del annos[i]

    @classmethod
    def get_common_role(cls, roles):
        for role in roles:
            if role in cls.COMMON_ROLE_DICT:
                return cls.COMMON_ROLE_DICT[role]

        return 'sequence_feature'

    @classmethod
    def __ask_swap_annotations(cls, anno, sub_part_anno, target_definition):
        if anno[4] is None:
            anno_ID = anno[3]
        else:
            anno_ID = anno[4]

        if sub_part_anno[4] is None:
            sub_part_anno_ID = sub_part_anno[3]
        else:
            sub_part_anno_ID = sub_part_anno[4]

        select_message = 'Annotation {an} ({ai}) and\nsub-part annotation {sp} ({si})\nat [{st}, {en}] in {td} appear \
to be nearly identical.\nRemove second annotation and link first to sub-part? \
(0=no,1=yes):\n'.format(an=anno[2],
                        ai=anno_ID,
                        sp=sub_part_anno[2],
                        si=sub_part_anno_ID,
                        st=anno[0],
                        en=anno[1],
                        td=target_definition.identity)

        selected_message = input(select_message)

        try:
            return int(selected_message.strip()) == 1
        except ValueError:
            return False

    @classmethod
    def __select_annotations(cls, doc, target_definition, annos, ask_user=True, canonical_library=None, keep_flat=True):
        kept_indices = []

        feature_messages = []

        for i in range(0, len(annos)):
            if annos[i][5] is None:
                if annos[i][4] is None:
                    feature_ID = annos[i][3]
                else:
                    feature_ID = annos[i][4]

                feature_role = cls.get_common_role(annos[i][6])

                if ask_user:
                    if annos[i][7] and len(annos[i][7]) > 0:
                        if len(feature_role) > 0:
                            feature_messages.append('{nx}: {id} ({fi}, {ro}) at [{st}, {en}]. {de}'.format(nx=str(i), id=annos[i][2],
                                fi=feature_ID, ro=feature_role, st=annos[i][0], en=annos[i][1], de=annos[i][7]))
                        else:
                            feature_messages.append('{nx}: {id} ({fi}) at [{st}, {en}]. {de}'.format(nx=str(i), id=annos[i][2],
                                fi=feature_ID, st=annos[i][0], en=annos[i][1], de=annos[i][7]))
                    elif len(feature_role) > 0:
                        feature_messages.append('{nx}: {id} ({fi}, {ro}) at [{st}, {en}]'.format(nx=str(i), id=annos[i][2],
                            fi=feature_ID, ro=feature_role, st=annos[i][0], en=annos[i][1]))
                    else:
                        feature_messages.append('{nx}: {id} ({fi}) at [{st}, {en}]'.format(nx=str(i), id=annos[i][2],
                            fi=feature_ID, st=annos[i][0], en=annos[i][1]))
                elif keep_flat:
                    kept_indices.append(i)
            else:
                feature_identity = target_definition.components.get(annos[i][5]).definition

                if ask_user:
                    feature_definition = doc.getComponentDefinition(feature_identity)
                    
                    if feature_definition.name is None:
                        feature_ID = feature_definition.displayId
                    else:
                        feature_ID = feature_definition.name

                    feature_role = cls.get_common_role(feature_definition.roles)

                    feature_description = feature_definition.description

                    if feature_description and len(feature_description) > 0:
                        if len(feature_role) > 0:
                            feature_messages.append('{nx}: {id} ({fi}, {ro}) at [{st}, {en}]. {de}'.format(nx=str(i), id=feature_identity,
                                fi=feature_ID, ro=feature_role, st=annos[i][0], en=annos[i][1], de=feature_description))
                        else:
                            feature_messages.append('{nx}: {id} ({fi}) at [{st}, {en}]. {de}'.format(nx=str(i), id=feature_identity,
                                fi=feature_ID, st=annos[i][0], en=annos[i][1], de=feature_description))
                    elif len(feature_role) > 0:
                        feature_messages.append('{nx}: {id} ({fi}, {ro}) at [{st}, {en}]'.format(nx=str(i), id=feature_identity,
                            fi=feature_ID, ro=feature_role, st=annos[i][0], en=annos[i][1]))
                    else:
                        feature_messages.append('{nx}: {id} ({fi}) at [{st}, {en}]'.format(nx=str(i), id=feature_identity,
                            fi=feature_ID, st=annos[i][0], en=annos[i][1]))
                elif not canonical_library or canonical_library.has_feature(feature_identity):
                    kept_indices.append(i)

        if ask_user:
            if target_definition.name is None:
                target_ID = target_definition.displayId
            else:
                target_ID = target_definition.name

            select_message = 'There appear to be redundant features in {pi}:\n{fm}\nPlease select which ones to \
remove if any (comma-separated list of indices, for example 0,2,5):\n'.format(fm='\n'.join(feature_messages),
                                                                              pi=target_ID)

            selected_message = input(select_message)

            try:
                selected_indices = [int(si.strip()) for si in selected_message.split(',')]
            except ValueError:
                selected_indices = []

            return set(selected_indices)
        else:
            return set(range(0, len(annos))).difference(set(kept_indices))

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

                if annos[i][5]:
                    target_definition.components.remove(annos[i][5])

                self.logger.debug('Removed %s at [%s, %s] in %s', feature_identity, annos[i][0], annos[i][1],
                              target_definition.identity)

                del annos[i]

    def __are_swappable(self, anno, sub_part_anno, target_definition):
        if not anno[5] and sub_part_anno[5] and anno[0] == sub_part_anno[0] and anno[1] == sub_part_anno[1]:
            feature_identity = target_definition.components.get(sub_part_anno[5]).definition

            feature_definition = self.feature_library.get_definition(feature_identity)

            return (not Feature.has_non_generic_role(anno[6]) or
                    anno[6] == set(feature_definition.roles))
        else:
            return False

    def __swap_annotations(self, anno, sub_part_anno, target_definition):
        seq_anno = target_definition.sequenceAnnotations.get(anno[2])

        seq_anno.roles = []
        seq_anno.component = sub_part_anno[5]

        target_definition.sequenceAnnotations.remove(sub_part_anno[2])

        self.logger.debug('Removed %s at [%s, %s] in %s', sub_part_anno[2], sub_part_anno[0], sub_part_anno[1],
                      target_definition.identity)
        self.logger.debug('Modified %s at [%s, %s] in %s to refer to %s', anno[2], anno[0], anno[1],
                      target_definition.identity, sub_part_anno[5])

    # @classmethod
    # def __get_annotations(cls, doc, comp_definition):
    #     cut_annos = [(sa.locations.getCut().at, sa.locations.getCut().at, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles)) for sa in comp_definition.sequenceAnnotations if len(sa.locations) == 1 and sa.locations[0].getTypeURI() == SBOL_CUT]
    #     annos = [(sa.locations.getRange().start, sa.locations.getRange().end, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles)) for sa in comp_definition.sequenceAnnotations if len(sa.locations) == 1 and sa.locations[0].getTypeURI() == SBOL_RANGE]
            
    #     annos.extend(cut_annos)

    #     for sub_comp in comp_definition.components:
    #         try:
    #             sub_definition = doc.getComponentDefinition(comp_definition.identity)
    #         except RuntimeError:
    #             sub_definition = None

    #         if sub_definition is not None:
    #             annos.extend(cls.__get_annotations(doc, sub_definition))

    #     return annos

    @classmethod
    def __get_flat_annotation_indices(cls, anno_group):
        flat_indices = []

        for i in range(0, len(anno_group)):
            if anno_group[i][5] is None:
                flat_indices.append(i)

        return flat_indices

    def clean(self, feature_library, annotated_features, annotating_features):
        self.logger.info('Cleaning up')

        sub_definitions = set()

        for annotated_feature in annotated_features:
            annotated_doc = feature_library.get_document(annotated_feature.identity)
            annotated_definition = annotated_doc.getComponentDefinition(annotated_feature.identity)

            sub_IDs = set()
            temp_sub_definitions = set()

            for comp in annotated_definition.components:
                sub_IDs.add(comp.displayId)
                temp_sub_definitions.add(comp.definition)
            for seq_anno in annotated_definition.sequenceAnnotations:
                sub_IDs.add(seq_anno.displayId)

            parent_sub_IDs = set()

            for parent_identity in annotated_definition.wasDerivedFrom:
                parent_doc = feature_library.get_document(parent_identity)

                if parent_doc:
                    parent_definition = parent_doc.getComponentDefinition(parent_identity)

                    for comp in parent_definition.components:
                        parent_sub_IDs.add(comp.displayId)
                    for seq_anno in parent_definition.sequenceAnnotations:
                        parent_sub_IDs.add(seq_anno.displayId)

            if len(sub_IDs) == len(parent_sub_IDs) and len(sub_IDs) == len(sub_IDs.intersection(parent_sub_IDs)):
                annotated_doc.componentDefinitions.remove(annotated_feature.identity)

                self.logger.debug('Removed %s from %s', annotated_feature.identity, annotated_doc.name)
            else:
                sub_definitions.update(temp_sub_definitions)

        for annotating_feature in annotating_features:
            if annotating_feature.identity not in sub_definitions:
                annotating_doc = feature_library.get_document(annotating_feature.identity)
                annotating_definition = annotating_doc.getComponentDefinition(annotating_feature.identity)

                annotating_doc.componentDefinitions.remove(annotating_feature.identity)
                for seq_identity in annotating_definition.sequences:
                    try:
                        annotating_doc.sequences.remove(seq_identity)
                    except ValueError:
                        pass

                self.logger.debug('Removed %s from %s', annotating_feature.identity, annotating_doc.name)

        self.logger.info('Finished cleaning up')

    def prune(self, target_library, cover_offset, min_target_length, ask_user=True, canonical_library=None,
              delete_flat=False, keep_flat=True, target_features=[], auto_swap=False, require_sequence=True):
        target_identities = set()
        for target_feature in target_features:
            target_identities.add(target_feature.identity)

        for target in target_library.features:
            if ((not require_sequence or self.__has_min_length(target, min_target_length))
                    and (len(target_identities) == 0 or target.identity in target_identities)):
                self.logger.info('Pruning %s', target.identity)

                target_doc = target_library.get_document(target.identity)

                target_definition = target_doc.getComponentDefinition(target.identity)

                cut_annos = [(sa.locations.getCut().at, sa.locations.getCut().at, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles), sa.description) for sa in target_definition.sequenceAnnotations if len(sa.locations) == 1 and sa.locations[0].getTypeURI() == sbol2.SBOL_CUT]
                annos = [(sa.locations.getRange().start, sa.locations.getRange().end, sa.identity, sa.displayId, sa.name, sa.component, set(sa.roles), sa.description) for sa in target_definition.sequenceAnnotations if len(sa.locations) == 1 and sa.locations[0].getTypeURI() == sbol2.SBOL_RANGE]
                
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

                if delete_flat:
                    for anno_group in grouped_annos:
                        flat_indices = self.__get_flat_annotation_indices(anno_group)

                        self.__remove_annotations(flat_indices, anno_group, target_definition)

                for anno_group in grouped_annos:
                    if len(anno_group) > 1:
                        selected_indices = self.__select_annotations(target_doc, target_definition, anno_group,
                            ask_user, canonical_library, keep_flat)

                        self.__remove_annotations(selected_indices, anno_group, target_definition)

                for anno_group in grouped_annos:
                    if ask_user or auto_swap:
                        if len(anno_group) == 2:
                            if (self.__are_swappable(anno_group[0], anno_group[1], target_definition) and
                                    (auto_swap or
                                    self.__ask_swap_annotations(anno_group[0], anno_group[1], target_definition))):
                                self.__swap_annotations(anno_group[0], anno_group[1], target_definition)
                            elif (self.__are_swappable(anno_group[1], anno_group[0], target_definition) and
                                    (auto_swap or
                                    self.__ask_swap_annotations(anno_group[1], anno_group[0], target_definition))):
                                self.__swap_annotations(anno_group[1], anno_group[0], target_definition)

                for anno_group in grouped_annos:
                    if len(anno_group) > 1:
                        redundant_URIs = []

                        for anno in anno_group:
                            if anno[5]:
                                redundant_URIs.append(target_definition.components.get(anno[5]).definition)
                            else:
                                redundant_URIs.append(anno[2])

                        self.logger.debug('Detected potentially redundant sub-parts in %s: %s.',
                                          target_definition.identity, 
                                          str(redundant_URIs))

                self.logger.info('Finished pruning %s', target.identity)

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    # Common arguments
    parser.add_argument('-n', '--namespace')
    parser.add_argument('-t', '--target_files', nargs='*', default=[])
    parser.add_argument('-o', '--output_files', nargs='*', default=[])
    parser.add_argument('-s', '--output_suffix', nargs='?', default='')
    parser.add_argument('-p', '--in_place', action='store_true')
    parser.add_argument('-m', '--min_target_length', nargs='?', default=2000)
    parser.add_argument('-mo', '--minimal_output', action='store_true')
    parser.add_argument('-ni', '--non_interactive', action='store_true')
    parser.add_argument('-l', '--log_file', nargs='?', default='')
    parser.add_argument('-v', '--validate', action='store_true')

    # Sequence annotation arguments
    parser.add_argument('-f', '--feature_files', nargs='*', default=[])
    parser.add_argument('-M', '--min_feature_length', nargs='?', default=40)
    
    parser.add_argument('-na', '--no_annotation', action='store_true')
    parser.add_argument('-e', '--extend_features', action='store_true')
    parser.add_argument('-xs', '--extension_suffix', nargs='?', default='')
    parser.add_argument('-x', '--extension_threshold', nargs='?', default=0.05)

    # Annotation pruning arguments
    parser.add_argument('-c', '--cover_offset', nargs='?', default=14)
    parser.add_argument('-r', '--deletion_roles', nargs='*', default=[])
    parser.add_argument('-d', '--delete_flat', action='store_true')
    parser.add_argument('-np', '--no_pruning', action='store_true')
    parser.add_argument('-a', '--auto_swap', action='store_true')
    
    # parser.add_argument('-s', '--sbh_URL', nargs='?', default=None)
    # parser.add_argument('-u', '--username', nargs='?', default=None)
    # parser.add_argument('-p', '--password', nargs='?', default=None)
    # parser.add_argument('-F', '--feature_URLs', nargs='*', default=[])
    # parser.add_argument('-T', '--target_URLs', nargs='*', default=[])
    # parser.add_argument('-o', '--sbh_output_file', nargs='?', default=None)
    
    args = parser.parse_args(args)

    logger = logging.getLogger('synbict')
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s ; %(levelname)s ; %(message)s')

    console_handler.setFormatter(formatter)

    logger.addHandler(console_handler)

    if len(args.log_file) > 0:
        file_handler = logging.FileHandler(args.log_file, "w")
        file_handler.setLevel(logging.DEBUG)

        file_handler.setFormatter(formatter)

        logger.addHandler(file_handler)

    sbol2.setHomespace(args.namespace)
    sbol2.Config.setOption('validate', args.validate)
    sbol2.Config.setOption('sbol_typed_uris', False)

    target_files = []
    for target_file in args.target_files:
        if os.path.isdir(target_file):
            target_files.extend([os.path.join(target_file, tf) for tf in os.listdir(target_file) if
                                 os.path.isfile(os.path.join(target_file, tf)) and (tf.endswith('.xml') or
                                                                                    tf.endswith('.sbol') or
                                                                                    tf.endswith('.gb') or
                                                                                    tf.endswith('.genbank') or
                                                                                    tf.endswith('.fasta') or
                                                                                    tf.endswith('.faa') or
                                                                                    tf.endswith('.fa') or
                                                                                    tf.endswith('.fas') or
                                                                                    tf.endswith('.fsa'))])
        else:
            target_files.append(target_file)

    output_files = []
    for i in range(0, len(target_files)):
        if len(args.output_files) == 1 and os.path.isdir(args.output_files[0]):
            (target_file_path, target_filename) = os.path.split(target_files[i])
            (target_file_base, target_file_extension) = os.path.splitext(target_filename)

            if len(args.output_suffix) > 0:
                output_files.append(os.path.join(args.output_files[0], '_'.join([target_file_base, args.output_suffix + '.xml'])))
            else:
                output_files.append(os.path.join(args.output_files[0], target_file_base + '.xml'))
        elif i < len(args.output_files):
            output_files.append(args.output_files[i])
        else:
            (target_file_base, target_file_extension) = os.path.splitext(target_files[i])

            if len(args.output_suffix) > 0:
                output_files.append('_'.join([target_file_base, args.output_suffix + '.xml']))
            else:
                output_files.append(target_file_base + '.xml')

    feature_docs = []
    for feature_file in args.feature_files:
        feature_docs.append(load_sbol(feature_file))

    feature_library = FeatureLibrary(feature_docs)

    if args.extend_features or not args.no_annotation:
        feature_annotater = FeatureAnnotater(feature_library, int(args.min_feature_length))

    if args.extend_features:
        target_docs = []
        for target_file in target_files:
            target_docs.append(load_target_file(target_file))
        for i in range(len(target_files) - 1, -1, -1):
            if not target_docs[i]:
                del target_docs[i]
                del target_files[i]
                del output_files[i]

        target_library = FeatureLibrary(target_docs, False)

        if args.minimal_output:
            output_docs = [sbol2.Document() for i in range(0, len(target_library.docs))]
        else:
            output_docs = []

        output_library = FeatureLibrary(output_docs, False)

        feature_curator = FeatureCurator(target_library, output_library)
        feature_curator.extend_features(feature_annotater,
                                        int(args.min_target_length),
                                        float(args.extension_threshold))

        for extended_doc in feature_annotater.get_updated_documents():
            if len(args.extension_suffix) > 0:
                (extended_file_base, extended_file_extension) = os.path.splitext(extended_doc.name)

                extended_file = '_'.join([extended_file_base, args.extension_suffix]) + '.xml'
            else:
                extended_file = extended_doc.name

            logger.info('Writing %s', extended_file)

            extended_doc.write(extended_file)

        if not args.no_annotation:
            (annotated_features, annotating_features) = feature_curator.annotate_features(feature_annotater,
                                                                                          int(args.min_target_length),
                                                                                          args.in_place)

            if args.minimal_output:
                for i in range(0, len(output_library.docs)):
                    if len(output_library.docs[i].componentDefinitions) == 0:
                        logger.warning('Failed to annotate %s, possibly no constructs found with minimum length %s',
                                        target_files[i], args.min_target_length)
            else:
                for i in target_library.get_non_updated_indices():
                    logger.warning('Failed to annotate %s, possibly no constructs found with minimum length %s',
                                    target_files[i], args.min_target_length)
        else:
            annotated_features = []
            annotating_features = []

        if not args.no_pruning:
            feature_pruner = FeaturePruner(feature_library, set(args.deletion_roles))
            feature_curator.prune_features(feature_pruner,
                                           int(args.cover_offset),
                                           int(args.min_target_length),
                                           annotated_features,
                                           annotating_features,
                                           args.delete_flat,
                                           args.auto_swap,
                                           not args.non_interactive)

        if not args.no_annotation or not args.no_pruning:
            if len(output_docs) > 0:
                for i in range(0, len(output_docs)):
                    if sbol2.Config.getOption('validate') == True:
                        logger.info('Validating and writing %s', output_files[i])
                    else:
                        logger.info('Writing %s', output_files[i])

                    output_docs[i].write(output_files[i])
            else:
                for i in range(0, len(target_docs)):
                    if sbol2.Config.getOption('validate') == True:
                        logger.info('Validating and writing %s', output_files[i])
                    else:
                        logger.info('Writing %s', output_files[i])

                    target_docs[i].write(output_files[i])
    else:
        for i in range(0, len(target_files)):
            target_doc = load_target_file(target_files[i])

            if target_doc:
                target_library = FeatureLibrary([target_doc], False)

                if args.minimal_output:
                    output_docs = [sbol2.Document()]
                else:
                    output_docs = []

                output_library = FeatureLibrary(output_docs, False)

                feature_curator = FeatureCurator(target_library, output_library)

                if args.no_annotation:
                    annotated_features = []
                    annotating_features = []
                else:
                    (annotated_features, annotating_features) = feature_curator.annotate_features(feature_annotater,
                                                                                                  int(args.min_target_length),
                                                                                                  args.in_place)

                    if args.minimal_output:
                        if len(output_library.docs[i].componentDefinitions) == 0:
                            logger.warning('Failed to annotate %s, possibly no constructs found with minimum length %s',
                                            target_files[i], args.min_target_length)
                    elif len(target_library.get_non_updated_indices()) > 0:
                        logger.warning('Failed to annotate %s, possibly no constructs found with minimum length %s',
                                        target_files[i], args.min_target_length)


                if not args.no_pruning:
                    feature_pruner = FeaturePruner(feature_library, set(args.deletion_roles))
                    feature_curator.prune_features(feature_pruner,
                                                   int(args.cover_offset),
                                                   int(args.min_target_length),
                                                   annotated_features,
                                                   annotating_features,
                                                   args.delete_flat,
                                                   args.auto_swap,
                                                   not args.non_interactive)

                if not args.no_annotation or not args.no_pruning:
                    if len(output_docs) == 1:
                        if sbol2.Config.getOption('validate') == True:
                            logger.info('Validating and writing %s', output_files[i])
                        else:
                            logger.info('Writing %s', output_files[i])

                        output_docs[0].write(output_files[i])
                    else:
                        if sbol2.Config.getOption('validate') == True:
                            logger.info('Validating and writing %s', output_files[i])
                        else:
                            logger.info('Writing %s', output_files[i])

                        target_doc.write(output_files[i])

    logger.info('Finished curating')

if __name__ == '__main__':
    main()
