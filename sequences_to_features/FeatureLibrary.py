import sbol2
import logging
from .Feature import Feature
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
    def copy_sequence(cls, seq, source_doc, sink_doc, import_namespace=False, strip_prefixes=[]):
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

            cls.strip_origin_properties(seq_copy, strip_prefixes)
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

    # prop.startswith('http://wiki.synbiohub.org/wiki/Terms/synbiohub#')
    # prop.startswith('http://www.ncbi.nlm.nih.gov/genbank#')
    # prop.startswith('http://sbols.org/genBankConversion#')

    @classmethod
    def strip_origin_properties(cls, sbol_obj, other_prefixes):
        strip_props = []

        origin_prefixes = ['http://purl.org/dc/terms/created',
                           'http://purl.org/dc/terms/modified',
                           'http://purl.org/dc/terms/creator']

        strip_prefixes = tuple(origin_prefixes + other_prefixes)

        for prop in sbol_obj.properties:
            if prop.startswith(strip_prefixes):
                strip_props.append(prop)
        for strip_prop in strip_props:
            del sbol_obj.properties[strip_prop]

        sbol_obj.wasGeneratedBy = []

    @classmethod
    def copy_component_definition(cls, comp_definition, source_doc, sink_doc, import_namespace=False,
                                  min_seq_length=0, import_sequences=False, seq_elements=None,
                                  parent_definitions=[], parent_doc=None, make_variant=False,
                                  shallow_copy=False, strip_prefixes=[]):

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

                cls.strip_origin_properties(definition_copy, strip_prefixes)
                for sub_comp_copy in definition_copy.components:
                    cls.strip_origin_properties(sub_comp_copy, strip_prefixes)
                for anno_copy in definition_copy.sequenceAnnotations:
                    cls.strip_origin_properties(anno_copy, strip_prefixes)

                    for loc_copy in anno_copy.locations:
                        cls.strip_origin_properties(loc_copy, strip_prefixes)

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
                    seq_copy = cls.copy_sequence(seqs[0], source_doc, sink_doc, True, strip_prefixes)

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

                    cls.copy_sequence(seq, source_doc, sink_doc, False, strip_prefixes)

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
                                                                                shallow_copy=shallow_copy,
                                                                                strip_prefixes=strip_prefixes)

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