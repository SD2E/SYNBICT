import logging, re
from .FeaturePruner import FeaturePruner
from .FeatureLibrary import FeatureLibrary

import sbol2
from Bio.Seq import Seq
# run this after alignment, input is inline_matches, output is sbol
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
class FeatureAnnotatorSimple:
    def __init__(self, feature_library, inline_matches, rc_matches):
        self.feature_library = feature_library
        self.logger = logging.getLogger('synbict')
        self.inline_matches = inline_matches
        self.rc_matches = rc_matches
        
    @classmethod
    def __has_min_length(cls, feature, min_feature_length):
        return min_feature_length == 0 or len(feature.nucleotides) >= min_feature_length

    @classmethod
    def __create_similar_sub_component(cls, parent_definition, child_definition):
        i = 1
        while i > 0:
            try:
                tmp_id = re.sub(r"_v(\d+)", f"_v{i}", child_definition.displayId)
                sub_comp = parent_definition.components.create('_'.join([tmp_id,
                                                                         'comp']))
            except RuntimeError:
                sub_comp = None
            except NotUniqueError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    sub_comp = None
                else:
                    raise

            if sub_comp is None:
                match = re.search(r"_v(\d+)", child_definition.displayId)
                i = int(match.group(1)) + 1
            else:
                
                sub_comp.name = child_definition.name
                sub_comp.definition = child_definition.identity
                sub_comp.roleIntegration = None

                i = -1

        return sub_comp
    
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
    
    @classmethod
    def __create_similar_sequence_annotation(cls, parent_definition, child_definition, orientation, start, end,
                                     sub_comp_URI=None, parent_URI=None):
        i = 1

        while i > 0:
            try:
                tmp_id = re.sub(r"_v(\d+)", f"_v{i}", child_definition.displayId)
                seq_anno = parent_definition.sequenceAnnotations.create('_'.join([tmp_id,
                                                                                  'anno']))
                #print("created seq_anno with id: ", seq_anno.displayId)
            except RuntimeError:
                seq_anno = None
            except NotUniqueError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    seq_anno = None
                else:
                    raise

            if seq_anno is None:
                match = re.search(r"_v(\d+)", child_definition.displayId)
                i = int(match.group(1)) + 1
                
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
    
    # target_definition is 
    def process_feature_matches(self, target_doc, target_definition, feature_matches, orientation, target_length,
                                  copy_definitions=False, complete_matches=False, output_matches=False):
              
        output_match_list = []
        for feature_match in feature_matches:
            start = feature_match[1]
            end = feature_match[2]
            for feature in feature_match[0]: 
                if end - start < target_length or complete_matches:
                    feature_definition = self.feature_library.get_definition(feature.identity) 
                    if feature_definition.name is None:
                        feature_ID = feature_definition.displayId
                    else:
                        feature_ID = feature_definition.name
                    feature_doc = self.feature_library.get_document(feature.identity)
                    feature_role = FeaturePruner.get_common_role(feature_definition.roles)
                    target_nucleotides = target_definition.sequence.elements[start:end].upper()
                    rc_target_nucleotides = str(Seq(target_nucleotides).reverse_complement()).upper()
                    feature_seqs = FeatureLibrary.get_DNA_sequences(feature_definition, feature_doc)
                    feature_nucleotides = feature_seqs[0].elements.upper()
                    exact_match = (target_nucleotides == feature_nucleotides or rc_target_nucleotides == feature_nucleotides)
                    
                    if (exact_match):
                        sub_comp = self.__create_sub_component(target_definition, feature_definition)
                        self.__create_sequence_annotation(target_definition, feature_definition, orientation, start, end,
                                                        sub_comp.identity)
                        if copy_definitions: 
                            FeatureLibrary.copy_component_definition(feature_definition, feature_doc, target_doc)
                            
                    else:
                        # first copy: copy the variant componentDef into feature_doc, but the identity is wrong, which is duplicate of the other library, bug
                        # both compDef named "'http://seqimprove.synbiohub.org/AmpR_v1/1'"
                        variant_definition = FeatureLibrary.copy_component_definition(feature_definition,
                            feature_doc, feature_doc, import_namespace=True, import_sequences=True,
                            seq_elements=target_nucleotides, parent_definitions=[target_definition],
                            parent_doc=target_doc, make_variant=True, strip_prefixes=[])
                        #print("seq: ", variant_definition.sequence.elements) #sequence is correct, just doesn't create a new sbol:Sequence object
                        #because a duplicate identity happens when two Sequence object have duplicate identity (AmpR_v1/1)
                        
                        if variant_definition:
                            sub_identities = []
                            for sub_comp in variant_definition.components:
                                sub_identities.append(sub_comp.definition)
                        self.feature_library.update()
                        sub_comp = self.__create_similar_sub_component(target_definition, variant_definition)
                        #print("updated sub_comp id: ", sub_comp.displayId, sub_comp.identity)
                        self.__create_similar_sequence_annotation(target_definition, variant_definition, orientation, start, end,
                                                        sub_comp.identity)
                        FeatureLibrary.copy_component_definition(variant_definition, feature_doc, target_doc)

                    self.logger.debug('Annotated %s (%s, %s) at [%s, %s] in %s',
                                      feature_definition.identity,
                                      feature_ID,
                                      feature_role,
                                      start,
                                      end,
                                      target_definition.identity)
        return output_match_list 

    def annotate(self, inline_matches, rc_matches, target_library, min_target_length, in_place=True, output_library=None, complete_matches=False,
                 strip_prefixes=[], output_matches=False, exact_match=True):
        annotated_identities = []
        output_match_lists = []
        for target in target_library.features:
            if self.__has_min_length(target, min_target_length):
                self.logger.info('Annotating %s', target.identity)
                
                if len(inline_matches) > 0 or len(rc_matches) > 0:
                    target_doc = target_library.get_document(target.identity)

                    target_definition = target_doc.getComponentDefinition(target.identity)

                    doc_index = target_library.get_document_index(target.identity)
                    
                    definition_copy = target_definition

                    if definition_copy:
                        copy_definitions = (not output_library or doc_index >= len(output_library.docs))
                        if(in_place):
                            copy_definitions = False

                        output_match_list = self.process_feature_matches(target_doc,
                                                                           definition_copy,
                                                                           inline_matches,
                                                                           sbol2.SBOL_ORIENTATION_INLINE,
                                                                           len(target.nucleotides),
                                                                           copy_definitions,
                                                                           complete_matches=complete_matches,
                                                                           output_matches=output_matches)
                        if(len(rc_matches) > 0):
                            output_match_list.extend(self.process_feature_matches(target_doc,
                                                                                    definition_copy,
                                                                                    rc_matches,
                                                                                    sbol2.SBOL_ORIENTATION_REVERSE_COMPLEMENT,
                                                                                    len(target.nucleotides),
                                                                                    copy_definitions,
                                                                                    complete_matches=complete_matches,
                                                                                    output_matches=output_matches))
                        output_match_lists.append(output_match_list)

                        annotated_identities.append(definition_copy.identity)
                    else:
                        self.logger.warning('%s was not annotated because its version could not be incremented.',
                                        target.identity)

                self.logger.info('Finished annotating %s', target.identity)
        if output_matches:
            return annotated_identities, output_match_lists
        else:
            return annotated_identities
