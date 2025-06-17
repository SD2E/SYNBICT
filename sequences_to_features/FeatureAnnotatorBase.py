import logging
from sequences_to_features import FeaturePruner
from FeatureExtractor import Feature_No_Sequence
from sequences_to_features import FeatureLibrary
import argparse
import os
import sys
import requests
import json

from Bio.Seq import Seq
from Bio import Align
import sbol2

from sequences_to_features.sequences_to_features import NotUniqueError

# run this after alignment, input is inline_matches, output is sbol
class FeatureAnnotatorSimple:
    def __init__(self, inline_matches, rc_matches):
        self.logger = logging.getLogger('synbict')
        self.inline_matches = inline_matches
        self.rc_matches = rc_matches
        
    def _has_min_length(self, feature, min_length):
        #Check if feature has minimum required nucleotide length.
        return len(feature.nucleotides) >= min_length

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

    # need change this
    @classmethod
    def __process_feature_matches(self, target_doc, target_definition, feature_matches, orientation, target_length,
                                  copy_definitions=True, complete_matches=False):
        for feature_match in feature_matches:
            start = feature_match[2]
            end = feature_match[3]
            for feature in feature_match[0]: # is it list
                if end - start < target_length or complete_matches:
                    feature_definition = feature_match[0][0]

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

    # this function become the encapsulate insert the inline_matches or rc_matches to sbol
    # only test until this function, other functions is not inside the test
    def annotate(self, inline_matches,  rc_matches, target_library, min_target_length, in_place=False, output_library=None, complete_matches=False,
                 strip_prefixes=[]):
        annotated_identities = []

        for target in target_library.features:
            if self.__has_min_length(target, min_target_length):
                self.logger.info('Annotating %s', target.identity)

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
                                                                                       shallow_copy=True,
                                                                                       strip_prefixes=strip_prefixes)
                        else:
                            definition_copy = FeatureLibrary.copy_component_definition(target_definition,
                                                                                       target_doc,
                                                                                       output_doc, True,
                                                                                       min_target_length,
                                                                                       shallow_copy=True,
                                                                                       strip_prefixes=strip_prefixes)
                    elif in_place:
                        definition_copy = target_definition
                    else:
                        definition_copy = FeatureLibrary.copy_component_definition(target_definition, target_doc,
                                                                                   target_doc, True,
                                                                                   min_target_length,
                                                                                   strip_prefixes=strip_prefixes)
                    #self, target_doc, target_definition, feature_matches, orientation, target_length, copy_definitions=True, complete_matches=False
                    if definition_copy:
                        self.__process_feature_matches(target_doc, definition_copy, inline_matches,
                            len(target.nucleotides),
                            copy_definitions=(not output_library or doc_index >= len(output_library.docs)),
                            complete_matches=complete_matches)
                        self.__process_feature_matches(target_doc, definition_copy, rc_matches,
                            len(target.nucleotides), len(target.nucleotides) + 1,
                            (not output_library or doc_index >= len(output_library.docs)),
                            complete_matches=complete_matches)

                        annotated_identities.append(definition_copy.identity)
                    else:
                        self.logger.warning('%s was not annotated because its version could not be incremented.',
                                        target.identity)

                self.logger.info('Finished annotating %s', target.identity)

        return annotated_identities
