from importlib import simple
import logging
import argparse
import os
from pydoc import doc
import sys
import requests
import json

from Bio.Seq import Seq
from Bio import Align
import sbol2
from Feature import Feature
from FeatureLibrary import FeatureLibrary 
from FeaturePruner import FeaturePruner
from flashtext import KeywordProcessor
from Annotator import SAMFeatureMapper, TableFeatureMapper
from FeatureAnnotatorBase import FeatureAnnotatorSimple
from FeatureExtractor import FeatureExtractor
from BwaAligner import BwaAligner
from BlastAligner import BlastAligner
from Minimap2Aligner import Minimap2Aligner

# import time

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

    doc.name = non_sbol_file

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')

    return doc

class FeatureCurator():

    def __init__(self, target_library, output_library=None):
        self.target_library = target_library
        self.output_library = output_library

        self.logger = logging.getLogger('synbict')

    def annotate_features(self, feature_annotater, min_target_length, in_place=False, complete_matches=False,
                          strip_prefixes=[]):
        # start_time = time.clock()

        annotated_identities = feature_annotater.annotate(self.target_library, min_target_length, in_place,
                                                          self.output_library, complete_matches, strip_prefixes)

        # self.logger.info('Annotation Time: ' + str(time.clock() - start_time))

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

    def extend_features(self, feature_annotater, min_target_length, extension_threshold, strip_prefixes=[]):
        # start_time = time.clock()

        feature_annotater.extend_features_by_name(self.target_library,
                                                  min_target_length,
                                                  extension_threshold,
                                                  strip_prefixes)

        # self.logger.info('Extension Time: ' + str(time.clock() - start_time))

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
                                  rc_factor=0, copy_definitions=True, complete_matches=False, output_matches=False):
        output_match_list = []

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
                if len(feature.nucleotides) < target_length or complete_matches:
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
                    
                    if output_matches:
                        output_match_list.append({'feature_identity': feature_definition.identity,
                                                  'feature_ID':feature_ID,
                                                  'role': feature_role,
                                                  'start': start,
                                                  'end': end,
                                                  'orientation': orientation,
                                                  'target_identity': target_definition.identity})

                    self.logger.debug('Annotated %s (%s, %s) at [%s, %s] in %s',
                                      feature_definition.identity,
                                      feature_ID,
                                      feature_role,
                                      start,
                                      end,
                                      target_definition.identity)

        return output_match_list

    def extend_features_by_name(self, target_library, min_target_length, mismatch_threshold, strip_prefixes=[]):
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

                        target_nucleotides = target.nucleotides[anno_start - 1:anno_end].upper()
                        rc_target_nucleotides = str(Seq(target_nucleotides).reverse_complement()).upper()

                        inline_elements = ' '.join(target_nucleotides)
                        rc_elements = ' '.join(rc_target_nucleotides)

                        if not inline_elements in self.feature_matcher and not rc_elements in self.feature_matcher:
                            feature_definitions = self.feature_library.get_definitions_by_name(seq_anno.name)

                            for feature_definition in feature_definitions:
                                if set(seq_anno.roles) == set(feature_definition.roles):
                                    feature_doc = self.feature_library.get_document(feature_definition.identity)

                                    feature_seqs = FeatureLibrary.get_DNA_sequences(feature_definition, feature_doc)
                                    
                                    if len(feature_seqs) > 0:
                                        feature_nucleotides = feature_seqs[0].elements.upper()

                                        score = aligner.score(target_nucleotides, feature_nucleotides)
                                        rc_score = aligner.score(rc_target_nucleotides, feature_nucleotides)

                                        self.logger.debug('%s score %s', seq_anno.name, str(score))
                                        self.logger.debug('%s rc score %s', seq_anno.name, str(rc_score))
                                        self.logger.debug('target nucleotides: %s', target_nucleotides)
                                        self.logger.debug('feature nucleotides: %s', feature_nucleotides)

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
                                                parent_doc=target_doc, make_variant=True, strip_prefixes=strip_prefixes)

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

    def annotate_raw_sequences(self, raw_seqs, comp_IDs=[], min_target_length=0, complete_matches=False,
                               strip_prefixes=[]):
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

            self.annotate(target_library, min_target_length, True, complete_matches=complete_matches,
                          strip_prefixes=strip_prefixes)

        if len(annotated_comps) == 1:
            return annotated_comps[0]
        else:
            return annotated_comps

    def annotate(self, target_library, min_target_length, in_place=False, output_library=None, complete_matches=False,
                 strip_prefixes=[], output_matches=False):
        annotated_identities = []
        output_match_lists = []

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

                    if definition_copy:
                        copy_definitions = (not output_library or doc_index >= len(output_library.docs))

                        output_match_list = self.__process_feature_matches(target_doc,
                                                                           definition_copy,
                                                                           inline_matches,
                                                                           sbol2.SBOL_ORIENTATION_INLINE,
                                                                           len(target.nucleotides),
                                                                           copy_definitions=copy_definitions,
                                                                           complete_matches=complete_matches,
                                                                           output_matches=output_matches)
                        output_match_list.extend(self.__process_feature_matches(target_doc,
                                                                                definition_copy,
                                                                                rc_matches,
                                                                                sbol2.SBOL_ORIENTATION_REVERSE_COMPLEMENT,
                                                                                len(target.nucleotides),
                                                                                len(target.nucleotides) + 1,
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

def curate(feature_library, target_library, output_library, output_files, extend_features, no_annotation,
           min_feature_length, min_target_length, extension_threshold, extension_suffix, in_place, minimal_output,
           no_pruning, deletion_roles, cover_offset, delete_flat, auto_swap, non_interactive, logger,
           complete_matches=False, strip_prefixes=[], flashtext_mapping=True, bwa_mapping=False, minimap2_mapping=False,
           blastn_mapping=False, exact_match=False, build_index=False):
    if extend_features or not no_annotation:
        feature_annotater = FeatureAnnotater(feature_library, min_feature_length)
    
    feature_curator = FeatureCurator(target_library, output_library)

    if extend_features:
        feature_curator.extend_features(feature_annotater,
                                        min_target_length,
                                        extension_threshold,
                                        strip_prefixes)

        for extended_doc in feature_annotater.get_updated_documents():
            (extended_file_base, extended_file_extension) = os.path.splitext(extended_doc.name)

            if len(extension_suffix) > 0:
                extended_file = '_'.join([extended_file_base, extension_suffix]) + '.xml'
            else:
                extended_file = extended_file_base + '.xml'

            logger.info('Writing %s', extended_file)

            extended_doc.write(extended_file)

    if no_annotation:
        annotated_features = []
        annotating_features = []
    else:
        if(not flashtext_mapping):
            doc = target_library.docs[0] 
            index_prefix = 'test'
            if(bwa_mapping):
                bwa = BwaAligner(index_prefix)
                output_sam_path = 'aligned.sam'
                bwa.align(doc, output_sam_path, exact_match)
                mapper = SAMFeatureMapper('aligned.sam')
                inline_matches, rc_matches = mapper.extract_matches(min_feature_length, exact_match)

            elif(minimap2_mapping):
                minimap2 = Minimap2Aligner(index_prefix)
                output_sam_path = 'aligned.sam'
                minimap2.align(doc, output_sam_path, exact_match)
                mapper = SAMFeatureMapper('aligned.sam')
                inline_matches, rc_matches = mapper.extract_matches(min_feature_length, exact_match)

            elif(blastn_mapping):
                blast = BlastAligner(index_prefix)
                blast.align(doc, output_sam_path, exact_match) # True for exact match
                mapper = TableFeatureMapper('aligned.txt')
                inline_matches, rc_matches = mapper.extract_matches(min_feature_length, exact_match)

            simple = FeatureAnnotatorSimple(feature_library, inline_matches, rc_matches)
            # this is a different annnotate function, belong to FeatureAnnotatorSimple class
            simple.annotate(inline_matches, rc_matches, target_library, min_feature_length, in_place=True, output_library=output_library, output_matches=False)#True, in_place=True
        else:
            (annotated_features, annotating_features) = feature_curator.annotate_features(feature_annotater,
                                                                                        min_target_length,
                                                                                        in_place,
                                                                                        complete_matches,
                                                                                        strip_prefixes)

        if minimal_output:
            for i in range(0, len(output_library.docs)):
                if len(output_library.docs[i].componentDefinitions) == 0:
                    logger.warning('Failed to annotate %s, possibly no constructs found with minimum length %s',
                                   target_library.docs[i].name, min_target_length)
        else:
            for i in target_library.get_non_updated_indices():
                logger.warning('Failed to annotate %s, possibly no constructs found with minimum length %s',
                               target_library.docs[i].name, min_target_length)

    if not no_pruning:
        feature_pruner = FeaturePruner(feature_library, set(deletion_roles))

        feature_curator.prune_features(feature_pruner,
                                       cover_offset,
                                       min_target_length,
                                       annotated_features,
                                       annotating_features,
                                       delete_flat,
                                       auto_swap,
                                       not non_interactive)

    if not no_annotation or not no_pruning:
        if len(output_library.docs) > 0:
            for i in range(0, len(output_library.docs)):
                if sbol2.Config.getOption('validate') == True:
                    logger.info('Validating and writing %s', output_files[i])
                else:
                    logger.info('Writing %s', output_files[i])

                output_library.docs[i].write(output_files[i])
        else:
            for i in range(0, len(target_library.docs)):
                if sbol2.Config.getOption('validate') == True:
                    logger.info('Validating and writing %s', output_files[i])
                else:
                    logger.info('Writing %s', output_files[i])

                target_library.docs[i].write(output_files[i])

def download_sequences(doc, synbiohub):
    for comp_definition in doc.componentDefinitions:
        for seq_URI in comp_definition.sequences:
            download_sequence = False

            try:
                doc.getSequence(seq_URI)
            except RuntimeError:
                download_sequence = True
            except NotFoundError as exc:
                if is_sbol_not_found(exc):
                    download_sequence = True

            if download_sequence:
                try:
                    synbiohub.pull(seq_URI, doc)
                except NotFoundError as exc:
                    if is_sbol_not_found(exc):
                        logger.warning('Unable to download sequence %s, seq_URI)')
                    else:
                        raise
# run it the first time                  
def build_indexes(feature_docs):
    tmp = FeatureExtractor(feature_docs)
    fasta_path = 'test.fasta' #'/home/sophia/git_repo/SYNBICT/example/test.fasta'
    index_prefix = 'test'
    tmp.write_fasta(fasta_path)
    #tmp.write_metadata('test_metadata.json') # '/home/sophia/git_repo/SYNBICT/example/test_metadata.json'
    # build index for bwa
    tmp.build_index(fasta_path, index_prefix, 'bwa')
    # build index for blastn
    tmp.build_index(fasta_path, index_prefix, 'blast')
    # build index for minimap2
    tmp.build_index(fasta_path, index_prefix, 'minimap2')

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
    parser.add_argument('-cm', '--complete_matches', action='store_true')
    parser.add_argument('-sp', '--strip_prefixes', nargs='*', default=[])

    # Annotation pruning arguments
    parser.add_argument('-c', '--cover_offset', nargs='?', default=14)
    parser.add_argument('-r', '--deletion_roles', nargs='*', default=[])
    parser.add_argument('-d', '--delete_flat', action='store_true')
    parser.add_argument('-np', '--no_pruning', action='store_true')
    parser.add_argument('-a', '--auto_swap', action='store_true')
    
    parser.add_argument('-U', '--sbh_URL', nargs='?', default=None)
    parser.add_argument('-u', '--username', nargs='?', default=None)
    parser.add_argument('-w', '--password', nargs='?', default=None)
    parser.add_argument('-F', '--feature_URLs', nargs='*', default=[])
    parser.add_argument('-T', '--target_URLs', nargs='*', default=[])
    
    # Mapping method arguments
    parser.add_argument('-flashText', '--flashText_mapping', action='store_true')
    parser.add_argument('-bwa', '--bwa_mapping', action='store_true')
    parser.add_argument('-minimap2', '--minimap2_mapping', action='store_true')
    parser.add_argument('-blastn', '--blastn_mapping', action='store_true')
    parser.add_argument('-exact', '--exact_mapping', action='store_true')
    parser.add_argument('-bi', '--build_index', action='store_true')
    
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

    sbh_arg_types = []

    if args.username:
        sbh_arg_types.append('username')

    if args.password:
        sbh_arg_types.append('password')

    if args.sbh_URL:
        synbiohub = sbol2.PartShop(args.sbh_URL)

        if len(sbh_arg_types) == 2:
            try:
                synbiohub.login(args.username, args.password)
            except SBOLError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_BAD_HTTP_REQUEST:
                    logger.warning('Unable to log into SynBioHub instance with URL %s', args.sbh_URL)
        elif len(sbh_arg_types) == 1:
            logger.warning('SynBioHub %s was provided but %s is missing.',
                           sbh_arg_types[0],
                           {'username', 'password'}.difference(sbh_arg_types).pop())

    else:
        synbiohub = None

        if len(sbh_arg_types) > 0:
            logger.warning('SynBioHub %s were provided but a SynBioHub instance URL is missing.',
                           ', '.join(sbh_arg_types))

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

    for i in range(0, len(args.target_URLs)):
        if i + len(target_files) < len(args.output_files):
            output_files.append(args.output_files[i])
        else:
            output_files.append('_'.join(['output', str(i + len(target_files))]) + '.xml')

    feature_docs = []

    for feature_file in args.feature_files:
        feature_docs.append(load_sbol(feature_file))

    if synbiohub:
        for feature_URL in args.feature_URLs:
            feature_doc = sbol2.Document()

            try:
                synbiohub.pull(feature_URL, feature_doc)

                download_sequences(feature_doc, synbiohub)

                feature_docs.append(feature_doc)
            except NotFoundError as exc:
                if is_sbol_not_found(exc):
                    logger.warning('Unable to find feature URL %s at %s', feature_URL, sbh_URL)
                else:
                    raise

    feature_library = FeatureLibrary(feature_docs)
    # test here
    
    # build index for blastn, bwa, minimap2, this is pre-calculated for fast mode
    if args.build_index:
        build_indexes(feature_docs)
        logger.info('Finished building indexes')
        
    else:

        if args.extend_features:
            target_docs = []

            for target_file in target_files:
                target_docs.append(load_target_file(target_file))

            if synbiohub:
                for target_URL in args.target_URLs:
                    target_doc = sbol2.Document()

                    try:
                        synbiohub.pull(target_URL, target_doc)

                        download_sequences(target_doc, synbiohub)

                        target_docs.append(target_doc)
                    except NotFoundError as exc:
                        target_docs.append(None)

                        if is_sbol_not_found(exc):
                            logger.warning('Unable to find target URL %s at %s', target_URL, sbh_URL)
                        else:
                            raise

            filtered_output_files = [output_files[i] for i in range(0, len(target_docs)) if target_docs[i]]

            target_docs = [target_docs[i] for i in range(0, len(target_docs)) if target_docs[i]]

            target_library = FeatureLibrary(target_docs, False)

            if args.minimal_output:
                output_docs = [sbol2.Document() for i in range(0, len(target_library.docs))]
            else:
                output_docs = []

            output_library = FeatureLibrary(output_docs, False)

            curate(feature_library, target_library, output_library, filtered_output_files, args.extend_features,
                args.no_annotation, int(args.min_feature_length), int(args.min_target_length),
                float(args.extension_threshold), args.extension_suffix, args.in_place, args.minimal_output,
                args.no_pruning, args.deletion_roles, int(args.cover_offset), args.delete_flat, args.auto_swap,
                args.non_interactive, logger, args.complete_matches, args.strip_prefixes, args.flashText_mapping,
                        args.bwa_mapping, args.minimap2_mapping, args.blastn_mapping, args.exact_mapping, args.build_index)
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

                    curate(feature_library, target_library, output_library, [output_files[i]], args.extend_features,
                        args.no_annotation, int(args.min_feature_length), int(args.min_target_length),
                        float(args.extension_threshold), args.extension_suffix, args.in_place, args.minimal_output,
                        args.no_pruning, args.deletion_roles, int(args.cover_offset), args.delete_flat, args.auto_swap,
                        args.non_interactive, logger, args.complete_matches, args.strip_prefixes, args.flashText_mapping,
                        args.bwa_mapping, args.minimap2_mapping, args.blastn_mapping, args.exact_mapping, args.build_index)

            if synbiohub:
                for target_URL in args.target_URLs:
                    try:
                        target_doc = sbol2.Document()

                        synbiohub.pull(target_URL, target_doc)

                        download_sequences(target_doc, synbiohub)
                    except NotFoundError as exc:
                        if is_sbol_not_found(exc):
                            logger.warning('Unable to find target URL %s at %s', target_URL, sbh_URL)
                        else:
                            raise

                        target_doc = None

                    if target_doc:
                        target_library = FeatureLibrary([target_doc], False)

                        if args.minimal_output:
                            output_docs = [sbol2.Document()]
                        else:
                            output_docs = []

                        output_library = FeatureLibrary(output_docs, False)

                        curate(feature_library, target_library, output_library, [output_files[i]], args.extend_features,
                            args.no_annotation, int(args.min_feature_length), int(args.min_target_length),
                            float(args.extension_threshold), args.extension_suffix, args.in_place, args.minimal_output,
                            args.no_pruning, args.deletion_roles, int(args.cover_offset), args.delete_flat, args.auto_swap,
                            args.non_interactive, logger, args.complete_matches, args.strip_prefixes, args.flashText_mapping,
                        args.bwa_mapping, args.minimap2_mapping, args.blastn_mapping, args.exact_mapping, args.build_index)

        logger.info('Finished curating')

if __name__ == '__main__':
    main()
