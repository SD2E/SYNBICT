import math
import logging
import re
from .FeaturePruner import FeaturePruner
from .FeatureLibrary import FeatureLibrary

import sbol2
from Bio.Seq import Seq

# Set up the not found error for catching
try:
    NotFoundError = sbol2.SBOLError
except NameError:
    NotFoundError = RuntimeError

# Set up the not unique error for catching
try:
    NotUniqueError = sbol2.SBOLError
except NameError:
    NotUniqueError = RuntimeError

sbol2.Config.setOption('sbol_typed_uris', False)
sbol2.Config.setOption('sbol_compliant_uris', True)

# Compile regex patterns once at module level
_VARIANT_RE = re.compile(r'_v(\d+)')


class FeatureAnnotatorSimple:
    def __init__(self, feature_library, inline_matches, rc_matches):
        self.feature_library = feature_library
        self.logger = logging.getLogger('synbict')
        self.inline_matches = inline_matches
        self.rc_matches = rc_matches

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def __has_min_length(feature, min_feature_length):
        return min_feature_length == 0 or len(feature.nucleotides) >= min_feature_length

    @classmethod
    def __create_similar_sub_component(cls, parent_definition, child_definition):
        i = 1
        while True:
            try:
                tmp_id = _VARIANT_RE.sub(f'_v{i}', child_definition.displayId)
                sub_comp = parent_definition.components.create('_'.join([tmp_id, 'comp']))
            except (RuntimeError, NotUniqueError) as exc:
                if isinstance(exc, NotUniqueError) and exc.error_code() != sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    raise
                sub_comp = None

            if sub_comp is None:
                i += 1
                continue

            sub_comp.name = child_definition.name
            sub_comp.definition = child_definition.identity
            sub_comp.roleIntegration = None
            num = _VARIANT_RE.search(sub_comp.displayId).group(1)
            sub_comp.name = '_'.join(sub_comp.name.split('_')[:2] + [str(num)])
            return sub_comp

    @classmethod
    def __create_sub_component(cls, parent_definition, child_definition):
        i = 1
        while True:
            try:
                sub_comp = parent_definition.components.create(
                    '_'.join([child_definition.displayId, 'comp', str(i)])
                )
            except (RuntimeError, NotUniqueError) as exc:
                if isinstance(exc, NotUniqueError) and exc.error_code() != sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    raise
                sub_comp = None

            if sub_comp is None:
                i += 1
                continue

            sub_comp.name = child_definition.name
            sub_comp.definition = child_definition.identity
            sub_comp.roleIntegration = None
            return sub_comp

    @classmethod
    def __create_sequence_annotation(cls, parent_definition, child_definition, orientation,
                                     start, end, sub_comp_URI=None, parent_URI=None):
        i = 1
        while True:
            try:
                seq_anno = parent_definition.sequenceAnnotations.create(
                    '_'.join([child_definition.displayId, 'anno', str(i)])
                )
            except (RuntimeError, NotUniqueError) as exc:
                if isinstance(exc, NotUniqueError) and exc.error_code() != sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    raise
                seq_anno = None

            if seq_anno is None:
                i += 1
                continue

            seq_anno.name = child_definition.name
            seq_anno.description = child_definition.description
            if sub_comp_URI:
                seq_anno.component = sub_comp_URI
            if parent_URI:
                seq_anno.roles = seq_anno.roles + child_definition.roles
                seq_anno.wasDerivedFrom = seq_anno.wasDerivedFrom + [parent_URI]

            location = seq_anno.locations.createRange('_'.join([seq_anno.displayId, 'loc']))
            location.orientation = orientation
            location.start = start
            location.end = end
            return seq_anno

    @classmethod
    def __create_similar_sequence_annotation(cls, parent_definition, child_definition, orientation,
                                             start, end, sub_comp_URI=None, parent_URI=None):
        i = 1
        while True:
            try:
                seq_anno = parent_definition.sequenceAnnotations.create(
                    '_'.join([child_definition.displayId, 'anno', str(i)])
                )
            except (RuntimeError, NotUniqueError) as exc:
                if isinstance(exc, NotUniqueError) and exc.error_code() != sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    raise
                seq_anno = None

            if seq_anno is None:
                i += 1
                continue

            seq_anno.name = child_definition.name
            seq_anno.description = child_definition.description
            if sub_comp_URI:
                seq_anno.component = sub_comp_URI
            if parent_URI:
                seq_anno.roles = seq_anno.roles + child_definition.roles
                seq_anno.wasDerivedFrom = seq_anno.wasDerivedFrom + [parent_URI]

            location = seq_anno.locations.createRange('_'.join([seq_anno.displayId, 'loc']))
            location.orientation = orientation
            location.start = start
            location.end = end
            return seq_anno

    def create_component_definition(self, doc, base_display_id, name=None, roles=None,
                                    version=None, create_sequence=False,
                                    seq_elements=None, description=None):
        if roles is None:
            roles = []

        i = 1
        while True:
            display_id = base_display_id if i == 1 else f'{base_display_id}_{i}'
            try:
                cd = doc.componentDefinitions.create(display_id)
            except (RuntimeError, NotUniqueError) as exc:
                if isinstance(exc, NotUniqueError) and exc.error_code() != sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    raise
                cd = None

            if cd is None:
                i += 1
                continue

            cd.displayId = display_id
            if name is not None:
                cd.name = name
            cd.roles = list(roles)
            cd.description = description if description is not None else ''
            if version is not None:
                cd.version = str(version)

            if create_sequence:
                if seq_elements is None:
                    raise ValueError('seq_elements is required when create_sequence=True')
                j = 1
                while True:
                    seq_display_id = f'{display_id}_sequence' if j == 1 else f'{display_id}_sequence_{j}'
                    try:
                        seq = doc.sequences.create(seq_display_id)
                    except (RuntimeError, NotUniqueError) as exc:
                        if isinstance(exc, NotUniqueError) and exc.error_code() != sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                            raise
                        seq = None

                    if seq is None:
                        j += 1
                        continue

                    seq.elements = seq_elements
                    cd.sequences = [seq.identity]
                    break

            return cd

    # ------------------------------------------------------------------
    # Core processing
    # ------------------------------------------------------------------

    def _handle_unknown_feature(self, target_doc, target_definition, target_nucleotides,
                                orientation, start, end, match_type):
        """Create a hypothetical CDS or RNA component definition and annotate it."""
        if orientation != sbol2.SBOL_ORIENTATION_INLINE:
            target_nucleotides = str(Seq(target_nucleotides).reverse_complement()).upper()

        if match_type == 'CDS':
            new_compDef = self.create_component_definition(
                target_doc,
                base_display_id='hypothetical_cds',
                name='Hypothetical CDS',
                roles=[sbol2.SO_CDS],
                create_sequence=True,
                seq_elements=target_nucleotides,
                description='A hypothetical coding sequence',
            )
        else:
            new_compDef = self.create_component_definition(
                target_doc,
                base_display_id='hypothetical_rna',
                name='Hypothetical RNA',
                roles=['http://identifiers.org/so/SO:0001263'],
                create_sequence=True,
                seq_elements=target_nucleotides,
                description='A hypothetical RNA sequence',
            )

        sub_comp = self.__create_sub_component(target_definition, new_compDef)
        self.__create_sequence_annotation(target_definition, new_compDef, orientation, start, end, sub_comp.identity)

    def _handle_known_feature(self, target_doc, target_definition, feature, feature_definition,
                              feature_doc, orientation, start, end, identity_pct,
                              target_length, copy_definitions, complete_matches,
                              target_nucleotides):
        """Annotate a known feature — exact match or variant."""
        if end - start >= target_length and not complete_matches:
            return

        feature_role = FeaturePruner.get_common_role(feature_definition.roles)
        feature_ID = feature_definition.displayId if feature_definition.name is None else feature_definition.name

        rc_target_nucleotides = str(Seq(target_nucleotides).reverse_complement()).upper()

        feature_nucleotides = self.feature_library.get_feature(feature_definition.identity).nucleotides.upper()

        dna_exact = (target_nucleotides == feature_nucleotides or
                       rc_target_nucleotides == feature_nucleotides)
        # identity_pct == -1 means BWA DNA match — trust the aligner's inexact result
        # rather than re-checking nucleotides ourselves
        
        protein_match = math.isclose(float(identity_pct), 100.0, rel_tol=0.0, abs_tol=1e-6)
        exact_match = dna_exact and (identity_pct == -1 or protein_match)

        if exact_match:
            sub_comp = self.__create_sub_component(target_definition, feature_definition)
            self.__create_sequence_annotation(
                target_definition, feature_definition, orientation, start, end, sub_comp.identity
            )
            if copy_definitions:
                FeatureLibrary.copy_component_definition(feature_definition, feature_doc, target_doc)
        else:
            fid = str(feature_definition.identity)
            variant_definition = FeatureLibrary.copy_component_definition(
                feature_definition, feature_doc, target_doc, # <-- new, change to target_doc
                import_namespace=True, import_sequences=True,
                seq_elements=target_nucleotides,
                parent_definitions=[target_definition],
                parent_doc=target_doc, make_variant=True, strip_prefixes=[],
            )

            if variant_definition:
                variant_definition.wasDerivedFrom = [fid]
                num = _VARIANT_RE.search(variant_definition.displayId).group(1)
                variant_definition.name = f'{feature_definition.name}_variant_{num}'

                if identity_pct == -1:
                    variant_definition.description = (
                        f'Variant of part with the following description: {variant_definition.description}'
                    )
                elif protein_match:
                    variant_definition.description = (
                        f'Synonymous codon substitutions were introduced into the part described as follows: '
                        f'{variant_definition.description}. The resulting variant encodes a protein that is '
                        f'{identity_pct}% identical to the protein the original part coded for.'
                    )
                else:
                    variant_definition.description = (
                        f'Non-synonymous substitutions were introduced into the part described as follows: '
                        f'{variant_definition.description}. The resulting variant encodes a protein that is '
                        f'{identity_pct}% identical to the protein the original part coded for.'
                    )

                sub_comp = self.__create_similar_sub_component(target_definition, variant_definition)
                self.__create_similar_sequence_annotation(
                    target_definition, variant_definition, orientation, start, end, sub_comp.identity
                )
                #FeatureLibrary.copy_component_definition(variant_definition, feature_doc, target_doc) <-- new, deleted

        self.logger.debug(
            'Annotated %s (%s, %s) at [%s, %s] in %s',
            feature_definition.identity, feature_ID, feature_role,
            start, end, target_definition.identity,
        )

    def process_feature_matches(self, target_doc, target_definition, feature_matches,
                                orientation, target_length, copy_definitions=False,
                                complete_matches=False, output_matches=False):
        output_match_list = []
        seq_elements = target_definition.sequence.elements

        for feature_match in feature_matches:
            start = feature_match[1]
            end = feature_match[2]
            identity_pct = feature_match[3] if len(feature_match) > 3 else -1
            match_type = feature_match[4] if len(feature_match) > 3 else 'unknown'

            # Slice target nucleotides once per match
            target_nucleotides = seq_elements[start:end].upper()

            for feature in feature_match[0]:
                if feature.identity == ' ':
                    self._handle_unknown_feature(
                        target_doc, target_definition, target_nucleotides,
                        orientation, start, end, match_type,
                    )
                else:
                    feature_definition = self.feature_library.get_definition(feature.identity)
                    feature_doc = self.feature_library.get_document(feature.identity)

                    self._handle_known_feature(
                        target_doc, target_definition, feature,
                        feature_definition, feature_doc,
                        orientation, start, end, identity_pct,
                        target_length, copy_definitions, complete_matches,
                        target_nucleotides,
                    )

        # Moved outside the loop — only update once after all matches processed
        #self.feature_library.update()
        return output_match_list

    def annotate(self, inline_matches, rc_matches, target_library, min_target_length,
                 in_place=True, output_library=None, complete_matches=False,
                 strip_prefixes=[], output_matches=False, exact_match=True):
        annotated_identities = []
        output_match_lists = []

        for target in target_library.features:
            if not self.__has_min_length(target, min_target_length):
                continue

            self.logger.info('Annotating %s', target.identity)

            if not inline_matches and not rc_matches:
                self.logger.info('Finished annotating %s', target.identity)
                continue

            target_doc = target_library.get_document(target.identity)
            target_definition = target_doc.getComponentDefinition(target.identity)
            doc_index = target_library.get_document_index(target.identity)
            definition_copy = target_definition

            if not definition_copy:
                self.logger.warning(
                    '%s was not annotated because its version could not be incremented.',
                    target.identity,
                )
                continue

            copy_definitions = not in_place and (not output_library or doc_index >= len(output_library.docs))

            output_match_list = self.process_feature_matches(
                target_doc, definition_copy, inline_matches,
                sbol2.SBOL_ORIENTATION_INLINE,
                len(target.nucleotides), copy_definitions,
                complete_matches=complete_matches,
                output_matches=output_matches,
            )

            if rc_matches:
                output_match_list.extend(self.process_feature_matches(
                    target_doc, definition_copy, rc_matches,
                    sbol2.SBOL_ORIENTATION_REVERSE_COMPLEMENT,
                    len(target.nucleotides), copy_definitions,
                    complete_matches=complete_matches,
                    output_matches=output_matches,
                ))

            output_match_lists.append(output_match_list)
            annotated_identities.append(definition_copy.identity)

            self.logger.info('Finished annotating %s', target.identity)

        if output_matches:
            return annotated_identities, output_match_lists
        return annotated_identities
