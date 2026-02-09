import math
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
                num = re.search(r'_v(\d+)', sub_comp.displayId).group(1) # get the variant number
                out = "_".join(sub_comp.name.split('_')[:2] + [str(num)])
                sub_comp.name = out

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
    
    def create_component_definition(self, doc: sbol2.Document,
                                    base_display_id: str,
                                    name: str = None,
                                    roles=None,
                                    version: str = None,
                                    create_sequence: bool = False,
                                    seq_elements: str = None,
                                    description: str = None):
        """
        Create a new ComponentDefinition in `doc` with a unique displayId derived from base_display_id.

        Args:
            doc: SBOL2 Document to add the ComponentDefinition to
            base_display_id: base displayId to use (will append _cd_# if needed)
            name: human-readable name
            roles: list of role URIs (Sequence Ontology terms etc.)
            version: optional SBOL version string (depends on your SBOL homespace/versioning style)
            create_sequence: if True, also create/attach a Sequence
            seq_elements: DNA sequence string (required if create_sequence=True)
            description: optional description for the ComponentDefinition

        Returns:
            The created ComponentDefinition (sbol2.ComponentDefinition)
        """
        if roles is None:
            roles = []

        i = 1
        cd = None

        while i > 0:
            display_id = base_display_id if i == 1 else f"{base_display_id}_{i}"

            try:
                # Most pySBOL2 builds support doc.componentDefinitions.create(displayId)
                cd = doc.componentDefinitions.create(display_id)
            except RuntimeError:
                cd = None
            except NotUniqueError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    cd = None
                else:
                    raise

            if cd is None:
                i += 1
                continue

            # Fill metadata
            cd.displayId = display_id  # usually already set, but explicit is fine
            if name is not None:
                cd.name = name

            # Assign roles
            cd.roles = list(roles)
            
            # Assign description
            cd.description = description if description is not None else ""

            # Optional: set version if your objects use it
            if version is not None:
                cd.version = str(version)

            # Optional: create/attach a Sequence
            if create_sequence:
                if seq_elements is None:
                    raise ValueError("seq_elements is required when create_sequence=True")

                # Create a unique Sequence displayId too
                j = 1
                seq = None
                while j > 0:
                    seq_display_id = f"{display_id}_sequence" if j == 1 else f"{display_id}_sequence_{j}"
                    try:
                        seq = doc.sequences.create(seq_display_id)
                    except RuntimeError:
                        seq = None
                    except NotUniqueError as exc:
                        if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                            seq = None
                        else:
                            raise

                    if seq is None:
                        j += 1
                    else:
                        # Populate the sequence
                        seq.elements = seq_elements
                        # Attach
                        cd.sequences = [seq.identity]
                        j = -1

            i = -1

        return cd
    
    @classmethod
    def __create_similar_sequence_annotation(cls, parent_definition, child_definition, orientation, start, end,
                                     sub_comp_URI=None, parent_URI=None):
        i = 1

        while i > 0:
            try:
                tmp_id = re.sub(r"_v(\d+)", f"_v{i}", child_definition.displayId)
                seq_anno = parent_definition.sequenceAnnotations.create('_'.join([tmp_id,
                                                                                  'anno']))
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
            if len(feature_match) > 3:
                identity_pct = feature_match[3]
                type = feature_match[4]
            else:
                identity_pct = -1
                type = "unknown"
            
            for feature in feature_match[0]: 
                if(feature.identity == " "): 
                    # ref_name = " ", don't know the name; 1 hypothetical protein match; 2 hypothetical RNA match
                    #variant_definition.description = f"Hypothetical of part with the following description: {variant_definition.description}."
                    print("feature.identity: ", feature.identity)
                    print("feature type: ", type)
                    target_nucleotides = target_definition.sequence.elements[start:end].upper()
                    # assign DNA sequences from target instead from library, because it is not in the library
                    if(orientation != sbol2.SBOL_ORIENTATION_INLINE):
                        target_nucleotides = str(Seq(target_nucleotides).reverse_complement()).upper()
                    
                    if(type == "CDS"):
                        # create a new component definition for hypothetical protein or RNA
                        new_compDef = self.create_component_definition(target_doc,
                                base_display_id="hypothetical_cds",
                                name="Hypothetical CDS",
                                roles=[sbol2.SO_CDS],          # if you use SO terms in your codebase
                                create_sequence=True,
                                seq_elements=target_nucleotides,
                                description="A hypothetical coding sequence"
                            )
                    else:
                        new_compDef = self.create_component_definition(target_doc,
                                base_display_id="hypothetical_rna",
                                name="Hypothetical RNA",
                                roles=["http://identifiers.org/so/SO:0001263"],          # if you use SO terms in your codebase
                                create_sequence=True,
                                seq_elements=target_nucleotides,
                                description="A hypothetical RNA sequence"
                            )
                    
                    # Add the subcomponent and sequence annotation to the target definition
                    sub_comp = self.__create_sub_component(target_definition, new_compDef)
                    self.__create_sequence_annotation(target_definition, new_compDef, orientation, start, end, sub_comp.identity)
                    
                else:
                    if end - start < target_length or complete_matches:
                        feature_definition = self.feature_library.get_definition(feature.identity) 
                        feature_ID = feature_definition.displayId if feature_definition.name is None else feature_definition.name
                        
                        feature_doc = self.feature_library.get_document(feature.identity)
                        feature_role = FeaturePruner.get_common_role(feature_definition.roles)
                        target_nucleotides = target_definition.sequence.elements[start:end].upper()
                        rc_target_nucleotides = str(Seq(target_nucleotides).reverse_complement()).upper()
                        feature_seqs = FeatureLibrary.get_DNA_sequences(feature_definition, feature_doc)
                        feature_nucleotides = feature_seqs[0].elements.upper()
                        exact_match = (target_nucleotides == feature_nucleotides or rc_target_nucleotides == feature_nucleotides)
                        # translate into protein and compare
                        protein_match = math.isclose(float(identity_pct), 100.0, rel_tol=0.0, abs_tol=1e-6)
                        if (exact_match):
                            sub_comp = self.__create_sub_component(target_definition, feature_definition)
                            self.__create_sequence_annotation(target_definition, feature_definition, orientation, start, end,
                                                            sub_comp.identity)
                            if copy_definitions: 
                                FeatureLibrary.copy_component_definition(feature_definition, feature_doc, target_doc)                        
                        else:
                            # first copy: copy the variant componentDef into feature_doc, but the identity is wrong, which is duplicate of the other library, bug
                            # both compDef named "'http://seqimprove.synbiohub.org/AmpR_v1/1'"
                            
                            fid = str(feature_definition.identity)

                            variant_definition = FeatureLibrary.copy_component_definition(feature_definition,
                                feature_doc, feature_doc, import_namespace=True, import_sequences=True,
                                seq_elements=target_nucleotides, parent_definitions=[target_definition],
                                parent_doc=target_doc, make_variant=True, strip_prefixes=[])
                            print("variant_definition: ", variant_definition)
                            if variant_definition:
                                sub_identities = []
                                for sub_comp in variant_definition.components:
                                    print("sub_comp.definition: ", sub_comp.definition)
                                    sub_identities.append(sub_comp.definition)
                                
                                variant_definition.wasDerivedFrom = [fid]
                                num = re.search(r'_v(\d+)$', variant_definition.displayId).group(1)
                                variant_definition.name = f"{feature_definition.name}_variant_{num}"
                                
                                library_name = variant_definition.identity.split('/')[3]
                                sequence_id = variant_definition.sequence.identity
                                
                                if(identity_pct == -1): # DNA match
                                    variant_definition.description = f"Variant of part with the following description: {variant_definition.description}"
                                elif(protein_match): # exact protein match
                                    variant_definition.description = f"Synonymous codon substitutions were introduced into the part described as follows: {variant_definition.description}. The resulting variant encodes a protein that is {identity_pct}% identical to the protein the original part coded for."
                                else: # similar protein match
                                    variant_definition.description = f"Non-synonymous substitutions were introduced into the part described as follows: {variant_definition.description}. The resulting variant encodes a protein that is {identity_pct}% identical to the protein the original part coded for."
                                
                                self.feature_library.update()
                                sub_comp = self.__create_similar_sub_component(target_definition, variant_definition)
                                print("sub_comp: ", sub_comp)
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
