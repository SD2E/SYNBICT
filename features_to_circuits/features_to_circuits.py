import logging
import argparse
import os
import sys

from Bio.Seq import Seq
from Bio import Align
import sbol
from flashtext import KeywordProcessor
from sequences_to_features import Feature
from sequences_to_features import FeatureLibrary

def load_sbol(sbol_file):
    logging.info('Loading %s', sbol_file)

    doc = sbol.Document()
    doc.read(sbol_file)

    doc.name = sbol_file

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')
    doc.addNamespace('http://sbolstandard.org/gff3#', 'gff3')
    doc.addNamespace('http://cellocad.org/Terms/cello#', 'cello')

    logging.info('Finished loading %s', sbol_file)

    return doc

class FeatureAnnotation():

    def __init__(self, identity, displayId, start, end, definition, roles):
        self.identity = identity
        self.displayId = displayId
        self.start = start
        self.end = end
        self.definition = definition
        self.roles = set(roles)

    def __lt__(self, other):
        return (self.start, self.end) < (other.start, other.end)

class Circuit():

    def __init__(self, identity, features):
        self.identity = identity
        self.features = features

    def is_covered(self, feature_library):
        for feature in self.features:
            if not feature_library.has_feature(feature.identity):
                return False

        return True

    def get_feature_identities(self):
        return [feature.identity for feature in self.features]

class Sensor():

    def __init__(self, identity, sensor_element, sensed_element):
        self.identity = identity
        self.sensor_element = sensor_element
        self.sensed_element = sensed_element

class CircuitLibrary():

    SBO_TEMPLATE = 'http://identifiers.org/biomodels.sbo/SBO:0000645'
    SBO_NON_COVALENT_BINDING = 'http://identifiers.org/biomodels.sbo/SBO:0000177'

    def __init__(self, docs, require_sequence=False):
        self.circuits = []
        self.sensors = []
        self.docs = docs
        self.__feature_library = FeatureLibrary(docs, require_sequence)
        self.__updated_indices = set()
        self.__circuit_map = {}
        self.__sensor_map = {}
        self.__template_to_product = {}
        self.__activator_to_dna= {}
        self.__repressor_to_dna = {}
        self.__dna_to_dna_repression = {}
        self.__dna_to_dna_activation = {}

        logging.info('Loading circuits')

        for i in range(0, len(self.docs)):
            for mod_definition in self.docs[i].moduleDefinitions:
                if len(mod_definition.interactions) == 1:
                    intxn = mod_definition.interactions[0]

                    features = self.__extract_features(self.docs[i], mod_definition)

                    if len(features) > 0:
                        if sbol.SBO_STIMULATION in intxn.types:
                            stimulator = None
                            stimulated = None

                            for par in intxn.participations:
                                if sbol.SBO_STIMULATOR in par.roles:
                                    stimulator = mod_definition.functionalComponents.get(par.participant).definition
                                elif sbol.SBO_STIMULATED in par.roles:
                                    stimulated = mod_definition.functionalComponents.get(par.participant).definition

                            if stimulator is not None and stimulated is not None:
                                if stimulator not in self.__activator_to_dna:
                                    self.__activator_to_dna[stimulator] = []

                                self.__activator_to_dna[stimulator].append(stimulated)
                        elif sbol.SBO_INHIBITION in intxn.types:
                            inhibitor = None
                            inhibited = None

                            for par in intxn.participations:
                                if sbol.SBO_INHIBITOR in par.roles:
                                    inhibitor = mod_definition.functionalComponents.get(par.participant).definition
                                elif sbol.SBO_INHIBITED in par.roles:
                                    inhibited = mod_definition.functionalComponents.get(par.participant).definition

                            if inhibitor is not None and inhibited is not None:
                                if inhibitor not in self.__repressor_to_dna:
                                    self.__repressor_to_dna[inhibitor] = []

                                self.__repressor_to_dna[inhibitor].append(inhibited)
                        elif sbol.SBO_GENETIC_PRODUCTION in intxn.types:
                            template = None
                            product = None

                            for par in intxn.participations:
                                if sbol.SBO_PRODUCT in par.roles:
                                    product = mod_definition.functionalComponents.get(par.participant).definition
                                elif self.SBO_TEMPLATE in par.roles:
                                    template = mod_definition.functionalComponents.get(par.participant).definition

                            if template is not None and product is not None:
                                if template not in self.__template_to_product:
                                    self.__template_to_product[template] = []

                                self.__template_to_product[template].append(product)

                        for template in self.__template_to_product:
                            for product in self.__template_to_product[template]:
                                if product in self.__activator_to_dna:
                                    self.__dna_to_dna_activation[template] = list(self.__activator_to_dna[product])

                                if product in self.__repressor_to_dna:
                                    self.__dna_to_dna_repression[template] = list(self.__repressor_to_dna[product])

                        self.circuits.append(Circuit(mod_definition.identity, features))

                        self.__circuit_map[mod_definition.identity] = i
                    else:
                        if self.SBO_NON_COVALENT_BINDING in intxn.types:
                            sensor_element = None
                            sensed_element = None

                            for par in intxn.participations:
                                if sbol.SBO_REACTANT in par.roles:
                                    fc = mod_definition.functionalComponents.get(par.participant)

                                    reactant_definition = self.docs[i].componentDefinitions.get(fc.definition)

                                    if sbol.BIOPAX_PROTEIN in reactant_definition.types:
                                        sensor_element = reactant_definition.identity
                                    elif sbol.BIOPAX_SMALL_MOLECULE in reactant_definition.types:
                                        sensed_element = reactant_definition.identity

                            if sensor_element and sensed_element:
                                self.sensors.append(Sensor(mod_definition.identity, sensor_element, sensed_element))

                                self.__sensor_map[mod_definition.identity] = i

        logging.info('Finished loading circuits')

    def __extract_features(self, doc, mod_definition):
        features = []

        for func_comp in mod_definition.functionalComponents:
            if self.__feature_library.has_feature(func_comp.definition):
                features.append(self.__feature_library.get_feature(func_comp.definition))
            # else:
            #     logging.warning('Cannot retrieve definition %s for feature in %s',
            #             func_comp.definition, mod_definition.identity)

        for sub_mod in mod_definition.modules:
            try:
                sub_definition = doc.moduleDefinitions.get(sub_mod.definition)
            except:
                sub_definition = None

            if sub_definition is not None:
                sub_features = cls.__extract_features(doc, sub_definition)

                if len(sub_features) == 0:
                    logging.warning('%s was not loaded since its sub-circuit %s is incomplete or contains no DNA features',
                        mod_definition.identity, sub_mod.definition)

                    return []
                else:
                    features.extend(sub_features)

        if len(features) == 0:
            logging.warning('%s was not loaded since it contains no DNA features', mod_definition.identity)

        return features

    def get_document(self, identity):
        return self.docs[self.get_document_index(identity)]

    def get_document_index(self, identity):
        if identity in self.__circuit_map:
            return self.__circuit_map[identity]
        else:
            return -1

    def get_definition(self, identity):
        return self.get_document(identity).moduleDefinitions.get(identity)

    def get_repressed_by_template(self, template):
        if template in self.__dna_to_dna_repression:
            return set(self.__dna_to_dna_repression[template])
        else:
            return set()

    def get_activated_by_template(self, template):
        if template in self.__dna_to_dna_activation:
            return set(self.__dna_to_dna_activation[template])
        else:
            return set()

    def extend_circuits_by_name(self, mismatch_threshold):
        logging.info('Extending circuit library')

        aligner = Align.PairwiseAligner()
        aligner.match_score = 1
        aligner.mismatch_score = -2
        aligner.internal_gap_score = -2.5

        circuit_feature_identities = set()

        for circuit in self.circuits:
            circuit_feature_identities.update(circuit.get_feature_identities())

        built_circuits = []
        built_circuit_indices = []
        for feature in self.__feature_library.features:
            if feature.identity not in circuit_feature_identities:
                feature_definition = self.__feature_library.get_definition(feature.identity)

                feature_nucleotides = feature.nucleotides
                rc_feature_nucleotides = str(Seq(feature_nucleotides).reverse_complement())

                for circuit in self.circuits:
                    if len(circuit.features) == 1:
                        circuit_feature_definition = self.__feature_library.get_definition(circuit.features[0].identity)

                        if feature_definition.name in circuit_feature_definition.name:
                            score = aligner.score(feature_nucleotides, circuit.features[0].nucleotides)
                            rc_score = aligner.score(rc_feature_nucleotides, circuit.features[0].nucleotides)

                            if rc_score > score:
                                best_score = rc_score
                                best_nucleotides = rc_feature_nucleotides
                            else:
                                best_score = score
                                best_nucleotides = feature_nucleotides

                            if len(best_nucleotides) < len(circuit.features[0].nucleotides):
                                max_score = len(best_nucleotides)
                            else:
                                max_score = len(circuit.features[0].nucleotides)

                            if max_score - best_score < mismatch_threshold*max_score:
                                circuit_doc = self.get_document(circuit.identity)
                                circuit_definition = self.get_definition(circuit.identity)

                                definition_copy = self.copy_module_definition(circuit_definition,
                                                                              circuit_doc,
                                                                              circuit_doc,
                                                                              True)

                                feature_doc = self.__feature_library.get_document(feature_definition.identity)

                                self.make_variant_circuit_definition(circuit_doc,
                                                                     feature_doc,
                                                                     definition_copy, 
                                                                     circuit_feature_definition,
                                                                     feature_definition)

                                built_circuits.append(Circuit(definition_copy.identity, [feature]))

                                built_circuit_indices.append(self.get_document_index(circuit.identity))

        for i in range(0, len(built_circuits)):
            self.circuits.append(built_circuits[i])

            self.__circuit_map[built_circuits[i].identity] = built_circuit_indices[i]

            self.__updated_indices.add(built_circuit_indices[i])

            logging.debug('Extended circuit library with %s', built_circuits[i].identity)

        logging.info('Finished extending circuit library')

    def get_updated_documents(self):
        updated_docs = []

        for updated_index in self.__updated_indices:
            updated_docs.append(self.docs[updated_index])

        return updated_docs

    @classmethod
    def reidentify_SBOL(cls, sbol_obj, target_path, replacement_path):
        sbol_obj.identity = sbol_obj.identity.replace(target_path, replacement_path)
        sbol_obj.persistentIdentity = sbol_obj.persistentIdentity.replace(target_path, replacement_path)

        try:
            sbol_obj.participant = sbol_obj.participant.replace(target_path, replacement_path)
        except AttributeError:
            pass

        return sbol_obj.identity

    @classmethod
    def make_variant_circuit_definition(cls, circuit_doc, variant_doc, circuit_definition, feature_definition,
            variant_definition):
        circuit_doc.moduleDefinitions.remove(circuit_definition.identity)

        variant_seqs = FeatureLibrary.get_DNA_sequences(variant_definition, variant_doc)
        variant_amino_acids = str(Seq(variant_seqs[0].elements, generic_dna).translate(stop_to=True)).upper()

        variant_index = 1
        unique_flag = False
        
        while not unique_flag:
            variant_ID = '_'.join([circuit_definition.displayId, 'v' + str(variant_index)])

            split_identity = circuit_definition.identity.split('/')
            variant_identity = '/'.join(split_identity[:-2] + [variant_ID, split_identity[-1]])
            variant_p_identity = '/'.join(split_identity[:-2] + [variant_ID])

            original_identity = circuit_definition.identity
            original_ID = circuit_definition.displayId
            original_p_identity = circuit_definition.persistentIdentity

            circuit_definition.identity = variant_identity
            circuit_definition.displayId = variant_ID
            circuit_definition.persistentIdentity = variant_p_identity

            try:
                circuit_doc.moduleDefinitions.add(circuit_definition)

                unique_flag = True
            except RuntimeError:
                circuit_definition.identity = original_identity
                circuit_definition.displayId = original_ID
                circuit_definition.persistentIdentity = original_p_identity

                variant_index = variant_index + 1

                unique_flag = False

        for sub_mod in circuit_definition.modules:
            cls.reidentify_SBOL(sub_mod, original_p_identity, circuit_definition.persistentIdentity)

        for fc in circuit_definition.functionalComponents:
            if fc.definition == feature_definition.identity:
                fc.definition = variant_definition.identity

            cls.reidentify_SBOL(fc, original_p_identity, circuit_definition.persistentIdentity)

        for intxn in circuit_definition.interactions:
            cls.reidentify_SBOL(intxn, original_p_identity, circuit_definition.persistentIdentity)

            for par in intxn.participations:
                cls.reidentify_SBOL(par, original_p_identity, circuit_definition.persistentIdentity)

            if sbol.SBO_GENETIC_PRODUCTION in intxn.types:
                for par in intxn.participations:
                    if sbol.SBO_PRODUCT in par.roles:
                        product_fc = circuit_definition.functionalComponents.get(par.participant)

                        product_definition = circuit_doc.componentDefinitions.get(product_fc.definition)

                        product_seqs = FeatureLibrary.get_sequences(product_definition, circuit_doc)

                        if len(product_seqs) > 0:
                            amino_acids = product_seqs[0].elements.upper()

                            if amino_acids != variant_amino_acids:
                                variant_product_definition = FeatureLibrary.copy_component_definition(product_definition,
                                    circuit_doc, circuit_doc, import_namespace=True, import_sequences=True,
                                    seq_elements=variant_amino_acids, make_variant=True)

                                product_fc.definition = variant_product_definition.identity

    @classmethod
    def strip_non_copy_properties(cls, sbol_obj):
        try:
            sbol_obj.setPropertyValue('http://wiki.synbiohub.org/wiki/Terms/synbiohub#ownedBy', '')
        except LookupError:
            pass

        try:
            sbol_obj.setPropertyValue('http://wiki.synbiohub.org/wiki/Terms/synbiohub#topLevel', '')
        except LookupError:
            pass

        try:
            sbol_obj.setPropertyValue('http://purl.org/dc/terms/created', '')
        except LookupError:
            pass

        sbol_obj.wasGeneratedBy = []

    @classmethod
    def copy_module_definition(cls, mod_definition, source_doc, sink_doc, import_namespace=False, deep_copy=False):
        if import_namespace:
            mod_namespace = '/'.join(mod_definition.identity.split('/')[:-2])

            copy_version = 1

            mod_definition_copy = None

            while not mod_definition_copy:
                try:
                    mod_definition_copy = mod_definition.copy(sink_doc, mod_namespace, str(copy_version))
                except RuntimeError:
                    mod_definition_copy = None

                    copy_version = copy_version + 1
        else:
            try:
                return sink_doc.moduleDefinitions.get(mod_definition.identity)
            except RuntimeError:
                mod_definition_copy = mod_definition.copy(sink_doc)

        cls.strip_non_copy_properties(mod_definition_copy)
        for sub_mod_copy in mod_definition_copy.modules:
            cls.strip_non_copy_properties(sub_mod_copy)
        for func_comp_copy in mod_definition_copy.functionalComponents:
            cls.strip_non_copy_properties(func_comp_copy)
        for intxn_copy in mod_definition_copy.interactions:
            cls.strip_non_copy_properties(intxn_copy)

            for parti_copy in intxn_copy.participations:
                cls.strip_non_copy_properties(parti_copy)

        if deep_copy:
            for sub_mod in mod_definition.modules:
                sub_mod_definition = source_doc.moduleDefinitions.get(sub_mod.definition)

                sub_mod_definition_copy = cls.copy_module_definition(sub_mod_definition, source_doc, sink_doc, import_namespace, deep_copy)

                sub_mod_copy = mod_definition_copy.modules.get(sub_mod.displayId)

                if sub_mod_definition_copy is not None:
                    sub_mod_copy.definition = sub_mod_definition_copy.identity

            for func_comp in mod_definition.functionalComponents:
                sub_comp_definition = source_doc.componentDefinitions.get(func_comp.definition)

                sub_comp_definition_copy = FeatureLibrary.copy_component_definition(sub_comp_definition, source_doc, sink_doc)

                sub_comp_copy = mod_definition_copy.components.get(func_comp.displayId)

                if sub_comp_definition_copy is not None:
                    sub_comp_copy.definition = sub_comp_definition_copy.identity
        else:
            for sub_mod in mod_definition.modules:
                sub_mod_copy = mod_definition_copy.modules.get(sub_mod.displayId)

                sub_mod_copy.definition = sub_mod.definition

            for func_comp in mod_definition.functionalComponents:
                func_comp_copy = mod_definition_copy.functionalComponents.get(func_comp.displayId)

                func_comp_copy.definition = func_comp.definition

        return mod_definition_copy

class CircuitBuilder():

    NCIT_BIOCHEMICAL_PATHWAY = 'http://purl.obolibrary.org/obo/NCIT_C20633'

    def __init__(self, circuit_library):
        self.circuit_library = circuit_library

    def add_sensors(self, target_doc, circuit_definition, sub_circuits, input_identities=set(),
            sensor_index=0, species_index=0):
        product_identities = set()

        for sub_circuit in sub_circuits:
            sub_circuit_definition = target_doc.moduleDefinitions.get(sub_circuit.identity)

            for intxn in sub_circuit_definition.interactions:
                if sbol.SBO_GENETIC_PRODUCTION in intxn.types:
                    for par in intxn.participations:
                        if sbol.SBO_PRODUCT in par.roles:
                            fc = sub_circuit_definition.functionalComponents.get(par.participant)

                            product_identities.add(fc.definition)

        covered_sensors = [
                            sensor 
                            for sensor in self.circuit_library.sensors
                            if sensor.sensor_element in product_identities]

        k = 0

        for i in range(0, len(covered_sensors)):
            sub_mod = circuit_definition.modules.create('sub_circuit_' + str(sensor_index + i + 1))

            sub_mod.definition = covered_sensors[i].identity

            sensor_doc = self.circuit_library.get_document(covered_sensors[i].identity)

            sensor_definition = self.circuit_library.get_definition(covered_sensors[i].identity)

            for sub_fc in sensor_definition.functionalComponents:
                if sub_fc.definition in input_identities:
                    species_fc = circuit_definition.functionalComponents.create('circuit_species_' + str(species_index + k + 1))

                    if sub_fc.name:
                        species_fc.name = sub_fc.name
                    species_fc.definition = sub_fc.definition
                    species_fc.direction = sbol.SBOL_DIRECTION_IN

                    k = k + 1

            CircuitLibrary.copy_module_definition(sensor_definition, sensor_doc, target_doc)

    @classmethod
    def infer_transcription(cls, target_doc, circuit_definition, constructs, tx_threshold):
        for construct in constructs:
            construct_definition = target_doc.componentDefinitions.get(construct.identity)

            inline_annos = []
            rc_annos = []

            for seq_anno in construct_definition.sequenceAnnotations:
                if (seq_anno.component
                        and len(seq_anno.locations) == 1
                        and seq_anno.locations[0].getTypeURI() == sbol.SBOL_RANGE):
                    anno_range = seq_anno.locations.getRange()

                    sub_comp = construct_definition.components.get(seq_anno.component)
                    part_definition = target_doc.componentDefinitions.get(sub_comp.definition)

                    feature_anno = FeatureAnnotation(sub_comp.identity,
                                                     sub_comp.displayId,
                                                     anno_range.start,
                                                     anno_range.end,
                                                     part_definition.identity,
                                                     part_definition.roles)

                    if anno_range.orientation == sbol.SBOL_ORIENTATION_INLINE:
                        inline_annos.append(feature_anno)
                    elif anno_range.orientation == sbol.SBOL_ORIENTATION_REVERSE_COMPLEMENT:
                        rc_annos.append(feature_anno)

            inline_annos.sort()
            rc_annos.sort()
            
            for i in range(0, len(inline_annos) - 1):
                for j in range(i, len(inline_annos)):
                    tx_distance = inline_annos[j].start - inline_annos[i].end

                    if tx_distance > tx_threshold:
                        j = len(inline_annos)
                    elif (tx_distance > 0 
                            and sbol.SO_PROMOTER in inline_annos[i].roles
                            and sbol.SO_CDS in inline_annos[j].roles):
                        stimulation = circuit_definition.interactions.create('_stimulates_'.join([inline_annos[i].displayId, inline_annos[j].displayId]))
                        stimulation.types = [sbol.SBO_STIMULATION]

                        try:
                            fc1 = circuit_definition.functionalComponents.create(rc_annos[i].displayId)
                            fc1.definition = rc_annos[i].definition
                        except:
                            fc1 = circuit_definition.functionalComponents.get(rc_annos[i].displayId)

                        try:
                            fc2 = circuit_definition.functionalComponents.create(rc_annos[j].displayId)
                            fc2.definition = rc_annos[j].definition
                        except:
                            fc2 = circuit_definition.functionalComponents.get(rc_annos[j].displayId)

                        stimulator = stimulation.participations.create(inline_annos[i].displayId)
                        stimulator.roles = [sbol.SBO_STIMULATOR]
                        stimulator.participant = fc1.identity

                        stimulated = stimulation.participations.create(inline_annos[j].displayId)
                        stimulated.roles = [sbol.SBO_STIMULATED]
                        stimulated.participant = fc2.identity

                        logging.debug('Inferred promoter-CDS stimulation %s', stimulation.identity)

            for i in range(0, len(rc_annos) - 1):
                for j in range(i, len(rc_annos)):
                    tx_distance = rc_annos[j].start - rc_annos[i].end

                    if tx_distance > tx_threshold:
                        j = len(rc_annos)
                    elif (tx_distance > 0
                            and sbol.SO_CDS in rc_annos[i].roles
                            and sbol.SO_PROMOTER in rc_annos[j].roles):
                        stimulation = circuit_definition.interactions.create('_stimulates_'.join([rc_annos[j].displayId, rc_annos[i].displayId]))
                        stimulation.types = [sbol.SBO_STIMULATION]

                        try:
                            fc1 = circuit_definition.functionalComponents.create(rc_annos[i].displayId)
                            fc1.definition = rc_annos[i].definition
                        except:
                            fc1 = circuit_definition.functionalComponents.get(rc_annos[i].displayId)

                        try:
                            fc2 = circuit_definition.functionalComponents.create(rc_annos[j].displayId)
                            fc2.definition = rc_annos[j].definition
                        except:
                            fc2 = circuit_definition.functionalComponents.get(rc_annos[j].displayId)

                        stimulator = stimulation.participations.create(rc_annos[j].displayId)
                        stimulator.roles = [sbol.SBO_STIMULATOR]
                        stimulator.participant = fc2.identity

                        stimulated = stimulation.participations.create(rc_annos[i].displayId)
                        stimulated.roles = [sbol.SBO_STIMULATED]
                        stimulated.participant = fc1.identity
                        
                        logging.debug('Inferred promoter-CDS stimulation %s', stimulation.identity)


    def build(self, circuit_id, target_doc, constructs, version='1', tx_threshold=0,
            input_identities=set(), output_identities=set()):
        circuit_definition = sbol.ModuleDefinition(circuit_id, version)
        circuit_definition.roles = [self.NCIT_BIOCHEMICAL_PATHWAY]

        logging.info('Building %s', circuit_definition.identity)

        if len(constructs) > 0:
            construct_feature_identities = set()

            for construct in constructs:
                for sub_identity in construct.sub_identities:
                    construct_feature_identities.add(sub_identity)

            circuit_feature_identities = set()

            covered_circuits = []

            for circuit in self.circuit_library.circuits:
                covered_features = []

                for feature in circuit.features:
                    if feature.identity in construct_feature_identities:
                        covered_features.append(feature)

                if len(circuit.features) == len(covered_features):
                    covered_circuits.append(circuit)

                    for covered_feature in covered_features:
                        circuit_feature_identities.add(covered_feature.identity)

            if len(covered_circuits) > 0:
                covered_constructs = []

                for construct in constructs:
                    covered_count = 0

                    for sub_identity in construct.sub_identities:
                        if sub_identity in circuit_feature_identities:
                            covered_count = covered_count + 1

                    if covered_count > 0:
                        covered_constructs.append(construct)
                    else:
                        logging.warning('No sub-circuits found for construct %s', construct.identity)
            
                for i in range(0, len(covered_constructs)):
                    func_comp = circuit_definition.functionalComponents.create('construct_' + str(i + 1))

                    func_comp.definition = covered_constructs[i].identity

                k = 0

                for i in range(0, len(covered_circuits)):
                    sub_mod = circuit_definition.modules.create('sub_circuit_' + str(i + 1))

                    sub_mod.definition = covered_circuits[i].identity

                    sub_circuit_doc = self.circuit_library.get_document(covered_circuits[i].identity)

                    sub_circuit_definition = self.circuit_library.get_definition(covered_circuits[i].identity)

                    for sub_fc in sub_circuit_definition.functionalComponents:
                        if sub_fc.definition in input_identities:
                            species_fc = circuit_definition.functionalComponents.create('species_' + str(k + 1))

                            if sub_fc.name:
                                species_fc.name = sub_fc.name
                            species_fc.definition = sub_fc.definition
                            species_fc.direction = sbol.SBOL_DIRECTION_IN

                            k = k + 1
                        elif sub_fc.definition in output_identities:
                            species_fc = circuit_definition.functionalComponents.create('circuit_species_' + str(k + 1))
                            
                            if sub_fc.name:
                                species_fc.name = sub_fc.name
                            species_fc.definition = sub_fc.definition
                            species_fc.direction = sbol.SBOL_DIRECTION_OUT

                            k = k + 1

                    CircuitLibrary.copy_module_definition(sub_circuit_definition, sub_circuit_doc, target_doc)

                    logging.debug('Added sub-circuit %s', sub_circuit_definition.identity)

                self.add_sensors(target_doc, circuit_definition, covered_circuits, input_identities,
                    len(covered_circuits), k)

                self.infer_transcription(target_doc, circuit_definition, constructs, tx_threshold)

                target_doc.addModuleDefinition(circuit_definition)

                logging.info('Finished building %s', circuit_definition.identity)

                return True
            else:
                logging.error('Failed to build %s, no sub-circuits found for constructs', circuit_definition.identity)

                return False
        else:
            logging.error('Failed to build %s, no constructs found with minimum length', circuit_definition.identity)

            return False

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--namespace')
    parser.add_argument('-t', '--target_files', nargs='*', default=[])
    parser.add_argument('-c', '--sub_circuit_files', nargs='+')
    parser.add_argument('-i', '--circuit_IDs', nargs='*', default=[])
    parser.add_argument('-o', '--output_files', nargs='*', default=[])
    parser.add_argument('-m', '--min_target_length', nargs='?', default='2000')
    parser.add_argument('-l', '--log_file', nargs='?', default='')
    parser.add_argument('-cv', '--circuit_version', nargs='?', default='1')
    parser.add_argument('-v', '--validate', action='store_true')
    parser.add_argument('-s', '--circuit_suffix', nargs='?', default='')
    parser.add_argument('-e', '--extend_sub_circuits', action='store_true')
    parser.add_argument('-x', '--extension_threshold', nargs='?', default='0.05')
    parser.add_argument('-d', '--tx_threshold', nargs='?', default='200')
    parser.add_argument('-nb', '--no_build', action='store_true')
    parser.add_argument('-ii', '--input_identities', nargs='*', default=[])
    parser.add_argument('-oi', '--output_identities', nargs='*', default=[])
    
    args = parser.parse_args(args)

    if len(args.log_file) > 0:
        logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s ; %(levelname)s ; %(message)s',
                        datefmt='%m-%d-%y %H:%M',
                        filename=args.log_file,
                        filemode='w')

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)

    console_formatter = logging.Formatter('%(levelname)s ; %(message)s')
    
    console_handler.setFormatter(console_formatter)
    
    logging.getLogger('').addHandler(console_handler)

    sbol.setHomespace(args.namespace)
    sbol.Config.setOption('validate', args.validate)
    sbol.Config.setOption('sbol_typed_uris', False)

    circuit_docs = []
    for circuit_file in args.sub_circuit_files:
        circuit_docs.append(load_sbol(circuit_file))

    if args.extend_sub_circuits:
        circuit_library = CircuitLibrary(circuit_docs, True)

        circuit_library.extend_circuits_by_name(float(args.extension_threshold))

        for extended_doc in circuit_library.get_updated_documents():
            (extended_file_base, extended_file_extension) = os.path.splitext(extended_doc.name)
            extended_file = extended_file_base + '_extended.xml'

            logging.info('Writing %s', extended_file)

            extended_doc.write(extended_file)
    else:
        circuit_library = CircuitLibrary(circuit_docs)

    if not args.no_build:
        circuit_builder = CircuitBuilder(circuit_library)

        circuit_memo = set()
        circuit_index = 0

        target_files = []
        for target_file in args.target_files:
            if os.path.isdir(target_file):
                target_files.extend([os.path.join(target_file, tf) for tf in os.listdir(target_file) if
                                     os.path.isfile(os.path.join(target_file, tf))])
            else:
                target_files.append(target_file)

        for i in range (0, len(target_files)):
            target_doc = load_sbol(target_files[i])
            
            target_library = FeatureLibrary([target_doc], False)

            if i < len(args.circuit_IDs):
                circuit_ID = args.circuit_IDs[i]
            else:
                (circuit_ID, file_extension) = os.path.splitext(os.path.basename(target_files[i]))

                if len(args.circuit_suffix) > 0:
                    circuit_ID = '_'.join([circuit_ID, args.circuit_suffix])
                else:
                    circuit_ID + '_circuit'

            unique_ID = circuit_ID

            while unique_ID in circuit_memo:
                circuit_index = circuit_index + 1

                unique_ID = circuit_ID + str(circuit_index)
            
            circuit_memo.add(unique_ID)

            build_success = circuit_builder.build(
                                                    circuit_ID,
                                                    target_doc,
                                                    target_library.get_features(int(args.min_target_length), True),
                                                    args.circuit_version,
                                                    int(args.tx_threshold),
                                                    set(args.input_identities),
                                                    set(args.output_identities))

            if build_success:
                if len(args.output_files) == 1 and os.path.isdir(args.output_files[0]):
                    (target_file_path, target_filename) = os.path.split(target_files[i])
                    (target_file_base, target_file_extension) = os.path.splitext(target_filename)

                    if len(args.circuit_suffix) > 0:
                        output_file = os.path.join(args.output_files[0], '_'.join([target_file_base, args.circuit_suffix + target_file_extension]))
                    else:
                        output_file = os.path.join(args.output_files[0], target_file_base + target_file_extension)

                elif i < len(args.output_files):
                    output_file = args.output_files[i]
                else:
                    (target_file_base, target_file_extension) = os.path.splitext(target_files[i])

                    if len(args.circuit_suffix) > 0:
                        output_file = '_'.join([target_file_base, args.circuit_suffix + target_file_extension])
                    else:
                        output_file = target_file_base + target_file_extension

                if sbol.Config.getOption('validate') == True:
                    logging.info('Validating and writing %s', output_file)
                else:
                    logging.info('Writing %s', output_file)

                target_doc.write(output_file)

                logging.info('Finished writing %s', output_file)

    logging.info('Finished curating')


class CircuitIdentityError(Exception):

    def __init__(self, circuit_identity, feature_ID):
        self.circuit_identity = circuit_identity
        self.feature_ID = feature_ID

    def __str__(self):
        return "Circuit identity {ci} does not contain the feature {fi}.".format(ci=self.circuit_identity, fi=self.feature_ID)


if __name__ == '__main__':
    main()
