import logging
import argparse
import os
import sys

from Bio.Seq import Seq
from Bio import Align
import sbol2
from flashtext import KeywordProcessor
from sequences_to_features import Feature
from sequences_to_features import FeatureLibrary

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

    logger.info('Finished loading %s', sbol_file)

    return doc

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

class FeatureAnnotation():

    def __init__(self, start, end, definition=None, roles=[], display_ID=None,
                 name=None, comp_name=None, comp_ID=None, orientation=None,
                 comp_identity=None, identity=None):
        self.start = start
        self.end = end
        self.definition = definition
        self.roles = set(roles)
        self.display_ID = display_ID
        self.name = name
        self.comp_ID = comp_ID
        self.comp_name = comp_name
        self.orientation = orientation
        self.comp_identity = comp_identity
        self.identity = identity

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

    def __init__(self, identity, sensor_element, sensed_element, sensor_type):
        self.identity = identity
        self.sensor_element = sensor_element
        self.sensed_element = sensed_element
        self.sensor_types = {sensor_type}

    def add_sensor_type(self, sensor_type):
        self.sensor_types.add(sensor_type)

    def get_abstraction_types(self):
        return list(self.sensor_types.difference({CircuitLibrary.SBO_NON_COVALENT_BINDING}))

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
        self.logger = logging.getLogger('synbict')

        self.logger.info('Loading circuits')

        for i in range(0, len(self.docs)):
            for mod_definition in self.docs[i].moduleDefinitions:
                # if len(mod_definition.interactions) == 1:
                features = self.__extract_features(self.docs[i], mod_definition)

                if len(features) > 0:
                    for intxn in mod_definition.interactions:
                        # intxn = mod_definition.interactions[0]

                        if sbol2.SBO_STIMULATION in intxn.types:
                            stimulator = None
                            stimulated = None

                            for par in intxn.participations:
                                if sbol2.SBO_STIMULATOR in par.roles:
                                    stimulator = mod_definition.functionalComponents.get(par.participant).definition
                                elif sbol2.SBO_STIMULATED in par.roles:
                                    stimulated = mod_definition.functionalComponents.get(par.participant).definition

                            if stimulator and stimulated:
                                if stimulator not in self.__activator_to_dna:
                                    self.__activator_to_dna[stimulator] = []

                                self.__activator_to_dna[stimulator].append(stimulated)
                        elif sbol2.SBO_INHIBITION in intxn.types:
                            inhibitor = None
                            inhibited = None

                            for par in intxn.participations:
                                if sbol2.SBO_INHIBITOR in par.roles:
                                    inhibitor = mod_definition.functionalComponents.get(par.participant).definition
                                elif sbol2.SBO_INHIBITED in par.roles:
                                    inhibited = mod_definition.functionalComponents.get(par.participant).definition

                            if inhibitor and inhibited:
                                if inhibitor not in self.__repressor_to_dna:
                                    self.__repressor_to_dna[inhibitor] = []

                                self.__repressor_to_dna[inhibitor].append(inhibited)
                        elif sbol2.SBO_GENETIC_PRODUCTION in intxn.types:
                            template = None
                            product = None

                            for par in intxn.participations:
                                if sbol2.SBO_PRODUCT in par.roles:
                                    product = mod_definition.functionalComponents.get(par.participant).definition
                                elif self.SBO_TEMPLATE in par.roles:
                                    template = mod_definition.functionalComponents.get(par.participant).definition

                            if template and product:
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
                    for intxn in mod_definition.interactions:
                        if self.SBO_NON_COVALENT_BINDING in intxn.types:
                            sensor_element = None
                            sensed_element = None

                            for par in intxn.participations:
                                if sbol2.SBO_REACTANT in par.roles:
                                    reactant = mod_definition.functionalComponents.get(par.participant).definition

                                    reactant_definition = self.docs[i].componentDefinitions.get(reactant)

                                    if sbol2.BIOPAX_PROTEIN in reactant_definition.types:
                                        sensor_element = reactant_definition.identity
                                    elif sbol2.BIOPAX_SMALL_MOLECULE in reactant_definition.types:
                                        sensed_element = reactant_definition.identity

                            if sensor_element and sensed_element:
                                self.sensors.append(Sensor(mod_definition.identity,
                                                           sensor_element,
                                                           sensed_element,
                                                           self.SBO_NON_COVALENT_BINDING))

                                self.__sensor_map[mod_definition.identity] = i

        self.__abstract_sensors()

        self.logger.info('Finished loading circuits')

    def __abstract_sensors(self):
        for sensor in self.sensors:
            if self.SBO_NON_COVALENT_BINDING in sensor.sensor_types:
                if sensor.sensor_element in self.__activator_to_dna:
                    sensor.add_sensor_type(sbol2.SBO_STIMULATION)
                elif sensor.sensor_element in self.__repressor_to_dna:
                    sensor.add_sensor_type(sbol2.SBO_INHIBITION)

    def __extract_features(self, doc, mod_definition):
        features = []

        for func_comp in mod_definition.functionalComponents:
            if self.__feature_library.has_feature(func_comp.definition):
                features.append(self.__feature_library.get_feature(func_comp.definition))
            # else:
            #     self.logger.warning('Cannot retrieve definition %s for feature in %s',
            #             func_comp.definition, mod_definition.identity)

        for sub_mod in mod_definition.modules:
            try:
                sub_definition = doc.moduleDefinitions.get(sub_mod.definition)
            except:
                sub_definition = None

            if sub_definition:
                sub_features = cls.__extract_features(doc, sub_definition)

                if len(sub_features) == 0:
                    self.logger.warning('%s was not loaded since its sub-circuit %s is incomplete or contains no DNA features',
                        mod_definition.identity, sub_mod.definition)

                    return []
                else:
                    features.extend(sub_features)

        if len(features) == 0:
            self.logger.warning('%s was not loaded since it contains no DNA features', mod_definition.identity)

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

    def extend_circuits_by_name(self, mismatch_threshold, strip_prefixes=[]):
        self.logger.info('Extending circuit library')

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

                            # if len(best_nucleotides) < len(circuit.features[0].nucleotides):
                            #     max_score = len(best_nucleotides)
                            # else:
                            #     max_score = len(circuit.features[0].nucleotides)

                            if len(feature_nucleotides) < len(circuit.features[0].nucleotides):
                                max_score = len(feature_nucleotides)
                            else:
                                max_score = len(circuit.features[0].nucleotides)

                            if max_score - best_score < mismatch_threshold*max_score:
                                circuit_doc = self.get_document(circuit.identity)
                                circuit_definition = self.get_definition(circuit.identity)

                                definition_copy = self.copy_module_definition(circuit_definition,
                                                                              circuit_doc,
                                                                              circuit_doc,
                                                                              True,
                                                                              strip_prefixes=strip_prefixes)

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

            self.logger.debug('Extended circuit library with %s', built_circuits[i].identity)

        self.logger.info('Finished extending circuit library')

    def get_updated_documents(self):
        updated_docs = []

        for updated_index in self.__updated_indices:
            updated_docs.append(self.docs[updated_index])

        return updated_docs

    @classmethod
    def reidentify_SBOL(cls, sbol_obj, target_path, replacement_path):
        sbol_obj.identity = sbol_obj.identity.replace(target_path, replacement_path)
        split_identity = sbol_obj.identity.split('/')
        sbol_obj.identity = '/'.join(split_identity[:-1] + ['1'])

        sbol_obj.version = '1'

        sbol_obj.persistentIdentity = sbol_obj.persistentIdentity.replace(target_path, replacement_path)

        try:
            sbol_obj.participant = sbol_obj.participant.replace(target_path, replacement_path)
            split_par_identity = sbol_obj.participant.split('/')
            sbol_obj.participant = '/'.join(split_par_identity[:-1] + ['1'])
        except AttributeError:
            pass

        return sbol_obj.identity

    @classmethod
    def make_variant_circuit_definition(cls, circuit_doc, variant_doc, circuit_definition, feature_definition,
            variant_definition):
        circuit_doc.moduleDefinitions.remove(circuit_definition.identity)

        variant_seqs = FeatureLibrary.get_DNA_sequences(variant_definition, variant_doc)
        variant_amino_acids = str(Seq(variant_seqs[0].elements).translate(to_stop=True)).upper()

        variant_index = 1
        unique_flag = False
        
        while not unique_flag:
            variant_ID = '_'.join([circuit_definition.displayId, 'v' + str(variant_index)])

            split_identity = circuit_definition.identity.split('/')
            variant_identity = '/'.join(split_identity[:-2] + [variant_ID, '1'])
            variant_p_identity = '/'.join(split_identity[:-2] + [variant_ID])

            original_identity = circuit_definition.identity
            original_ID = circuit_definition.displayId
            original_p_identity = circuit_definition.persistentIdentity
            original_version = circuit_definition.version

            circuit_definition.identity = variant_identity
            circuit_definition.displayId = variant_ID
            circuit_definition.persistentIdentity = variant_p_identity
            circuit_definition.version = '1'

            try:
                circuit_doc.moduleDefinitions.add(circuit_definition)

                unique_flag = True
            except RuntimeError:
                circuit_definition.identity = original_identity
                circuit_definition.displayId = original_ID
                circuit_definition.persistentIdentity = original_p_identity
                circuit_definition.version = original_version

                variant_index = variant_index + 1

                unique_flag = False
            except NotUniqueError as exc:
                if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                    circuit_definition.identity = original_identity
                    circuit_definition.displayId = original_ID
                    circuit_definition.persistentIdentity = original_p_identity
                    circuit_definition.version = original_version

                    variant_index = variant_index + 1

                    unique_flag = False
                else:
                    raise

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

            if sbol2.SBO_GENETIC_PRODUCTION in intxn.types:
                for par in intxn.participations:
                    if sbol2.SBO_PRODUCT in par.roles:
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

    # @classmethod
    # def strip_origin_properties(cls, sbol_obj):
    #     origin_props = []

    #     for prop in sbol_obj.properties:
    #         if (prop.startswith('http://wiki.synbiohub.org/wiki/Terms/synbiohub#')
    #                 or prop.startswith('http://purl.org/dc/terms/created')
    #                 or prop.startswith('http://purl.org/dc/terms/modified')
    #                 or prop.startswith('http://www.ncbi.nlm.nih.gov/genbank#')
    #                 or prop.startswith('http://sbols.org/genBankConversion#')):
    #             origin_props.append(prop)

    #     for origin_prop in origin_props:
    #         del sbol_obj.properties[origin_prop]

    #     sbol_obj.wasGeneratedBy = []

    @classmethod
    def copy_module_definition(cls, mod_definition, source_doc, sink_doc, import_namespace=False, deep_copy=False,
                               strip_prefixes=[]):
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
                except NotUniqueError as exc:
                    if exc.error_code() == sbol2.SBOLErrorCode.SBOL_ERROR_URI_NOT_UNIQUE:
                        mod_definition_copy = None

                        copy_version = copy_version + 1
                    else:
                        raise

            mod_definition_copy.wasDerivedFrom = [mod_definition.identity]
        else:
            try:
                return sink_doc.moduleDefinitions.get(mod_definition.identity)
            except RuntimeError:
                mod_definition_copy = mod_definition.copy(sink_doc)
            except NotFoundError as exc:
                if is_sbol_not_found(exc):
                    mod_definition_copy = mod_definition.copy(sink_doc)
                else:
                    raise

        FeatureLibrary.strip_origin_properties(mod_definition_copy, strip_prefixes)
        for sub_mod_copy in mod_definition_copy.modules:
            FeatureLibrary.strip_origin_properties(sub_mod_copy, strip_prefixes)
        for func_comp_copy in mod_definition_copy.functionalComponents:
            FeatureLibrary.strip_origin_properties(func_comp_copy, strip_prefixes)
        for intxn_copy in mod_definition_copy.interactions:
            FeatureLibrary.strip_origin_properties(intxn_copy, strip_prefixes)

            for parti_copy in intxn_copy.participations:
                FeatureLibrary.strip_origin_properties(parti_copy, strip_prefixes)

        if deep_copy:
            for sub_mod in mod_definition.modules:
                sub_mod_definition = source_doc.moduleDefinitions.get(sub_mod.definition)

                sub_mod_definition_copy = cls.copy_module_definition(sub_mod_definition, source_doc, sink_doc,
                                                                     import_namespace, deep_copy, strip_prefixes)

                sub_mod_copy = mod_definition_copy.modules.get(sub_mod.displayId)

                if sub_mod_definition_copy:
                    sub_mod_copy.definition = sub_mod_definition_copy.identity

            for func_comp in mod_definition.functionalComponents:
                sub_comp_definition = source_doc.componentDefinitions.get(func_comp.definition)

                sub_comp_definition_copy = FeatureLibrary.copy_component_definition(sub_comp_definition, source_doc, sink_doc)

                sub_comp_copy = mod_definition_copy.functionalComponents.get(func_comp.displayId)

                if sub_comp_definition_copy:
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
        self.logger = logging.getLogger('synbict')

    def add_sensors(self, target_doc, circuit_definition, sub_circuits, input_identities=set(),
            output_identities=set(), sensor_index=0, species_index=0):
        product_identities = set()

        for sub_circuit in sub_circuits:
            sub_circuit_definition = target_doc.moduleDefinitions.get(sub_circuit.identity)

            for intxn in sub_circuit_definition.interactions:
                if sbol2.SBO_GENETIC_PRODUCTION in intxn.types:
                    for par in intxn.participations:
                        if sbol2.SBO_PRODUCT in par.roles:
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

            abstraction_types = covered_sensors[i].get_abstraction_types()

            if len(abstraction_types) == 1:
                sensor_fc = self.get_circuit_species(covered_sensors[i].sensor_element, circuit_definition)

                if not sensor_fc:
                    sensor_fc = self.create_circuit_species(covered_sensors[i].sensor_element,
                                                sensor_doc,
                                                circuit_definition,
                                                species_index + k + 1,
                                                input_identities,
                                                output_identities)
                    
                    k = k + 1

                sensed_fc = self.get_circuit_species(covered_sensors[i].sensed_element, circuit_definition)

                if not sensed_fc:
                    sensed_fc = self.create_circuit_species(covered_sensors[i].sensed_element,
                                                sensor_doc,
                                                circuit_definition,
                                                species_index + k + 1,
                                                input_identities,
                                                output_identities)

                    k = k + 1
                
                if abstraction_types[0] == sbol2.SBO_STIMULATION:
                    sensor_intxn_ID = '_stimulates_'.join([sensed_fc.displayId, sensor_fc.displayId])
                elif abstraction_types[0] == sbol2.SBO_INHIBITION:
                    sensor_intxn_ID = '_inhibits_'.join([sensed_fc.displayId, sensor_fc.displayId])
                else:
                    sensor_intxn_ID = '_senses_'.join([sensor_fc.displayId, sensed_fc.displayId])

                sensor_intxn = circuit_definition.interactions.create(sensor_intxn_ID)
                sensor_intxn.types = abstraction_types

                if sensed_fc.name and sensor_fc.name:
                    if abstraction_types[0] == sbol2.SBO_STIMULATION:
                        sensor_intxn.name = ' stimulates '.join([sensed_fc.name, sensor_fc.name])
                    elif abstraction_types[0] == sbol2.SBO_INHIBITION:
                        sensor_intxn.name = ' inhibits '.join([sensed_fc.name, sensor_fc.name])
                    else:
                        sensor_intxn.name = ' senses '.join([sensor_fc.name, sensed_fc.name])

                sensor_par = sensor_intxn.participations.create(sensor_fc.displayId)
                sensed_par = sensor_intxn.participations.create(sensed_fc.displayId)

                if sensor_fc.name:
                    sensor_par.name = sensor_fc.name
                if sensed_fc.name:
                    sensed_par.name = sensed_fc.name

                sensor_par.participant = sensor_fc.identity
                sensed_par.participant = sensed_fc.identity

                if abstraction_types[0] == sbol2.SBO_STIMULATION:
                    sensor_par.roles = [sbol2.SBO_STIMULATED]
                    sensed_par.roles = [sbol2.SBO_STIMULATOR]
                elif abstraction_types[0] == sbol2.SBO_INHIBITION:
                    sensor_par.roles = [sbol2.SBO_INHIBITED]
                    sensed_par.roles = [sbol2.SBO_INHIBITOR]

                sensor_intxn.wasDerivedFrom = [covered_sensors[i].identity]
            else:
                sensor_definition = self.circuit_library.get_definition(covered_sensors[i].identity)

                CircuitLibrary.copy_module_definition(sensor_definition, sensor_doc, target_doc)

        return species_index + k

    @classmethod
    def get_circuit_species(cls, species_identity, circuit_definition):
        for species_fc in circuit_definition.functionalComponents:
            if species_fc.definition == species_identity:
                return species_fc

        return None

    @classmethod
    def create_circuit_species(cls, species_identity, species_doc, circuit_definition, species_index,
            input_identities=set(), output_identities=set()):
        species_ID = 'circuit_species_' + str(species_index)
        species_fc = circuit_definition.functionalComponents.create(species_ID)

        try:
            species_definition = species_doc.componentDefinitions.get(species_identity)

            if species_definition.name:
                species_fc.name = species_definition.name
        except RuntimeError:
            pass
        except NotFoundError as exc:
            if is_sbol_not_found(exc):
                pass
            else:
                raise

        species_fc.definition = species_identity

        if species_identity in input_identities:
            if species_identity in output_identities:
                species_fc.direction = sbol2.SBOL_DIRECTION_IN_OUT
            else:
                species_fc.direction = sbol2.SBOL_DIRECTION_IN
        elif species_identity in output_identities:
            species_fc.direction = sbol2.SBOL_DIRECTION_OUT

        return species_fc

    @classmethod
    def infer_devices(cls, target_doc, circuit_definition, constructs, circuit_feature_identities,
                      flanking_length):
        for construct in constructs:
            construct_definition = target_doc.componentDefinitions.get(construct.identity)

            annos = []

            for seq_anno in construct_definition.sequenceAnnotations:
                if (len(seq_anno.locations) == 1
                        and seq_anno.locations[0].getTypeURI() == sbol2.SBOL_RANGE):
                    anno_range = seq_anno.locations.getRange()

                    if seq_anno.component:
                        sub_comp = construct_definition.components.get(seq_anno.component)

                        feature_anno = FeatureAnnotation(anno_range.start,
                                                         anno_range.end,
                                                         sub_comp.definition,
                                                         display_ID=seq_anno.displayId,
                                                         name=seq_anno.name,
                                                         comp_ID=sub_comp.displayId,
                                                         comp_name=sub_comp.name,
                                                         orientation=anno_range.orientation)
                    else:
                        feature_anno = FeatureAnnotation(anno_range.start,
                                                         anno_range.end,
                                                         display_ID=seq_anno.displayId,
                                                         name=seq_anno.name,
                                                         orientation=anno_range.orientation)

                    annos.append(feature_anno)

            annos.sort()

            anno_group = []
            anno_groups = []

            device_end = -1

            for j in range(0, len(annos)):
                if device_end >= 0 and device_end + flanking_length <= annos[j].end:
                    anno_groups.append(anno_group)
                    anno_group = []
                    device_end = -1

                if annos[j].definition is not None and annos[j].definition in circuit_feature_identities:
                    if device_end < 0:
                        anno_group = [flanking_anno for flanking_anno in anno_group
                                      if flanking_anno.start + flanking_length >= annos[j].start]

                    device_end = annos[j].end

                anno_group.append(annos[j])

                if device_end >= 0 and j == len(annos) - 1:
                    anno_groups.append(anno_group)

            for i in range(0, len(anno_groups)):
                device_definition = sbol2.ComponentDefinition('_device_'.join([circuit_definition.displayId,
                                                                              str(i + 1)]),
                                                             sbol2.BIOPAX_DNA, '1')

                if circuit_definition.name:
                    device_definition.name = ' Device '.join([circuit_definition.name, str(i + 1)])

                inline_count = 0

                for anno in anno_groups[i]:
                    if anno.orientation == sbol2.SBOL_ORIENTATION_INLINE:
                        inline_count = inline_count + 1

                if inline_count/len(anno_groups[i]) < 0.5:
                    anno_groups[i].reverse()

                    rc_flag = True
                else:
                    rc_flag = False

                    for anno in anno_groups[i]:
                        if anno.orientation == sbol2.SBOL_ORIENTATION_INLINE:
                            anno.orientation = sbol2.SBOL_ORIENTATION_REVERSE_COMPLEMENT
                        else:
                            anno.orientation = sbol2.SBOL_ORIENTATION_INLINE

                if rc_flag:
                    device_start = anno_groups[i][0].end
                else:
                    device_start = anno_groups[i][0].start

                for anno in anno_groups[i]:
                    device_anno = device_definition.sequenceAnnotations.create(anno.display_ID)

                    device_anno.name = anno.name

                    anno_loc = device_anno.locations.createRange(device_anno.displayId + '_loc')
                    if rc_flag:
                        if anno.orientation == sbol2.SBOL_ORIENTATION_INLINE:
                            anno_loc.orientation = sbol2.SBOL_ORIENTATION_REVERSE_COMPLEMENT
                        else:
                            anno_loc.orientation = sbol2.SBOL_ORIENTATION_INLINE
                        anno_loc.start = -(anno.end - device_start - 1)
                        anno_loc.end = -(anno.start - device_start - 1)
                    else:
                        anno_loc.start = anno.start - device_start + 1
                        anno_loc.end = anno.end - device_start + 1
                        anno_loc.orientation = anno.orientation

                    if anno.comp_ID:
                        device_comp = device_definition.components.create(anno.comp_ID)

                        device_comp.name = anno.comp_name
                        device_comp.definition = anno.definition

                        device_anno.component = device_comp.identity

                target_doc.addComponentDefinition(device_definition)

                device_fc = circuit_definition.functionalComponents.create(device_definition.displayId)

                device_fc.name = device_definition.name
                device_fc.definition = device_definition.identity

    def infer_transcription(self, target_doc, circuit_definition, constructs, tx_threshold, species_index=0,
            input_identities=set(), output_identities=set()):
        k = species_index

        for construct in constructs:
            construct_definition = target_doc.componentDefinitions.get(construct.identity)

            inline_annos = []
            rc_annos = []

            for seq_anno in construct_definition.sequenceAnnotations:
                if (len(seq_anno.locations) == 1
                        and seq_anno.locations[0].getTypeURI() == sbol2.SBOL_RANGE):
                    anno_range = seq_anno.locations.getRange()

                    if seq_anno.component:
                        sub_comp = construct_definition.components.get(seq_anno.component)
                        try:
                            part_definition = target_doc.componentDefinitions.get(sub_comp.definition)

                            feature_anno = FeatureAnnotation(anno_range.start,
                                                             anno_range.end,
                                                             part_definition.identity,
                                                             part_definition.roles,
                                                             comp_identity=sub_comp.identity)
                        except NotFoundError as exc:
                            if is_sbol_not_found(exc):
                                feature_anno = FeatureAnnotation(anno_range.start,
                                                                 anno_range.end,
                                                                 roles=seq_anno.roles,
                                                                 identity=seq_anno.identity) 
                            else:
                                raise
                    else:
                        feature_anno = FeatureAnnotation(anno_range.start,
                                                         anno_range.end,
                                                         roles=seq_anno.roles,
                                                         identity=seq_anno.identity)

                    if anno_range.orientation == sbol2.SBOL_ORIENTATION_INLINE:
                        inline_annos.append(feature_anno)
                    elif anno_range.orientation == sbol2.SBOL_ORIENTATION_REVERSE_COMPLEMENT:
                        rc_annos.append(feature_anno)

            inline_annos.sort()
            rc_annos.sort()
            
            for i in range(0, len(inline_annos) - 1):
                for j in range(i, len(inline_annos)):
                    tx_distance = inline_annos[j].start - inline_annos[i].end

                    if tx_distance > tx_threshold:
                        j = len(inline_annos)
                    elif (tx_distance > 0 
                            and sbol2.SO_PROMOTER in inline_annos[i].roles
                            and sbol2.SO_CDS in inline_annos[j].roles):
                        
                        if not inline_annos[i].definition or not inline_annos[j].definition:
                            self.logger.debug('Inferred promoter-CDS stimulation between %s and %s in %s',
                                              (inline_annos[i].comp_identity + ' (' + inline_annos[i].definition  + ')'
                                                if inline_annos[i].definition else inline_annos[i].identity),
                                              (inline_annos[j].comp_identity + ' (' + inline_annos[j].definition  + ')'
                                                if inline_annos[j].definition else inline_annos[j].identity),
                                              construct_definition.identity)
                        else:
                            # fc1 = self.get_circuit_species(inline_annos[i].definition, circuit_definition)

                            # if not fc1:
                            fc1 = self.create_circuit_species(inline_annos[i].definition,
                                                              target_doc,
                                                              circuit_definition,
                                                              species_index + k + 1,
                                                              input_identities,
                                                              output_identities)

                            k = k + 1

                            # fc2 = self.get_circuit_species(inline_annos[j].definition, circuit_definition)

                            # if not fc2:
                            fc2 = self.create_circuit_species(inline_annos[j].definition,
                                                              target_doc,
                                                              circuit_definition,
                                                              species_index + k + 1,
                                                              input_identities,
                                                              output_identities)

                            k = k + 1

                            stimulation = circuit_definition.interactions.create('_stimulates_'.join([fc1.displayId, fc2.displayId]))
                            stimulation.types = [sbol2.SBO_STIMULATION]

                            if fc1.name and fc2.name:
                                stimulation.name = ' stimulates '.join([fc1.name, fc2.name])

                            stimulator = stimulation.participations.create(fc1.displayId)
                            stimulator.roles = [sbol2.SBO_STIMULATOR]
                            stimulator.participant = fc1.identity

                            stimulated = stimulation.participations.create(fc2.displayId)
                            stimulated.roles = [sbol2.SBO_STIMULATED]
                            stimulated.participant = fc2.identity

                            self.logger.debug('Inferred promoter-CDS stimulation %s between %s and %s in %s',
                                              stimulation.identity,
                                              inline_annos[i].comp_identity + ' (' + inline_annos[i].definition  + ')',
                                              inline_annos[j].comp_identity + ' (' + inline_annos[j].definition  + ')',
                                              construct_definition.identity)

            for i in range(0, len(rc_annos) - 1):
                for j in range(i, len(rc_annos)):
                    tx_distance = rc_annos[j].start - rc_annos[i].end

                    if tx_distance > tx_threshold:
                        j = len(rc_annos)
                    elif (tx_distance > 0
                            and sbol2.SO_CDS in rc_annos[i].roles
                            and sbol2.SO_PROMOTER in rc_annos[j].roles):

                        if not rc_annos[i].definition or not rc_annos[j].definition:
                            self.logger.debug('Inferred promoter-CDS stimulation between %s and %s in %s',
                                              (rc_annos[j].comp_identity + ' (' + rc_annos[j].definition  + ')'
                                                if rc_annos[j].definition else rc_annos[j].identity),
                                              (rc_annos[i].comp_identity + ' (' + rc_annos[i].definition  + ')'
                                                if rc_annos[i].definition else rc_annos[i].identity),
                                              construct_definition.identity)
                        else:
                            # fc1 = self.get_circuit_species(rc_annos[i].definition, circuit_definition)

                            # if not fc1:
                            fc1 = self.create_circuit_species(rc_annos[i].definition,
                                                              target_doc,
                                                              circuit_definition,
                                                              species_index + k + 1,
                                                              input_identities,
                                                              output_identities)

                            k = k + 1

                            # fc2 = self.get_circuit_species(rc_annos[j].definition, circuit_definition)

                            # if not fc2:
                            fc2 = self.create_circuit_species(rc_annos[j].definition,
                                                              target_doc,
                                                              circuit_definition,
                                                              species_index + k + 1,
                                                              input_identities,
                                                              output_identities)

                            k = k + 1

                            stimulation = circuit_definition.interactions.create('_stimulates_'.join([fc2.displayId, fc1.displayId]))
                            stimulation.types = [sbol2.SBO_STIMULATION]

                            if fc1.name and fc2.name:
                                stimulation.name = ' stimulates '.join([fc2.name, fc1.name])

                            stimulator = stimulation.participations.create(fc2.displayId)
                            stimulator.roles = [sbol2.SBO_STIMULATOR]
                            stimulator.participant = fc2.identity

                            stimulated = stimulation.participations.create(fc1.displayId)
                            stimulated.roles = [sbol2.SBO_STIMULATED]
                            stimulated.participant = fc1.identity
                            
                            self.logger.debug('Inferred promoter-CDS stimulation %s between %s and %s in %s',
                                              stimulation.identity,
                                              rc_annos[j].comp_identity + ' (' + rc_annos[j].definition  + ')',
                                              rc_annos[i].comp_identity + ' (' + rc_annos[i].definition  + ')',
                                              construct_definition.identity)

        return species_index + k

    @classmethod
    def map_circuit_species(cls, target_doc, circuit_definition, species_index=0):
        species_map = {}

        for species_fc in circuit_definition.functionalComponents:
            species_map[species_fc.definition] = species_fc.identity

        sub_species_map = {}

        for sub_circuit_mod in circuit_definition.modules:
            sub_circuit_definition = target_doc.moduleDefinitions.get(sub_circuit_mod.definition)

            for sub_species_fc in sub_circuit_definition.functionalComponents:
                if sub_species_fc.definition not in sub_species_map:
                    sub_species_map[sub_species_fc.definition] = []

                sub_species_map[sub_species_fc.definition].append((sub_circuit_mod.identity, sub_species_fc.identity))

        k = 0
        m = 0

        for species_identity in sub_species_map:
            if species_identity in species_map:
                local_identity = species_map[species_identity]
            elif len(sub_species_map[species_identity]) > 1:
                species_fc = circuit_definition.functionalComponents.create('_'.join(['circuit_species',
                                                                                      str(species_index + k + 1)]))

                k = k + 1

                species_fc.definition = species_identity

                local_identity = species_fc.identity
            else:
                local_identity = None

            if local_identity:
                for sub_species_pair in sub_species_map[species_identity]:
                    sub_circuit_mod = circuit_definition.modules.get(sub_species_pair[0])

                    maps_to = sub_circuit_mod.mapsTos.create('_'.join(['map', str(m + 1)]))

                    m = m + 1

                    maps_to.local = local_identity
                    maps_to.remote = sub_species_pair[1]
                    maps_to.refinement = sbol2.SBOL_REFINEMENT_VERIFY_IDENTICAL

        return species_index + k

    def build(self, circuit_id, target_doc, constructs, version='1', tx_threshold=0,
            input_identities=set(), output_identities=set(), infer_sensors=False,
            infer_devices=False, flanking_length=0):
        circuit_definition = sbol2.ModuleDefinition(circuit_id, version)
        circuit_definition.roles = [self.NCIT_BIOCHEMICAL_PATHWAY]

        self.logger.info('Building %s', circuit_definition.identity)

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
                        self.logger.warning('No sub-circuits found for construct %s', construct.identity)

                # if infer_devices:
                #     self.infer_devices(target_doc, circuit_definition, covered_constructs,
                #                        circuit_feature_identities, flanking_length)

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
                        if sub_fc.definition in input_identities or sub_fc.definition in output_identities:
                            species_fc = self.get_circuit_species(sub_fc.definition, circuit_definition)

                            if not species_fc:
                                self.create_circuit_species(sub_fc.definition,
                                                            sub_circuit_doc,
                                                            circuit_definition,
                                                            k + 1,
                                                            input_identities,
                                                            output_identities)

                                k = k + 1

                    CircuitLibrary.copy_module_definition(sub_circuit_definition, sub_circuit_doc, target_doc, deep_copy=True)

                    self.logger.debug('Added sub-circuit %s', sub_circuit_definition.identity)

                if infer_sensors:
                    k = self.add_sensors(target_doc, circuit_definition, covered_circuits, input_identities,
                        output_identities, len(covered_circuits), k)

                k = self.infer_transcription(target_doc, circuit_definition, constructs, tx_threshold, k,
                    input_identities, output_identities)

                self.map_circuit_species(target_doc, circuit_definition, k)

                target_doc.addModuleDefinition(circuit_definition)

                self.logger.info('Finished building %s', circuit_definition.identity)

                return True
            else:
                k = 0

                k = self.infer_transcription(target_doc, circuit_definition, constructs, tx_threshold, k,
                    input_identities, output_identities)

                if k > 0:
                    target_doc.addModuleDefinition(circuit_definition)

                    self.logger.info('Finished building %s', circuit_definition.identity)
                else:
                    self.logger.error('Failed to build %s, no sub-circuits or transcriptional relationships found for constructs', circuit_definition.identity)

                    return False
        else:
            self.logger.error('Failed to build %s, no constructs found with minimum length', circuit_definition.identity)

            return False

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()

    # Common arguments
    parser.add_argument('-n', '--namespace')
    parser.add_argument('-c', '--sub_circuit_files', nargs='*', default=[])
    parser.add_argument('-nb', '--no_build', action='store_true')
    parser.add_argument('-l', '--log_file', nargs='?', default='')
    parser.add_argument('-v', '--validate', action='store_true')

    # Circuit inference arguments
    parser.add_argument('-t', '--target_files', nargs='*', default=[])
    parser.add_argument('-i', '--circuit_IDs', nargs='*', default=[])
    parser.add_argument('-s', '--circuit_suffix', nargs='?', default='')
    parser.add_argument('-cv', '--circuit_version', nargs='?', default='1')
    parser.add_argument('-o', '--output_files', nargs='*', default=[])
    parser.add_argument('-os', '--output_suffix', nargs='?', default='')
    parser.add_argument('-ii', '--input_identities', nargs='*', default=[])
    parser.add_argument('-oi', '--output_identities', nargs='*', default=[])
    parser.add_argument('-m', '--min_target_length', nargs='?', default='2000')
    parser.add_argument('-is', '--infer_sensors', action='store_true')
    parser.add_argument('-d', '--tx_threshold', nargs='?', default='200')
    parser.add_argument('-f', '--flanking_length', nargs='?', default='200')
    parser.add_argument('-iv', '--infer_devices', action='store_true')
    parser.add_argument('-sp', '--strip_prefixes', nargs='*', default=[])

    # Sub-circuit library extension arguments
    parser.add_argument('-e', '--extend_sub_circuits', action='store_true')
    parser.add_argument('-xs', '--extension_suffix', nargs='?', default='')
    parser.add_argument('-x', '--extension_threshold', nargs='?', default='0.05')
    
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

    circuit_docs = []
    for circuit_file in args.sub_circuit_files:
        circuit_docs.append(load_sbol(circuit_file))

    if args.extend_sub_circuits:
        circuit_library = CircuitLibrary(circuit_docs, True)

        circuit_library.extend_circuits_by_name(float(args.extension_threshold), args.strip_prefixes)

        for extended_doc in circuit_library.get_updated_documents():
            (extended_file_base, extended_file_extension) = os.path.splitext(extended_doc.name)

            if len(args.extension_suffix) > 0:
                extended_file = '_'.join([extended_file_base, args.extension_suffix]) + '.xml'
            else:
                extended_file = extended_file_base + '_extended.xml'

            logger.info('Writing %s', extended_file)

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
                                     os.path.isfile(os.path.join(target_file, tf)) and tf.endswith('.xml')])
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

            unique_ID = circuit_ID

            while unique_ID in circuit_memo:
                circuit_index = circuit_index + 1

                unique_ID = circuit_ID + str(circuit_index)
            
            circuit_memo.add(unique_ID)

            build_success = circuit_builder.build(circuit_ID,
                                                  target_doc,
                                                  target_library.get_features(int(args.min_target_length), True),
                                                  args.circuit_version,
                                                  int(args.tx_threshold),
                                                  set(args.input_identities),
                                                  set(args.output_identities),
                                                  args.infer_sensors,
                                                  args.infer_devices,
                                                  int(args.flanking_length))

            if build_success:
                if len(args.output_files) == 1 and os.path.isdir(args.output_files[0]):
                    (target_file_path, target_filename) = os.path.split(target_files[i])
                    (target_file_base, target_file_extension) = os.path.splitext(target_filename)

                    if len(args.output_suffix) > 0:
                        output_file = os.path.join(args.output_files[0], '_'.join([target_file_base, args.output_suffix + target_file_extension]))
                    else:
                        output_file = os.path.join(args.output_files[0], target_file_base + target_file_extension)

                elif i < len(args.output_files):
                    output_file = args.output_files[i]
                else:
                    (target_file_base, target_file_extension) = os.path.splitext(target_files[i])

                    if len(args.output_suffix) > 0:
                        output_file = '_'.join([target_file_base, args.output_suffix + target_file_extension])
                    else:
                        output_file = target_file_base + target_file_extension

                if sbol2.Config.getOption('validate') == True:
                    logger.info('Validating and writing %s', output_file)
                else:
                    logger.info('Writing %s', output_file)

                target_doc.write(output_file)

                logger.info('Finished writing %s', output_file)

    logger.info('Finished curating')


class CircuitIdentityError(Exception):

    def __init__(self, circuit_identity, feature_ID):
        self.circuit_identity = circuit_identity
        self.feature_ID = feature_ID

    def __str__(self):
        return "Circuit identity {ci} does not contain the feature {fi}.".format(ci=self.circuit_identity, fi=self.feature_ID)


if __name__ == '__main__':
    main()
