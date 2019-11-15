import logging
import argparse
import os
import sys

from Bio.Seq import Seq
from sbol import *
from flashtext import KeywordProcessor
from sequences_to_features import Feature
from sequences_to_features import FeatureLibrary

def load_sbol(sbol_file):
    logging.info('Loading %s', sbol_file)

    doc = Document()
    doc.read(sbol_file)

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')

    logging.info('Finished loading %s', sbol_file)

    return doc

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

class CircuitLibrary():

    SBO_TEMPLATE = 'http://identifiers.org/biomodels.sbo/SBO:0000645'

    def __init__(self, docs):
        self.circuits = []
        self.docs = docs
        self.__circuit_map = {}
        self.__template_to_product = {}
        self.__activator_to_dna= {}
        self.__repressor_to_dna = {}
        self.__dna_to_dna_repression = {}
        self.__dna_to_dna_activation = {}

        logging.info('Loading circuits')

        for i in range(0, len(self.docs)):
            for mod_definition in self.docs[i].moduleDefinitions:
                features = self.__extract_features(self.docs[i], mod_definition)

                if len(features) > 0:
                    for intxn in mod_definition.interactions:
                        if SBO_STIMULATION in intxn.types:
                            stimulator = None
                            stimulated = None

                            for par in intxn.participations:
                                if SBO_STIMULATOR in par.roles:
                                    stimulator = mod_definition.functionalComponents.get(par.participant).definition
                                elif SBO_STIMULATED in par.roles:
                                    stimulated = mod_definition.functionalComponents.get(par.participant).definition

                            if stimulator is not None and stimulated is not None:
                                if stimulator not in self.__activator_to_dna:
                                    self.__activator_to_dna[stimulator] = []

                                self.__activator_to_dna[stimulator].append(stimulated)
                        elif SBO_INHIBITION in intxn.types:
                            inhibitor = None
                            inhibited = None

                            for par in intxn.participations:
                                if SBO_INHIBITOR in par.roles:
                                    inhibitor = mod_definition.functionalComponents.get(par.participant).definition
                                elif SBO_INHIBITED in par.roles:
                                    inhibited = mod_definition.functionalComponents.get(par.participant).definition

                            if inhibitor is not None and inhibited is not None:
                                if inhibitor not in self.__repressor_to_dna:
                                    self.__repressor_to_dna[inhibitor] = []

                                self.__repressor_to_dna[inhibitor].append(inhibited)
                        elif SBO_GENETIC_PRODUCTION in intxn.types:
                            template = None
                            product = None

                            for par in intxn.participations:
                                if SBO_PRODUCT in par.roles:
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

        logging.info('Finished loading circuits')

    @classmethod
    def __extract_features(cls, doc, mod_definition):
        features = []

        for func_comp in mod_definition.functionalComponents:
            try:
                comp_definition = doc.getComponentDefinition(func_comp.definition)
            except:
                comp_definition = None

            if comp_definition is None:
                logging.warning('%s was not loaded since its sub-ComponentDefinition %s was not found', mod_definition.identity, func_comp.definition)

                return []
            elif BIOPAX_DNA in comp_definition.types:
                features.append(Feature('', comp_definition.identity, set(comp_definition.roles),
                    parent_identities=comp_definition.wasDerivedFrom))

        for sub_mod in mod_definition.modules:
            try:
                sub_definition = doc.getModuleDefinition(sub_mod.definition)
            except:
                sub_definition = None

            if sub_definition is not None:
                sub_features = cls.__extract_features(doc, sub_definition)

                if len(sub_features) == 0:
                    logging.warning('%s was not loaded since its sub-ModuleDefinition %s is incomplete or contains no DNA features',
                        mod_definition.identity, sub_mod.definition)

                    return []
                else:
                    features.extend(sub_features)

        if len(features) == 0:
            logging.warning('%s was not loaded since it contains no DNA features', mod_definition.identity)

        return features

    def get_document(self, identity):
        return self.docs[self.__circuit_map[identity]]

    def get_definition(self, identity):
        return self.get_document(identity).getModuleDefinition(identity)

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

    @classmethod
    def copy_module_definition(cls, mod_definition, source_doc, sink_doc, is_recursive=True):
        definition_copy = mod_definition.copy(sink_doc)

        if is_recursive:
            for sub_mod in mod_definition.modules:
                try:
                    sub_mod_definition = source_doc.getModuleDefinition(sub_mod.definition)

                    cls.copy_module_definition(sub_mod_definition, source_doc, sink_doc)
                except RuntimeError:
                    pass

            for func_comp in mod_definition.functionalComponents:
                try:
                    sub_comp_definition = source_doc.getComponentDefinition(func_comp.definition)

                    FeatureLibrary.copy_component_definition(sub_comp_definition, source_doc, sink_doc)
                except RuntimeError:
                    pass

class CircuitBuilder():

    def __init__(self, circuit_library):
        self.circuit_library = circuit_library

    def build(self, circuit_id, target_doc, target_library, min_target_length, version='1'):
        circuit_definition = ModuleDefinition(circuit_id, version)

        logging.info('Building %s', circuit_definition.identity)

        constructs = target_library.get_features(min_target_length)

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
            
                for i in range(0, len(covered_constructs)):
                    func_comp = circuit_definition.functionalComponents.create('construct_' + str(i + 1))

                    func_comp.definition = covered_constructs[i].identity

                for i in range(0, len(covered_circuits)):
                    sub_mod = circuit_definition.modules.create('sub_circuit_' + str(i + 1))

                    sub_mod.definition = covered_circuits[i].identity

                    sub_circuit_doc = self.circuit_library.get_document(covered_circuits[i].identity)

                    sub_circuit_definition = self.circuit_library.get_definition(covered_circuits[i].identity)

                    CircuitLibrary.copy_module_definition(sub_circuit_definition, sub_circuit_doc, target_doc)

                target_doc.addModuleDefinition(circuit_definition)

                logging.info('Finished building %s', circuit_definition.identity)
            else:
                logging.warning('Failed to build %s, no sub-circuits found for constructs', circuit_definition.identity)
        else:
            logging.warning('Failed to build %s, no constructs found with minimum length %s', circuit_definition.identity, str(min_target_length))

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--namespace')
    parser.add_argument('-t', '--target_files', nargs='+')
    parser.add_argument('-c', '--sub_circuit_files', nargs='+')
    parser.add_argument('-i', '--circuit_IDs', nargs='*', default=[])
    parser.add_argument('-o', '--output_files', nargs='*', default=[])
    parser.add_argument('-m', '--min_target_length', nargs='?', default=2000)
    parser.add_argument('-l', '--curation_log', nargs='?', default='')
    parser.add_argument('-v', '--version', nargs='?', default='1')
    parser.add_argument('-x', '--validate', action='store_true')
    
    args = parser.parse_args(args)

    
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    if len(args.curation_log) > 0:
        logging.basicConfig(level=logging.DEBUG, filename=args.curation_log, filemode='w',
                                format='%(levelname)s : %(message)s')
    else:
        logging.basicConfig(level=logging.DEBUG, format='%(levelname)s : %(message)s')

    setHomespace(args.namespace)
    Config.setOption('validate', args.validate)
    Config.setOption('sbol_typed_uris', False)

    circuit_docs = []
    for circuit_file in args.sub_circuit_files:
        circuit_docs.append(load_sbol(circuit_file))
    circuit_library = CircuitLibrary(circuit_docs)

    circuit_builder = CircuitBuilder(circuit_library)

    circuit_memo = set()
    circuit_index = 0

    for i in range (0, len(args.target_files)):
        target_doc = load_sbol(args.target_files[i])
        
        target_library = FeatureLibrary([target_doc], require_sequence=False)

        if i < len(args.circuit_IDs):
            circuit_ID = args.circuit_IDs[i]
        else:
            (circuit_ID, file_extension) = os.path.splitext(os.path.basename(args.target_files[i]))

            circuit_ID = circuit_ID + '_circuit'

        unique_ID = circuit_ID

        while unique_ID in circuit_memo:
            circuit_index = circuit_index + 1

            unique_ID = circuit_ID + str(circuit_index)
        
        circuit_memo.add(unique_ID)

        circuit_builder.build(circuit_ID, target_doc, target_library, int(args.min_target_length),
            args.version)

        if i < len(args.output_files):
            output_file = args.output_files[i]
        else:
            (target_file_base, file_extension) = os.path.splitext(args.target_files[i])
            output_file = target_file_base + '_circuit' + file_extension

        if Config.getOption('validate') == True:
            logging.info('Validating and writing %s', output_file)
        else:
            logging.info('Writing %s', output_file)

        target_doc.write(output_file)

        logging.info('Finished writing %s', output_file)

    logging.info('Finished curating')

if __name__ == '__main__':
    main()
