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
    print('Loading ' + sbol_file + '...')

    doc = Document()
    doc.read(sbol_file)

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')

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

        print('Loading circuits...\n')

        for i in range(0, len(self.docs)):
            for mod_definition in self.docs[i].moduleDefinitions:
                features = []

                is_complete = True

                features = self.__extract_features(self.docs[i], mod_definition)

                if len(features) > 0 and is_complete:
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
                else:
                    logging.warning('%s is incomplete and was not loaded.', mod_definition.identity)

        print(self.__template_to_product)
        print(self.__activator_to_dna)
        print(self.__repressor_to_dna)
        print(self.__dna_to_dna_repression)
        print(self.__dna_to_dna_activation)

    @classmethod
    def __extract_features(cls, doc, mod_definition):
        features = []

        for func_comp in mod_definition.functionalComponents:
            try:
                comp_definition = doc.getComponentDefinition(func_comp.definition)
            except:
                comp_definition = None

            if comp_definition is None:
                logging.warning('%s not found.', func_comp.definition)

                is_complete = False
            elif BIOPAX_DNA in comp_definition.types:
                dna_seqs = []

                for seq_URI in comp_definition.sequences:
                    try:
                        seq = doc.getSequence(seq_URI)
                    except RuntimeError:
                        seq = None

                    if seq is not None and seq.encoding == SBOL_ENCODING_IUPAC:
                        dna_seqs.append(seq)

                if len(dna_seqs) > 0:
                    features.append(Feature(dna_seqs[0].elements, comp_definition.identity, set(comp_definition.roles),
                        comp_definition.wasDerivedFrom))
                else:
                    logging.warning('DNA sequence for %s not found.', comp_definition.identity)

        for sub_mod in mod_definition.modules:
            try:
                sub_definition = doc.getModuleDefinition(sub_mod.definition)
            except:
                sub_definition = None

            if sub_definition is not None:
                features.extend(cls.__extract_features(doc, sub_definition))

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

class CircuitBuilder():

    def __init__(self, circuit_library):
        self.circuit_library = circuit_library

    @classmethod
    def __move_component_definition(cls, source_doc, sink_doc, component_definition):
        try:
            sink_doc.addComponentDefinition(component_definition)
        except:
            pass

        for sub_comp in component_definition.components:
            try:
                sub_definition = source_doc.getComponentDefinition(sub_comp.definition)

                cls.__move_component_definition(source_doc, sink_doc, sub_definition)
            except RuntimeError:
                pass

    @classmethod
    def __move_module_definition(cls, source_doc, sink_doc, module_definition):
        try:
            sink_doc.addModuleDefinition(module_definition)
        except:
            pass

        for sub_mod in module_definition.modules:
            try:
                sub_definition = source_doc.getModuleDefinition(sub_mod.definition)

                cls.__move_module_definition(source_doc, sink_doc, sub_definition)
            except RuntimeError:
                pass

        for func_comp in module_definition.functionalComponents:
            try:
                sub_definition = source_doc.getComponentDefinition(func_comp.definition)

                cls.__move_component_definition(source_doc, sink_doc, sub_definition)
            except RuntimeError:
                pass

    def build(self, circuit_id, target_doc, target_library, min_feature_length):
        covered_circuits = []

        for circuit in self.circuit_library.circuits:
            if circuit.is_covered(target_library):
                covered_circuits.append(circuit)

        if len(covered_circuits) > 0:
            circuit_definition = ModuleDefinition(circuit_id, '1')

            print('Building ' + circuit_definition.identity + '...')

            features = target_library.get_features(min_feature_length)

            for i in range(0, len(features)):
                func_comp = circuit_definition.functionalComponents.create('construct_' + str(i))

                func_comp.definition = features[i].identity

            for i in range(0, len(covered_circuits)):
                sub_mod = circuit_definition.modules.create('sub_circuit_' + str(i))

                sub_mod.definition = covered_circuits[i].identity

                sub_circuit_doc = self.circuit_library.get_document(covered_circuits[i].identity)

                sub_circuit_definition = self.circuit_library.get_definition(covered_circuits[i].identity)

                self.__move_module_definition(sub_circuit_doc, target_doc, sub_circuit_definition)

            target_doc.addModuleDefinition(circuit_definition)

            logging.info('Finished building %s.\n', circuit_definition.identity)
        else:
            logging.info('Failed to build %s.\n', circuit_definition.identity)

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--namespace')
    parser.add_argument('-i', '--circuit_id')
    parser.add_argument('-t', '--target_files', nargs='+')
    parser.add_argument('-c', '--circuit_files', nargs='*', default=[])
    parser.add_argument('-m', '--min_feature_length', nargs='?', default=1000)
    parser.add_argument('-l', '--curation_log', nargs='?', default='')
    parser.add_argument('-v', '--validate', action='store_true')
    
    args = parser.parse_args(args)

    if len(args.curation_log) > 0:
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)

        logging.basicConfig(level=logging.INFO, filename=args.curation_log, filemode='w',
                            format='%(levelname)s : %(message)s')

    setHomespace(args.namespace)
    Config.setOption('validate', args.validate)
    Config.setOption('sbol_typed_uris', False)

    circuit_docs = []
    for circuit_file in args.circuit_files:
        circuit_docs.append(load_sbol(circuit_file))
    circuit_library = CircuitLibrary(circuit_docs)

    circuit_builder = CircuitBuilder(circuit_library)

    for i in range (0, len(args.target_files)):
        target_doc = load_sbol(args.target_files[i])
        
        target_library = FeatureLibrary([target_doc])

        circuit_builder.build(args.circuit_id, target_doc, target_library, args.min_feature_length)

        (target_file_base, file_extension) = os.path.splitext(args.target_files[i])
        target_doc.write('_'.join([target_file_base, 'circuit']) + file_extension)

    print('Finished curating.')

if __name__ == '__main__':
    main()
