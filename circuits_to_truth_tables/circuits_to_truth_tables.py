import logging
import argparse
import os
import sys
import csv
from itertools import product

import sbol

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

class CircuitSpecies():

    def __init__(self, identity, canonical_identity, is_canonical=False):
        self.identity = identity
        self.canonical_identity = canonical_identity

    def get_identity(self, is_canonical=True):
        if is_canonical:
            return self.canonical_identity
        else:
            return self.identity

class LogicGate():

    AND = 'AND'
    NAND = 'NAND'
    OR = 'OR'
    NOR = 'NOR'
    NOT = 'NOT'
    YES = 'YES'

    def __init__(self, gate_output, gate_inputs, gate_type):
        if len(gate_inputs) == 0:
            pass
        elif len(gate_inputs) > 2:
            pass

        self.__gate_output = gate_output
        self.__gate_inputs = gate_inputs
        self.__gate_type = gate_type

    def calculate_output_value(self, input_map, is_canonical=True):
        for gate_input in self.__gate_inputs:
            if gate_input.get_identity(is_canonical) not in input_map:
                pass

        if len(self.__gate_inputs) == 1:
            input_1_value = input_map[self.__gate_inputs[0].get_identity(is_canonical)]

            if self.__gate_type == self.YES:
                return input_1_value
            elif self.__gate_type == self.NOT:
                return ~(input_1_value - 2)
            else:
                return -1
        elif len(self.__gate_inputs) == 2:
            input_1_value = input_map[self.__gate_inputs[0].get_identity(is_canonical)]
            input_2_value = input_map[self.__gate_inputs[1].get_identity(is_canonical)]

            if self.__gate_type == self.OR:
                return input_1_value | input_2_value
            elif self.__gate_type == self.NOR:
                return ~(input_1_value | input_2_value - 2)
            elif self.__gate_type == self.AND:
                return input_1_value & input_2_value
            elif self.__gate_type == self.NAND:
                return ~(input_1_value & input_2_value - 2)
            else:
                return -1
        else:
            return -1

    def get_inputs(self):
        return self.__gate_inputs

    def get_output(self):
        return self.__gate_output

    def serialize(self, is_canonical=True):
        if len(self.__gate_inputs) == 1:
            input_1_identity = self.__gate_inputs[0].get_identity(is_canonical)
            output_identity = self.__gate_output.get_identity(is_canonical)

            return ' '.join([self.__gate_type, input_1_identity, '->', output_identity])
        elif len(self.__gate_inputs) == 2:
            input_1_identity = self.__gate_inputs[0].get_identity(is_canonical)
            input_2_identity = self.__gate_inputs[1].get_identity(is_canonical)

            return ' '.join([input_1_identity, self.__gate_type, input_2_identity, '->', self.__gate_output])
        else:
            return ''

class LogicCircuit():

    SBO_TEMPLATE = 'http://identifiers.org/biomodels.sbo/SBO:0000645'
    SBO_CODING = 'http://identifiers.org/biomodels.sbo/SBO:0000335'
    SBO_REGULATORY = 'http://identifiers.org/biomodels.sbo/SBO:0000369'
    SBO_CONTROL = 'http://identifiers.org/biomodels.sbo/SBO:0000168'
    SBO_MODIFIER = 'http://identifiers.org/biomodels.sbo/SBO:0000019'
    SBO_MODIFIED = 'http://identifiers.org/biomodels.sbo/SBO:0000644'

    def __init__(self, circuit_def=None, circuit_doc=None, product_to_gates={}, is_canonical=True):
        self.__product_to_gates = product_to_gates

        self.__circuit_inputs = []
        self.__circuit_outputs = []
        self.__species_dict = {}

        if circuit_def and circuit_doc:
            production_map = self.__build_production_map(circuit_def, circuit_doc, is_canonical=is_canonical)
            repression_map = self.__build_repression_map(circuit_def, circuit_doc, is_canonical=is_canonical)
            activation_map = self.__build_activation_map(circuit_def, circuit_doc, is_canonical=is_canonical)
            transcription_map = self.__distill_transcription_map(circuit_def, circuit_doc, activation_map, is_canonical=is_canonical)
            
            # print('Production')
            # for key in production_map:
            #     for producer in production_map[key]:
            #         print(producer.get_identity(is_canonical) + ' : ' + key)
            # print('Repression')
            # for key in repression_map:
            #     for repressors in repression_map[key]:
            #         for repressor in repressors:
            #             print(repressor.get_identity(is_canonical) + ' : ' + key)
            # print('Activation')
            # for key in activation_map:
            #     for activators in activation_map[key]:
            #         for activator in activators:
            #             print(activator.get_identity(is_canonical) + ' : ' + key)
            # print('Transcription')
            # for key in transcription_map:
            #     for transcriber in transcription_map[key]:
            #         print(transcriber.get_identity(is_canonical) + ' : ' + key)

            for fc in circuit_def.functionalComponents:
                if fc.direction == sbol.SBOL_DIRECTION_IN:
                    if is_canonical:
                        self.__circuit_inputs.append(self.__species_dict[fc.definition])
                    else:
                        self.__circuit_inputs.append(self.__species_dict[fc.identity])
                elif fc.direction == sbol.SBOL_DIRECTION_OUT:
                    if is_canonical:
                        self.__circuit_outputs.append(self.__species_dict[fc.definition])
                    else:
                        self.__circuit_outputs.append(self.__species_dict[fc.identity])

            self.__product_to_gates.update(self.__build_gate_map(production_map,
                                                                 repression_map,
                                                                 activation_map,
                                                                 transcription_map,
                                                                 is_canonical))

        if len(self.__circuit_inputs) == 0 or len(self.__circuit_outputs) == 0:
            gate_input_dict = {
                                gate_input.get_identity(is_canonical) : gate_input
                                for product_identity in self.__product_to_gates
                                for gate in self.__product_to_gates[product_identity]
                                for gate_input in gate.get_inputs()
                              }

            gate_output_dict = {
                                product_identity : self.__species_dict[product_identity]
                                for product_identity in self.__product_to_gates
                               }

            if len(self.__circuit_inputs) == 0:
                self.__circuit_inputs.extend(
                                                [
                                                    gate_input_dict[gate_input_identity]
                                                    for gate_input_identity in gate_input_dict
                                                    if gate_input_identity not in gate_output_dict
                                                ]
                                            )

            if len(self.__circuit_outputs) == 0:
                self.__circuit_outputs.extend(
                                                [
                                                    gate_output_dict[gate_output_identity]
                                                    for gate_output_identity in gate_output_dict
                                                    if gate_output_identity not in gate_input_dict
                                                ]
                                            )

    def get_intermediates(self, is_canonical=True):
        circuit_input_identities = {
                                    circuit_input.get_identity(is_canonical)
                                    for circuit_input in self.__circuit_inputs
                                   }

        circuit_output_identities = {
                                        circuit_output.get_identity(is_canonical)
                                        for circuit_output in self.__circuit_outputs
                                    }

        circuit_intermediates = [
                                    self.__product_to_gates[product_identity]
                                    for product_identity in self.__product_to_gates
                                    if (
                                        product_identity not in circuit_input_identities
                                        and product_identity not in circuit_output_identities
                                    )
                                ]

        return circuit_intermediates

    def serialize_truth_table(self, truth_table, is_canonical=True):
        truth_csv = []

        circuit_intermediates = self.get_intermediates(is_canonical)

        truth_csv.append(
                            [
                                'input'
                                for i in range(0, len(self.__circuit_inputs))
                            ] +
                            [
                                'intermediate'
                                for i in range(0, len(circuit_intermediates))
                            ] +
                            [
                                'output'
                                for i in range(0, len(self.__circuit_outputs))
                            ]
                        )

        circuit_species = self.__circuit_inputs + circuit_intermediates + self.__circuit_outputs

        truth_csv.append(
                            [
                                cs.get_identity(is_canonical)
                                for cs in circuit_species
                            ]
                        )

        for truth_row in truth_table:
            truth_csv.append(
                [
                    truth_row[cs.get_identity(is_canonical)]
                    for cs in circuit_species
                ]
            )

        return truth_csv

    def serialize(self):
        return '\n'.join(
                            [
                                gate.serialize() 
                                for product_identity in self.__product_to_gates
                                for gate in self.__product_to_gates[product_identity]
                            ]
                        )

    def __build_gate_map(self, production_map, repression_map, activation_map, transcription_map, is_canonical=True):
        product_to_gates = {}

        for product_identity in production_map:
            for template in production_map[product_identity]:
                activation_flag = False
                repression_flag = False
                cooperativity_flag = False

                gate_input_dict = {}
                promoters = []

                template_identity = template.get_identity(is_canonical)

                if template_identity in transcription_map:
                    promoters.extend(transcription_map[template_identity])

                    for promoter in transcription_map[template_identity]:
                        promoter_identity = promoter.get_identity(is_canonical)

                        if promoter_identity in activation_map:
                            for activators in activation_map[promoter_identity]:
                                if len(activators) == 2:
                                    cooperativity_flag = True

                                for activator in activators:
                                    gate_input_dict[activator.get_identity(is_canonical)] = activator

                            activation_flag = True
                        elif promoter_identity in repression_map:
                            for repressors in repression_map[promoter_identity]:
                                if len(repressors) == 2:
                                    cooperativity_flag = True

                                for repressor in repressors:
                                    gate_input_dict[repressor.get_identity(is_canonical)] = repressor

                            repression_flag = True
                elif template_identity in activation_map:
                    for activators in activation_map[template_identity]:
                        if len(activators) == 2:
                            cooperativity_flag = True

                        for activator in activators:
                            gate_input_dict[activator.get_identity(is_canonical)] = activator

                    activation_flag = True
                elif template_identity in repression_map:
                    for repressors in repression_map[template_identity]:
                        if len(repressors) == 2:
                            cooperativity_flag = True

                        for repressor in repressors:
                            gate_input_dict[repressor.get_identity(is_canonical)] = repressor

                    repression_flag = True

                if not repression_flag and not activation_flag:
                    pass
                elif repression_flag and activation_flag:
                    pass
                elif len(promoters) > 2:
                    pass
                else:
                    if product_identity not in product_to_gates:
                        product_to_gates[product_identity] = []

                    product = self.__species_dict[product_identity]

                    gate_inputs = [gate_input_dict[gate_input_identity] for gate_input_identity in gate_input_dict]

                    if len(gate_inputs) == 1:
                        if activation_flag:
                            product_to_gates[product_identity].append(LogicGate(product, gate_inputs, LogicGate.YES))
                        else:
                            product_to_gates[product_identity].append(LogicGate(product, gate_inputs, LogicGate.NOT))
                    else:
                        if activation_flag:
                            if len(promoters) < 2:
                                if cooperativity_flag:
                                    product_to_gates[product_identity] = LogicGate(product, gate_inputs, LogicGate.AND)
                                else:
                                    product_to_gates[product_identity] = LogicGate(product, gate_inputs, LogicGate.OR)
                            else:
                                product_to_gates[product_identity] = LogicGate(product, gate_inputs, LogicGate.OR)
                        else:
                            if len(promoters) < 2:
                                if cooperativity_flag:
                                    product_to_gates[product_identity] = LogicGate(product, gate_inputs, LogicGate.NAND)
                                else:
                                    product_to_gates[product_identity] = LogicGate(product, gate_inputs, LogicGate.NOR)
                            else:
                                product_to_gates[product_identity] = LogicGate(product, gate_inputs, LogicGate.NAND)

        return product_to_gates

    def __build_production_map(self, circuit_def, circuit_doc, production_map={}, is_canonical=True):
        for intxn in circuit_def.interactions:
            if sbol.SBO_GENETIC_PRODUCTION in intxn.types:
                regulator_dict = {}

                product = None
                producer = None

                for parti in intxn.participations:
                    if sbol.SBO_PRODUCT in parti.roles:
                        product_fc = circuit_def.functionalComponents.get(parti.participant)
                        product = CircuitSpecies(product_fc.identity, product_fc.definition)
                    elif self.SBO_TEMPLATE in parti.roles or self.SBO_CODING in parti.roles:
                        producer_fc = circuit_def.functionalComponents.get(parti.participant)
                        producer = CircuitSpecies(producer_fc.identity, producer_fc.definition)
                    elif self.SBO_REGULATORY in parti.roles:
                        regulator_fc = circuit_def.functionalComponents.get(parti.participant)
                        regulator = CircuitSpecies(regulator_fc.identity, regulator_fc.definition)
                        regulator_dict[regulator.get_identity(is_canonical)] = regulator

                if product:
                    product_identity = product.get_identity(is_canonical)

                    self.__species_dict[product_identity] = product

                    if len(regulator_dict) > 0:
                        production_map[product_identity] = [
                                                            regulator_dict[regulator_identity]
                                                            for regulator_identity in regulator_dict
                                                           ]

                        self.__species_dict.update(regulator_dict)
                    elif producer:
                        production_map[product_identity] = [producer]

                        self.__species_dict[producer.get_identity(is_canonical)] = producer

        for sub_circuit in circuit_def.modules:
            sub_circuit_def = circuit_doc.moduleDefinitions.get(sub_circuit.definition)

            self.__build_production_map(sub_circuit_def, circuit_doc, production_map, is_canonical)

        return production_map

    def __build_activation_map(self, circuit_def, circuit_doc, activation_map={}, is_canonical=True):
        for intxn in circuit_def.interactions:
            if sbol.SBO_STIMULATION in intxn.types:
                activator_dict = {}
                activated = None

                for parti in intxn.participations:
                    if sbol.SBO_STIMULATOR in parti.roles:
                        activator_fc = circuit_def.functionalComponents.get(parti.participant)
                        activator = CircuitSpecies(activator_fc.identity, activator_fc.definition)
                        activator_dict[activator.get_identity(is_canonical)] = activator
                    elif sbol.SBO_STIMULATED in parti.roles:
                        activated_fc = circuit_def.functionalComponents.get(parti.participant)
                        activated = CircuitSpecies(activated_fc.identity, activated_fc.definition)

                if len(activator_dict) > 0 and activated:
                    activated_identity = activated.get_identity(is_canonical)

                    if activated_identity not in activation_map:
                        activation_map[activated_identity] = []

                    activation_map[activated_identity].append(
                                                                [
                                                                    activator_dict[activator_identity]
                                                                    for activator_identity in activator_dict
                                                                ]
                                                             )

                    self.__species_dict[activated_identity] = activated
                    self.__species_dict.update(activator_dict)

        for sub_circuit in circuit_def.modules:
            sub_circuit_def = circuit_doc.moduleDefinitions.get(sub_circuit.definition)

            self.__build_activation_map(sub_circuit_def, circuit_doc, activation_map, is_canonical)

        return activation_map

    def __build_repression_map(self, circuit_def, circuit_doc, repression_map={}, is_canonical=True):
        for intxn in circuit_def.interactions:
            if sbol.SBO_INHIBITION in intxn.types:
                repressor_dict = {}
                repressed = None

                for parti in intxn.participations:
                    if sbol.SBO_INHIBITOR in parti.roles:
                        repressor_fc = circuit_def.functionalComponents.get(parti.participant)
                        repressor = CircuitSpecies(repressor_fc.identity, repressor_fc.definition)
                        repressor_dict[repressor.get_identity(is_canonical)] = repressor
                    elif sbol.SBO_INHIBITED in parti.roles:
                        repressed_fc = circuit_def.functionalComponents.get(parti.participant)
                        repressed = CircuitSpecies(repressed_fc.identity, repressed_fc.definition)

                if len(repressor_dict) > 0 and repressed:
                    repressed_identity = repressed.get_identity(is_canonical)

                    if repressed_identity not in repression_map:
                        repression_map[repressed_identity] = []

                    repression_map[repressed_identity].append(
                                                                [
                                                                    repressor_dict[repressor_identity]
                                                                    for repressor_identity in repressor_dict
                                                                ]
                                                             )

                    self.__species_dict[repressed_identity] = repressed
                    self.__species_dict.update(repressor_dict)

        for sub_circuit in circuit_def.modules:
            sub_circuit_def = circuit_doc.moduleDefinitions.get(sub_circuit.definition)

            self.__build_repression_map(sub_circuit_def, circuit_doc, repression_map, is_canonical)

        return repression_map

    def __distill_transcription_map(self, circuit_def, circuit_doc, activation_map, transcription_map={}, is_canonical=True):
        deleted_keys = []

        for activated_identity in activation_map:
            for activators in activation_map[activated_identity]:
                if len(activators) == 1:
                    activated_def = circuit_doc.componentDefinitions.get(self.__species_dict[activated_identity].canonical_identity)
                    activator_def = circuit_doc.componentDefinitions.get(activators[0].canonical_identity)

                    if (sbol.BIOPAX_DNA in activated_def.types and sbol.SO_CDS in activated_def.roles
                            and sbol.BIOPAX_DNA in activator_def.types and sbol.SO_PROMOTER in activator_def.roles):
                        if activated_identity not in transcription_map:
                            transcription_map[activated_identity] = []

                        transcription_map[activated_identity].append(activators[0])

                        deleted_keys.append(activated_identity)

        for deleted_key in deleted_keys:
            del activation_map[deleted_key]

        for sub_circuit in circuit_def.modules:
            sub_circuit_def = circuit_doc.moduleDefinitions.get(sub_circuit.definition)

            self.__distill_transcription_map(sub_circuit_def, circuit_doc, activation_map, transcription_map, is_canonical)

        return transcription_map

    def compute_truth_table(self, is_canonical=True):
        truth_table = []

        circuit_input_ranges = [[0, 1] for i in range(0, len(self.__circuit_inputs))]

        for circuit_input_values in product(*circuit_input_ranges):
            truth_table.append(
                                {
                                    self.__circuit_inputs[i].get_identity(is_canonical) : circuit_input_values[i] 
                                    for i in range(0, len(self.__circuit_inputs))
                                }
                              )

        for truth_row in truth_table:
            for circuit_output in self.__circuit_outputs:
                circuit_output_identity = circuit_output.get_identity(is_canonical)

                if circuit_output_identity in self.__product_to_gates:
                    for output_gate in self.__product_to_gates[circuit_output_identity]:
                        truth_row[circuit_output_identity] = self.__calculate_row_value(output_gate, truth_row)

        return truth_table

    def __calculate_row_value(self, gate, truth_row, is_canonical=True):
        for gate_input in gate.get_inputs():
            gate_input_identity = gate_input.get_identity(is_canonical)

            if gate_input_identity not in truth_row and gate_input_identity in self.__product_to_gates:
                input_gates = self.__product_to_gates[gate_input_identity]

                if len(input_gates) == 1:
                    truth_row[gate_input_identity] = self.__calculate_row_value(input_gates[0], truth_row, is_canonical)
                else:
                    truth_row[gate_input_identity] = 0

                    for input_gate in input_gates:
                        if self.__calculate_row_value(input_gate, truth_row, is_canonical) == 1:
                            truth_row[gate_input_identity] = 1

        return gate.calculate_output_value(truth_row, is_canonical)

def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--namespace')
    parser.add_argument('-t', '--target_files', nargs='*', default=[])
    parser.add_argument('-o', '--output_files', nargs='*', default=[])
    parser.add_argument('-l', '--log_file', nargs='?', default='')
    parser.add_argument('-v', '--validate', action='store_true')
    parser.add_argument('-s', '--table_suffix', nargs='?', default='')
    
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

    for i in range (0, len(args.target_files)):
        target_doc = load_sbol(args.target_files[i])

        sub_network_identities = set()
        network_identities = set()

        for mod_def in target_doc.moduleDefinitions:
            if 'http://purl.obolibrary.org/obo/NCIT_C20633' in mod_def.roles:
                network_identities.add(mod_def.identity)

                for mod in mod_def.modules:
                    sub_network_identities.add(mod.definition)

        network_identities.difference_update(sub_network_identities)

        for network_identity in network_identities:
            logic_circuit = LogicCircuit(target_doc.moduleDefinitions.get(network_identity), target_doc)

            truth_table = logic_circuit.compute_truth_table()

            print(logic_circuit.serialize())

            for csv_row in logic_circuit.serialize_truth_table(truth_table):
                print(csv_row)
        
        # if len(args.output_files) == 1 and os.path.isdir(args.output_files[0]):
        #     (target_file_path, target_filename) = os.path.split(target_files[i])
        #     (target_file_base, target_file_extension) = os.path.splitext(target_filename)

        #     if len(args.circuit_suffix) > 0:
        #         output_file = os.path.join(args.output_files[0], '_'.join([target_file_base, args.circuit_suffix + target_file_extension]))
        #     else:
        #         output_file = os.path.join(args.output_files[0], target_file_base + target_file_extension)

        # elif i < len(args.output_files):
        #     output_file = args.output_files[i]
        # else:
        #     (target_file_base, target_file_extension) = os.path.splitext(target_files[i])

        #     if len(args.circuit_suffix) > 0:
        #         output_file = '_'.join([target_file_base, args.circuit_suffix + target_file_extension])
        #     else:
        #         output_file = target_file_base + target_file_extension

        #     if sbol.Config.getOption('validate') == True:
        #         logging.info('Validating and writing %s', output_file)
        #     else:
        #         logging.info('Writing %s', output_file)

        # target_doc.write(output_file)

    #     logging.info('Finished writing %s', output_file)

    logging.info('Finished generating truth tables')


# class CircuitIdentityError(Exception):

#     def __init__(self, circuit_identity, feature_ID):
#         self.circuit_identity = circuit_identity
#         self.feature_ID = feature_ID

#     def __str__(self):
#         return "Circuit identity {ci} does not contain the feature {fi}.".format(ci=self.circuit_identity, fi=self.feature_ID)


if __name__ == '__main__':
    main()
