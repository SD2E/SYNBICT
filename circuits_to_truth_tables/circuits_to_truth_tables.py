import logging
import argparse
import os
import sys
import csv
import json
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

    NOT = 'https://www.openmath.org/cd/logic1#not'
    NOR = 'https://www.openmath.org/cd/logic1#nor'
    XNOR = 'https://www.openmath.org/cd/logic1#xnor'
    AND = 'https://www.openmath.org/cd/logic1#and'
    OR = 'https://www.openmath.org/cd/logic1#or'
    XOR = 'https://www.openmath.org/cd/logic1#xor'
    NAND = 'https://www.openmath.org/cd/logic1#nand'
    YES = 'http://synbict.org/logic#yes'
    UNDEFINED = 'http://synbict.org/logic#logic_function'

    GATE_TYPES = {NOT, NOR, XNOR, AND, OR, XOR, NAND, YES, UNDEFINED}

    def __init__(self, gate_output, gate_inputs, gate_type):
        if len(gate_inputs) == 0:
            pass
        elif len(gate_inputs) > 2:
            pass

        self.__gate_output = gate_output
        self.__gate_inputs = gate_inputs

        if gate_type in self.GATE_TYPES:
            self.__gate_type = gate_type
        else:
            raise GateTypeError(gate_type)

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
                return ~((input_1_value & input_2_value) - 2)
            else:
                return -1
        else:
            return -1

    def get_inputs(self):
        return self.__gate_inputs

    def get_output(self):
        return self.__gate_output

    def serialize_logic(self, display_ID=False):
        if display_ID:
            input_identities = [
                                 gi.get_identity().split('/')[-2]
                                 for gi in self.__gate_inputs
                                 if len(gi.get_identity().split('/')) > 1
                                ]
        else:
            input_identities = [
                                 gi.get_identity()
                                 for gi in self.__gate_inputs
                                ]

        if display_ID:
            if len(self.__gate_output.get_identity().split('/')) > 1:
                output_identity = self.__gate_output.get_identity().split('/')[-2]
            else:
                output_identity = ''
        else:
            output_identity = self.__gate_output.get_identity()

        if len(input_identities) == len(self.__gate_inputs) and len(output_identity) > 0:
            if len(input_identities) == 1:
                if self.__gate_type == self.YES:
                    # return output_identity + '=(' + input_identities[0] + ')'
                    return input_identities[0]
                elif self.__gate_type == self.NOT:
                    # return output_identity + '=~(' + input_identities[0] + ')'
                    return '~' + input_identities[0]
                else:
                    return '-'
            elif len(input_identities) == 2:
                if self.__gate_type == self.OR:
                    # return output_identity + '=(' + '|'.join([ii for ii in input_identities]) + ')'
                    return '|'.join([ii for ii in input_identities])
                elif self.__gate_type == self.NOR:
                    # return output_identity + '=~(' + '|'.join([ii for ii in input_identities]) + ')'
                    return '~(' + '|'.join([ii for ii in input_identities]) + ')'
                elif self.__gate_type == self.AND:
                    # return output_identity + '=(' + '&'.join([ii for ii in input_identities]) + ')'
                    return '&'.join([ii for ii in input_identities])
                elif self.__gate_type == self.NAND:
                    # return output_identity + '=~(' + '&'.join([ii for ii in input_identities]) + ')'
                    return '~(' + '&'.join([ii for ii in input_identities])
                else:
                    return '-'
        else:
            return '-' 

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

    def __init__(self, circuit_def=None, circuit_doc=None, is_canonical=True, infer_io=False, open_math=''):
        self.__circuit_inputs = []
        self.__circuit_outputs = []
        self.__circuit_input_intermediates = []
        self.__circuit_output_intermediates = []
        self.__species_dict = {}

        if circuit_def and circuit_doc:
            production_map = self.__build_production_map(circuit_def, circuit_doc, is_canonical)
            (repression_map, repression_logic_map) = self.__build_repression_maps(circuit_def, circuit_doc, is_canonical, open_math)
            (activation_map, activation_logic_map) = self.__build_activation_maps(circuit_def, circuit_doc, is_canonical, open_math)

            transcription_map = self.__distill_transcription_map(circuit_def, circuit_doc, activation_map, is_canonical)
            positive_induction_map = self.__distill_induction_map(circuit_def, circuit_doc, activation_map, is_canonical)
            negative_induction_map = self.__distill_induction_map(circuit_def, circuit_doc, repression_map, is_canonical)

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

            self.__product_to_gates = self.__build_gate_map(production_map,
                                                            repression_map,
                                                            activation_map,
                                                            transcription_map,
                                                            positive_induction_map,
                                                            negative_induction_map,
                                                            repression_logic_map,
                                                            activation_logic_map,
                                                            is_canonical,
                                                            open_math)

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
                inferred_inputs = [
                                    gate_input_dict[gate_input_identity]
                                    for gate_input_identity in gate_input_dict
                                    if gate_input_identity not in gate_output_dict
                                  ]

                if infer_io:
                    self.__circuit_inputs.extend(inferred_inputs)
                else:
                    self.__circuit_input_intermediates.extend(inferred_inputs)

            if len(self.__circuit_outputs) == 0:
                inferred_outputs = [
                                    gate_output_dict[gate_output_identity]
                                    for gate_output_identity in gate_output_dict
                                    if gate_output_identity not in gate_input_dict
                                   ]

                if infer_io:
                    self.__circuit_outputs.extend(inferred_outputs)
                else:
                    self.__circuit_output_intermediates.extend(inferred_outputs)

    def is_complete(self):
        return ((len(self.__circuit_inputs) > 0 or len(self.__circuit_input_intermediates) > 0)
                and (len(self.__circuit_outputs) > 0 or len(self.__circuit_output_intermediates) > 0))

    def get_intermediates(self, is_canonical=True):
        circuit_input_identities = {
                                        circuit_input.get_identity(is_canonical)
                                        for circuit_input in self.__circuit_inputs
                                   }
        circuit_input_identities.update(
                                            {
                                                cii.get_identity(is_canonical)
                                                for cii in self.__circuit_input_intermediates
                                            }
                                        )

        circuit_output_identities = {
                                        circuit_output.get_identity(is_canonical)
                                        for circuit_output in self.__circuit_outputs
                                    }
        circuit_output_identities.update(
                                            {
                                                cio.get_identity(is_canonical)
                                                for cio in self.__circuit_output_intermediates
                                            }
                                        )

        circuit_intermediates = [
                                    self.__species_dict[product_identity]
                                    for product_identity in self.__product_to_gates
                                    if (
                                        product_identity not in circuit_input_identities
                                        and product_identity not in circuit_output_identities
                                    )
                                ]

        return circuit_intermediates

    def serialize_truth_table_to_json(self, truth_table, is_canonical=True, verbose=False):
        if len(self.__circuit_input_intermediates) > 0:
            table_inputs = []

            table_intermediates = [
                                    cii
                                    for cii in self.__circuit_input_intermediates
                                  ]
        else:
            table_inputs = [
                            ci
                            for ci in self.__circuit_inputs
                           ]

            table_intermediates = []

        table_intermediates.extend(self.get_intermediates(is_canonical))

        if len(self.__circuit_output_intermediates) > 0:
            table_outputs = []

            table_intermediates.extend([
                                         cio
                                         for cio in self.__circuit_output_intermediates
                                       ])
        else:
            table_outputs = [
                             co
                             for co in self.__circuit_outputs
                            ]

        if verbose:                
            truth_json = {
                          'inputs': len(table_inputs),
                          'outputs': len(table_outputs),
                          'variable_map': {},
                          'logic_map': {},
                          'variable_display_map' : {},
                          'logic_debug_map': {}
                         }
        else:
            truth_json = {
                          'inputs': len(table_inputs),
                          'outputs': len(table_outputs),
                          'variable_map': {}
                         }

        truth_json.update({
                            'input_column_' + str(i) : []
                            for i in range(0, len(table_inputs))
                          })

        truth_json['variable_map'].update({
                                            table_inputs[i].get_identity(is_canonical) : 'input_column_' + str(i)
                                            for i in range(0, len(table_inputs))
                                          })

        if verbose:
            truth_json['variable_display_map'].update({
                                                        'input_column_' + str(i) : table_inputs[i].get_identity(is_canonical).split('/')[-2]
                                                        for i in range(0, len(table_inputs))
                                                      })

        if verbose:
            truth_json.update({
                                'intermediate_column_' + str(i) : []
                                for i in range(0, len(table_intermediates))
                              })

            truth_json['variable_map'].update({
                                                table_intermediates[i].get_identity(is_canonical) : 'intermediate_column_' + str(i)
                                                for i in range(0, len(table_intermediates))
                                              })

            truth_json['variable_display_map'].update({
                                                        'intermediate_column_' + str(i) : table_intermediates[i].get_identity(is_canonical).split('/')[-2]
                                                        for i in range(0, len(table_intermediates))
                                                      })

        truth_json.update({
                            'output_column_' + str(i) : []
                            for i in range(0, len(table_outputs))
                          })

        truth_json['variable_map'].update({
                                            table_outputs[i].get_identity(is_canonical) : 'output_column_' + str(i)
                                            for i in range(0, len(table_outputs))
                                          })

        if verbose:
            truth_json['variable_display_map'].update({
                                                        'output_column_' + str(i) : table_outputs[i].get_identity(is_canonical).split('/')[-2]
                                                        for i in range(0, len(table_outputs))
                                                      })

            for table_input in table_inputs:
                input_identity = table_input.get_identity(is_canonical)

                truth_json['logic_map'][input_identity] = input_identity

                truth_json['logic_debug_map'][input_identity.split('/')[-2]] = input_identity.split('/')[-2]

            for table_intermediate in table_intermediates:
                intermediate_identity = table_intermediate.get_identity(is_canonical)

                if intermediate_identity in self.__product_to_gates:
                    truth_json['logic_map'][intermediate_identity] = '|'.join([
                                                                                ga.serialize_logic()
                                                                                for ga in self.__product_to_gates[intermediate_identity]
                                                                              ])

                    truth_json['logic_debug_map'][intermediate_identity.split('/')[-2]] = '|'.join([
                                                                                                     ga.serialize_logic(True)
                                                                                                     for ga in self.__product_to_gates[intermediate_identity]
                                                                                                   ])

            for table_output in table_outputs:
                output_identity = table_output.get_identity(is_canonical)

                if output_identity in self.__product_to_gates:
                    truth_json['logic_map'][output_identity] = '|'.join([
                                                                          ga.serialize_logic()
                                                                          for ga in self.__product_to_gates[output_identity]
                                                                        ])

                    truth_json['logic_debug_map'][output_identity.split('/')[-2]] = '|'.join([
                                                                                               ga.serialize_logic(True)
                                                                                               for ga in self.__product_to_gates[output_identity]
                                                                                             ])

        for table_input in table_inputs:
            input_identity = table_input.get_identity(is_canonical)

            truth_json[truth_json['variable_map'][input_identity]] = [
                                                                    truth_row[input_identity]
                                                                    for truth_row in truth_table
                                                                  ]

        if verbose:
            for table_intermediate in table_intermediates:
                intermediate_identity = table_intermediate.get_identity(is_canonical)

                truth_json[truth_json['variable_map'][intermediate_identity]] = [
                                                                        truth_row[intermediate_identity]
                                                                        for truth_row in truth_table
                                                                      ]

        for table_output in table_outputs:
            output_identity = table_output.get_identity(is_canonical)

            truth_json[truth_json['variable_map'][output_identity]] = [
                                                                    truth_row[output_identity]
                                                                    for truth_row in truth_table
                                                                  ]

        return truth_json

    def serialize_truth_table_to_csv(self, truth_table, is_canonical=True, verbose=False):
        truth_csv = []

        if len(self.__circuit_input_intermediates) > 0:
            table_inputs = []

            table_intermediates = [
                                    cii
                                    for cii in self.__circuit_input_intermediates
                                  ]
        else:
            table_inputs = [
                            ci
                            for ci in self.__circuit_inputs
                           ]

            table_intermediates = []

        table_intermediates.extend(self.get_intermediates(is_canonical))

        if len(self.__circuit_output_intermediates) > 0:
            table_outputs = []

            table_intermediates.extend([
                                         cio
                                         for cio in self.__circuit_output_intermediates
                                       ])
        else:
            table_outputs = [
                             co
                             for co in self.__circuit_outputs
                            ]

        
        table_input_headers = [
                                'input_column_' + str(i)
                                for i in range(0, len(table_inputs))
                              ]

        table_intermediate_headers = [
                                        'intermediate_column_' + str(i)
                                        for i in range(0, len(table_intermediates))
                                     ]

        table_output_headers = [
                                'output_column_' + str(i)
                                for i in range(0, len(table_outputs))
                               ]

        if verbose:
            truth_csv.append(table_input_headers + table_intermediate_headers + table_output_headers)
        
            circuit_species = table_inputs + table_intermediates + table_outputs
        else:
            truth_csv.append(table_input_headers + table_output_headers)

            circuit_species = table_inputs + table_outputs

        species_identities = [
                                cs.get_identity(is_canonical)
                                for cs in circuit_species
                             ]

        truth_csv.append(species_identities)

        species_logic = []

        for i in range(0, len(species_identities)):
            if species_identities[i] in self.__product_to_gates:
                species_logic.append('|'.join([ga.serialize_logic()
                                               for ga in self.__product_to_gates[species_identities[i]]]))
            else:
                species_logic.append(species_identities[i])

        truth_csv.append(species_logic)

        species_IDs = [
                        si.split('/')[-2]
                        for si in species_identities
                        if len(si.split('/')) > 1
                      ]

        if len(species_IDs) == len(species_identities):
            truth_csv.append(species_IDs)

            species_ID_logic = []

            for i in range(0, len(species_identities)):
                if species_identities[i] in self.__product_to_gates:
                    species_ID_logic.append('|'.join([ga.serialize_logic(True)
                                                      for ga in self.__product_to_gates[species_identities[i]]]))
                else:
                    species_ID_logic.append(species_identities[i].split('/')[-2])

            truth_csv.append(species_ID_logic)

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

    def __build_gate_map(self, production_map, repression_map, activation_map, transcription_map,
            positive_induction_map, negative_induction_map, repression_logic_map, activation_logic_map,
            is_canonical=True, open_math=''):
        product_to_gates = {}

        if len(production_map) == 0:
            for activated_identity in activation_map:
                if activated_identity not in product_to_gates:
                    product_to_gates[activated_identity] = []

                activated = self.__species_dict[activated_identity]

                for i in range(0, len(activation_map[activated_identity])):
                    activators = activation_map[activated_identity][i]

                    if activated_identity in activation_logic_map and i < len(activation_logic_map[activated_identity]):
                        gate_type = activation_logic_map[activated_identity][i]
                    else:
                        gate_type = LogicGate.UNDEFINED

                    if gate_type == LogicGate.UNDEFINED:
                        gate_type = self.__infer_gate_type(len(activators), 0, True, False, len(activators) > 1)

                    product_to_gates[activated_identity].append(LogicGate(activated, activators, gate_type))

            for repressed_identity in repression_map:
                if repressed_identity not in product_to_gates:
                    product_to_gates[repressed_identity] = []

                repressed = self.__species_dict[repressed_identity]

                for i in range(0, len(repression_map[repressed_identity])):
                    repressors = repression_map[repressed_identity][i]

                    if repressed_identity in repression_logic_map and i < len(repression_logic_map[repressed_identity]):
                        gate_type = repression_logic_map[repressed_identity][i]
                    else:
                        gate_type = LogicGate.UNDEFINED

                    if gate_type == LogicGate.UNDEFINED:
                        gate_type = self.__infer_gate_type(len(repressors), 0, False, True, len(repressors) > 1)

                    product_to_gates[repressed_identity].append(LogicGate(repressed, repressors, gate_type))
        else:
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
                                    if len(activators) > 1:
                                        cooperativity_flag = True

                                    for activator in activators:
                                        if activator.get_identity(is_canonical) in negative_induction_map:
                                            inducers = negative_induction_map[activator.get_identity(is_canonical)]

                                            for inducer in inducers:
                                                gate_input_dict[inducer.get_identity(is_canonical)] = inducer

                                            repression_flag = True
                                        elif activator.get_identity(is_canonical) in positive_induction_map:
                                            inducers = positive_induction_map[activator.get_identity(is_canonical)]

                                            for inducer in inducers:
                                                gate_input_dict[inducer.get_identity(is_canonical)] = inducer

                                            activation_flag = True
                                        else:
                                            gate_input_dict[activator.get_identity(is_canonical)] = activator

                                            activation_flag = True
                            elif promoter_identity in repression_map:
                                for repressors in repression_map[promoter_identity]:
                                    if len(repressors) > 1:
                                        cooperativity_flag = True

                                    for repressor in repressors:
                                        if repressor.get_identity(is_canonical) in negative_induction_map:
                                            inducers = negative_induction_map[repressor.get_identity(is_canonical)]

                                            for inducer in inducers:
                                                gate_input_dict[inducer.get_identity(is_canonical)] = inducer

                                            activation_flag = True
                                        elif repressor.get_identity(is_canonical) in positive_induction_map:
                                            inducers = positive_induction_map[activator.get_identity(is_canonical)]

                                            for inducer in inducers:
                                                gate_input_dict[inducer.get_identity(is_canonical)] = inducer

                                            repression_flag = True
                                        else:
                                            gate_input_dict[repressor.get_identity(is_canonical)] = repressor

                                            repression_flag = True
                    elif template_identity in activation_map:
                        for activators in activation_map[template_identity]:
                            if len(activators) > 1:
                                cooperativity_flag = True

                            for activator in activators:
                                gate_input_dict[activator.get_identity(is_canonical)] = activator

                        activation_flag = True
                    elif template_identity in repression_map:
                        for repressors in repression_map[template_identity]:
                            if len(repressors) > 1:
                                cooperativity_flag = True

                            for repressor in repressors:
                                gate_input_dict[repressor.get_identity(is_canonical)] = repressor

                        repression_flag = True

                    gate_type = self.__infer_gate_type(len(gate_input_dict), len(promoters), repression_flag,
                                                       activation_flag, cooperativity_flag)

                    if gate_type != LogicGate.UNDEFINED:
                        gate_inputs = [gate_input_dict[gate_input_identity] for gate_input_identity in gate_input_dict]

                        if product_identity not in product_to_gates:
                            product_to_gates[product_identity] = []

                        product = self.__species_dict[product_identity]

                        product_to_gates[product_identity].append(LogicGate(product, gate_inputs, gate_type))

        return product_to_gates

    @classmethod
    def __infer_gate_type(cls, num_inputs, num_promoters, is_repression, is_activation, is_cooperative):
        if num_inputs == 0:
            # print('No inputs')

            return LogicGate.UNDEFINED
        elif num_inputs > 2:
            # print('Too many inputs')

            return LogicGate.UNDEFINED
        elif not is_repression and not is_activation:
            # print('No regulation')

            return LogicGate.UNDEFINED
        elif is_repression and is_activation:
            # print('Clashing regulation')

            return LogicGate.UNDEFINED
        elif num_promoters > 2:
            # print('Too many promoters')

            return LogicGate.UNDEFINED
        elif num_inputs == 1:
            if is_activation:
                return LogicGate.YES
            else:
                return LogicGate.NOT
        elif is_activation:
                if num_promoters < 2 and is_cooperative:
                    return LogicGate.AND
                else:
                    return LogicGate.OR
        elif num_promoters < 2 and not is_cooperative:
            return LogicGate.NOR
        else:
            return LogicGate.NAND


    def __build_production_map(self, circuit_def, circuit_doc, is_canonical=True):
        production_map = {}

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

            production_map.update(self.__build_production_map(sub_circuit_def, circuit_doc, is_canonical))

        return production_map

    def __build_activation_maps(self, circuit_def, circuit_doc, is_canonical=True, open_math=''):
        activation_map = {}
        activation_logic_map = {}

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

                    if len(open_math) > 0:
                        try:
                            activation_logic = intxn.getPropertyValue(open_math)
                        except LookupError:
                            activation_logic = LogicGate.UNDEFINED

                        if activated_identity not in activation_logic_map:
                            activation_logic_map[activated_identity] = []

                        activation_logic_map[activated_identity].append(activation_logic)

                    self.__species_dict[activated_identity] = activated
                    self.__species_dict.update(activator_dict)

        for sub_circuit in circuit_def.modules:
            sub_circuit_def = circuit_doc.moduleDefinitions.get(sub_circuit.definition)

            (sub_activation_map, sub_activation_logic_map) = self.__build_activation_maps(sub_circuit_def, circuit_doc, is_canonical)

            activation_map.update(sub_activation_map)
            sub_activation_map.update(sub_activation_logic_map)

        return (activation_map, activation_logic_map)

    def __build_repression_maps(self, circuit_def, circuit_doc, is_canonical=True, open_math=''):
        repression_map = {}
        repression_logic_map = {}

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

                    if len(open_math) > 0:
                        try:
                            repression_logic = intxn.getPropertyValue(open_math)
                        except LookupError:
                            repression_logic = LogicGate.UNDEFINED

                        if repressed_identity not in repression_logic_map:
                            repression_logic_map[repressed_identity] = []

                        repression_logic_map[repressed_identity].append(repression_logic)

                    self.__species_dict[repressed_identity] = repressed
                    self.__species_dict.update(repressor_dict)

        for sub_circuit in circuit_def.modules:
            sub_circuit_def = circuit_doc.moduleDefinitions.get(sub_circuit.definition)

            (sub_repression_map, sub_repression_logic_map) = self.__build_repression_maps(sub_circuit_def, circuit_doc, is_canonical)

            repression_map.update(sub_repression_map)
            repression_logic_map.update(sub_repression_logic_map)

        return (repression_map, repression_logic_map)

    def __distill_transcription_map(self, circuit_def, circuit_doc, activation_map, is_canonical=True):
        transcription_map = {}

        deleted_keys = set()

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

                        deleted_keys.add(activated_identity)

        for deleted_key in deleted_keys:
            del activation_map[deleted_key]

        for sub_circuit in circuit_def.modules:
            sub_circuit_def = circuit_doc.moduleDefinitions.get(sub_circuit.definition)

            transcription_map.update(self.__distill_transcription_map(sub_circuit_def, circuit_doc, activation_map, is_canonical))

        return transcription_map

    def __distill_induction_map(self, circuit_def, circuit_doc, regulation_map, is_canonical=True):
        induction_map = {}

        deleted_keys = []

        for regulated_identity in regulation_map:
            for regulators in regulation_map[regulated_identity]:
                if len(regulators) == 1:
                    regulated_def = circuit_doc.componentDefinitions.get(self.__species_dict[regulated_identity].canonical_identity)
                    regulator_def = circuit_doc.componentDefinitions.get(regulators[0].canonical_identity)

                    if sbol.BIOPAX_PROTEIN in regulated_def.types and sbol.BIOPAX_SMALL_MOLECULE in regulator_def.types:
                        if regulated_identity not in induction_map:
                            induction_map[regulated_identity] = []

                        induction_map[regulated_identity].append(regulators[0])

                        deleted_keys.append(regulated_identity)

        for deleted_key in deleted_keys:
            del regulation_map[deleted_key]

        for sub_circuit in circuit_def.modules:
            sub_circuit_def = circuit_doc.moduleDefinitions.get(sub_circuit.definition)

            induction_map.update(self.__distill_induction_map(sub_circuit_def, circuit_doc, regulation_map, is_canonical))

        return induction_map

    def compute_truth_table(self, is_canonical=True):
        truth_table = []

        if len(self.__circuit_inputs) > 0:
            table_inputs = [ci for ci in self.__circuit_inputs]
        else:
            table_inputs = [ci for ci in self.__circuit_input_intermediates]

        table_input_ranges = [[0, 1] for i in range(0, len(table_inputs))]

        for table_input_values in product(*table_input_ranges):
            truth_table.append(
                                {
                                    table_inputs[i].get_identity(is_canonical) : table_input_values[i] 
                                    for i in range(0, len(table_inputs))
                                }
                              )

        if len(self.__circuit_outputs) > 0:
            table_outputs = [co for co in self.__circuit_outputs]
        else:
            table_outputs = [co for co in self.__circuit_output_intermediates]

        for truth_row in truth_table:
            for table_output in table_outputs:
                table_output_identity = table_output.get_identity(is_canonical)

                if table_output_identity in self.__product_to_gates:
                    for output_gate in self.__product_to_gates[table_output_identity]:
                        truth_row[table_output_identity] = self.__calculate_row_value(output_gate, truth_row)

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
    parser.add_argument('-io', '--infer_io', action='store_true')
    parser.add_argument('-m', '--open_math', nargs='?', default='')
    parser.add_argument('-c', '--csv', action='store_true')
    parser.add_argument('-j', '--json', action='store_true')
    parser.add_argument('-b', '--verbose', action='store_true')
    
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

    target_files = []
    for target_file in args.target_files:
        if os.path.isdir(target_file):
            target_files.extend([os.path.join(target_file, tf) for tf in os.listdir(target_file) if
                                 os.path.isfile(os.path.join(target_file, tf)) and tf.endswith('.xml')])
        else:
            target_files.append(target_file)

    for i in range (0, len(target_files)):
        target_doc = load_sbol(target_files[i])

        sub_network_identities = set()
        network_identities = set()

        for mod_def in target_doc.moduleDefinitions:
            if 'http://purl.obolibrary.org/obo/NCIT_C20633' in mod_def.roles:
                network_identities.add(mod_def.identity)

                for mod in mod_def.modules:
                    sub_network_identities.add(mod.definition)

        network_identities.difference_update(sub_network_identities)

        if len(network_identities) == 0:
            logging.warning('Cannot generate truth tables for %s since it lacks root ModuleDefinitions with role "biochemical pathway" (URI for NCIT term is http://purl.obolibrary.org/obo/NCIT_C20633).', target_files[i])
        else:
            for network_identity in network_identities:
                logic_circuit = LogicCircuit(target_doc.moduleDefinitions.get(network_identity), target_doc, infer_io=args.infer_io,
                                             open_math=args.open_math)

                if logic_circuit.is_complete():
                    truth_table = logic_circuit.compute_truth_table()

                    # print(logic_circuit.serialize())

                    # for csv_row in logic_circuit.serialize_truth_table_to_csv(truth_table):
                        # print(csv_row)
                
                    if len(args.output_files) == 1 and os.path.isdir(args.output_files[0]):
                        (target_file_path, target_filename) = os.path.split(target_files[i])
                        (target_file_base, target_file_extension) = os.path.splitext(target_filename)

                        if len(args.table_suffix) > 0:
                            output_file = os.path.join(args.output_files[0], '_'.join([target_file_base, args.table_suffix]))
                        else:
                            output_file = os.path.join(args.output_files[0], target_file_base)

                    elif i < len(args.output_files):
                        output_file = args.output_files[i]
                    else:
                        (target_file_base, target_file_extension) = os.path.splitext(target_files[i])

                        if len(args.table_suffix) > 0:
                            output_file = '_'.join([target_file_base, args.table_suffix])
                        else:
                            output_file = target_file_base

                    logging.info('Writing %s', output_file)

                    if args.csv or not args.json:
                        with open(output_file + '.csv', 'w', newline='\n') as csv_handle:
                            tt_writer = csv.writer(csv_handle, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

                            for csv_row in logic_circuit.serialize_truth_table_to_csv(truth_table, verbose=args.verbose):
                                tt_writer.writerow(csv_row)

                            logging.info('Finished writing %s', output_file + '.csv')

                    if args.json:
                        with open(output_file + '.json', 'w', newline='\n') as json_handle:
                            json.dump(logic_circuit.serialize_truth_table_to_json(truth_table, verbose=args.verbose), json_handle, indent=4)

                        logging.info('Finished writing %s', output_file + '.json')
                else:
                    logging.warning('Cannot generate truth table for %s since it lacks inputs and/or outputs.', network_identity)

    logging.info('Finished generating truth tables')


# class CircuitIdentityError(Exception):

#     def __init__(self, circuit_identity, feature_ID):
#         self.circuit_identity = circuit_identity
#         self.feature_ID = feature_ID

#     def __str__(self):
#         return "Circuit identity {ci} does not contain the feature {fi}.".format(ci=self.circuit_identity, fi=self.feature_ID)


class GateTypeError(Exception):

    def __init__(self, gate_type):
        self.gate_type = gate_type

    def __str__(self):
        return "Gate type {gt} is not recognized as a URI for a valid OpenMath term or SYNBICT term.".format(gt=self.gate_type)


if __name__ == '__main__':
    main()
