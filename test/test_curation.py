import unittest
import os
from sbol import *

from sequences_to_features import *
from features_to_circuits import *

class CurationTests(unittest.TestCase):

    def test_curate_nand_circuit(self):
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

        HOMESPACE = 'http://synbict.org'
        VERSION = '1'

        MIN_TARGET_LENGTH = 2000

        setHomespace(HOMESPACE)
        Config.setOption('sbol_typed_uris', False)

        cello_doc = load_sbol(os.path.join(__location__, 'cello_library.xml'))
        feature_library = FeatureLibrary([cello_doc])

        target_doc = load_sbol(os.path.join(__location__, 'genetic_nand.xml'))
        target_construct_library = FeatureLibrary([target_doc], True)

        feature_annotater = FeatureAnnotater(feature_library, 40)
        feature_annotater.annotate(target_doc, target_construct_library.features, MIN_TARGET_LENGTH, VERSION)

        feature_pruner = FeaturePruner(feature_library)
        feature_pruner.prune(target_doc, target_construct_library.features, 14, MIN_TARGET_LENGTH, False, feature_library)

        circuit_library = CircuitLibrary([cello_doc])
        target_device_library = FeatureLibrary([target_doc], require_sequence=False)

        circuit_ID = 'nand_circuit'

        circuit_builder = CircuitBuilder(circuit_library)
        circuit_builder.build(circuit_ID, target_doc, target_device_library, MIN_TARGET_LENGTH)

        nand_circuit = target_doc.getModuleDefinition('/'.join([HOMESPACE, circuit_ID]))

        nand_devices = {
            '/'.join([HOMESPACE, 'Strain_4_MG1655_Genomic_NAND_Circuit', VERSION])
        }

        devices = set()

        for device in nand_circuit.functionalComponents:
            devices.add(device.definition)

        device_intersection = nand_devices.intersection(devices)
        device_difference = nand_devices.difference(devices)

        nand_sub_circuits = {
            'https://synbiohub.programmingbiology.org/public/Cello_Parts/YFP_protein_production/1',
            'https://synbiohub.programmingbiology.org/public/Cello_Parts/PhlF_pPhlF_repression/1',
            'https://synbiohub.programmingbiology.org/public/Cello_Parts/IcaRA_protein_production/1',
            'https://synbiohub.programmingbiology.org/public/Cello_Parts/PhlF_protein_production/1',
            'https://synbiohub.programmingbiology.org/public/Cello_Parts/LacI_pTac_repression/1',
            'https://synbiohub.programmingbiology.org/public/Cello_Parts/TetR_protein_production/1',
            'https://synbiohub.programmingbiology.org/public/Cello_Parts/IcaRA_pIcaRA_repression/1',
            'https://synbiohub.programmingbiology.org/public/Cello_Parts/LacI_protein_production/1'
        }

        sub_circuits = set()

        for sub_circuit in nand_circuit.modules:
            sub_circuits.add(sub_circuit.definition)

        sub_circuit_intersection = nand_sub_circuits.intersection(sub_circuits)
        sub_circuit_difference = nand_sub_circuits.difference(sub_circuits)

        self.assertEqual(len(device_intersection), len(nand_devices),
            "Inferred circuit is missing expected devices: {mi}".format(mi=', '.join(device_difference)))
        self.assertEqual(len(sub_circuit_intersection), len(nand_sub_circuits),
            "Inferred circuit is missing expected sub-circuits: {mi}".format(mi=', '.join(sub_circuit_difference)))

if __name__ == '__main__':
    unittest.main()