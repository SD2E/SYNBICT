# SYNBICT
Synthetic Biology Curation Tools

![SYNBICT architecture diagram](synbict_architecture_diagram.png)

## Installation instructions

This project depends on Python 3.

To install, run this command after changing directories to SYNBICT:

`python setup.py install`

If you want to visualize circuits, you need to install biopython, matplotlib and the fork of dnaplotlib at https://github.com/nroehner/dnaplotlib.

`pip install matplotlib`


## Testing

`python sequences_to_features.py -n http://synbict.org -t ../test/genetic_nand.xml -f ../test/cello_library.xml`

`python features_to_circuits.py -n http://synbict.org -i test_gate -t ../test/genetic_nand_annotated.xml -c ../test/cello_library.xml`

`python circuit_visualization.py -c ../test/genetic_nand_annotated_circuit.xml -f ../test/genetic_nand_annotated_circuit.xml`

## sequences_to_features.py

sequences_to_features.py annotates an SBOL file.

For input file [NAME].xml, its output is [NAME]\_annotated.xml in the same directory as the input file.

### Arguments for sequences_to_features.py

Argument | Short Arg | Type | Description | Example
---- | --- | --- | --- | ---
`--namespace` | `-n` | `String` | **Required**. Annotation namespace that you own or that you are reasonably certain is only used by you. | http://mynamespace.org
`--target_files` | `-t` | `String` | **Required**. List of paths to SBOL files containing sequences to curate. | mytargets_1.xml mytargets_2.xml
`--feature_files` | `-f` | `String` | **Optional**. List of paths to SBOL files containing sequence features with which to annotate. Default is an empty list. | myfeatures_1.xml myfeatures_2.xml
`--output_files` | `-o` | `String` | **Optional**. List of paths to output annotated SBOL files. By default uses the target file paths postfixed with \"annotated\". | mytargets_1_annotated.xml mytargets_2_annotated.xml
`--version` | `-v` | `String` | **Optional**. Version given to annotated sequences. Default is 1. | 1
`--log_file` | `-l` | `String` | **Optional**. Path to log file to populate with more verbose curation history. Default is to not generate a log file. | mycuration.log
`--min_target_length` | `-m` | `Integer` | **Optional**. Minimum length that a sequence must be to curate. Default is 2000 bp. | 2000
`--min_feature_length` | `-M` | `Integer` | **Optional**. Minimum length that a sequence feature must be to use for annotation. Default is 40 bp. | 40
`--cover_offset` | `-c` | `Integer` | **Optional**. Maximum distance between the start/end of one sequence feature and the start/end of another sequence feature for pruning to initiate. Default is 14 bp. | 14
`--roles` | `-r` | `String` | **Optional**. List of URIs for Sequence Ontology term to filter sequence features and remove any without at least one such term as a role. Default is an empty list. | http://identifiers.org/so/SO:0000167 http://identifiers.org/so/SO:0000316
`--delete_flat_annotations` | `-d` | `Boolean` | **Optional**. If included, flat sequence annotations that do not refer to a sub-part will be automatically deleted. Default is to not automatically delete flat sequence annotations. | -d
`--validate` | `-v` | `Boolean` | **Optional**. If included, output SBOL files will be be validated. Default is to not validate output files. | -v

## features\_to\_circuits.py 

Argument | Short Arg | Type | Description | Example
--- | --- | --- | --- | ---
`--namespace` | `-n` | `String` | **Required**. Annotation namespace that you own or that you are reasonably certain is only used by you. | http://mynamespace.org
`--target_files` | `-t` | `String` | **Required**. List of paths to SBOL files containing DNA features from which to infer genetic circuits. | mytargets_1_annotated.xml mytargets_2_annotated.xml
`--sub_circuit_files` | `-c` | `String` | **Required**. List of paths to SBOL files containing sub-circuits to use in inferring genetic circuits. Default is an empty list. | subcircuits_1.xml subcircuits_2.xml
`--circuit_IDs` | `-i` | `String` | **Optional**. List of IDs given to the inferred genetic circuits (one per target file). By default uses the names of the target files postfixed with \"circuit\". | mytargets_1_annotated_circuit mytargets_2_annotated_circuit
`--output_files` | `-o` | `String` | **Optional**. List of paths to output circuit SBOL files. By default uses the target file paths postfixed with \"circuit\". | mytargets_1_annotated_circuit.xml mytargets_2_annotated_circuit.xml
`--version` | `-v` | `String` | **Optional**. Version given to the inferred genetic circuits. Default is 1. | 1
`--log_file` | `-l` | `String` | **Optional**. Path to log file to populate with more verbose curation history. Default is to not generate a log file. | mycuration.log
`--min_target_length` | `-m` | `Integer` | **Optional**. Minimum length that a DNA feature must be to include it as a genetic device in the inferred circuit. Default is 2000 bp. | 2000
`--validate` | `-x` | `Boolean` | **Optional**. If included, output SBOL files will be be validated. Default is to not validate output files. | -x


-i circuit name (required)
-c parts collection (required)

-n http://foo.bar -i bob 
-t ~/tmp/cpc/Cello_Parts_collection/Strain_3_MG1655_Genomic_IcaR_Gate_annotated.xml 
-c ~/tmp/cpc/Cello_Parts_collection/Cello_Parts_collection.xml 

## circuit_visualization.py

### Arguments
 -c file to visualize, which must have at least one ModuleDefinition
   -h, --help            show this help message and exit
  -c CIRCUIT_FILE, --circuit_file CIRCUIT_FILE
  -f [FEATURE_FILES [FEATURE_FILES ...]], --feature_files [FEATURE_FILES [FEATURE_FILES ...]]
  -l [CURATION_LOG], --curation_log [CURATION_LOG]
  -m [MIN_FEATURES], --min_features [MIN_FEATURES]
  -v, --validate
