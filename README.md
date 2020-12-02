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

sequences_to_features.py annotates and prunes sequences in target SBOL, GenBank, or FASTA files.

### Common arguments for sequences_to_features.py

Argument | Short Arg | Type | Description | Example
---- | --- | --- | --- | ---
`--namespace` | `-n` | `String` | **Required**. Namespace that you own or that you are reasonably certain is only used by you. | http://mynamespace.org
`--target_files` | `-t` | `String` | **Required**. List of input files and/or directories containing sequences to curate. Accepted file formats include SBOL XML, FASTA, and GenBank. | targets1.xml targets2.fa
`--output_files` | `-o` | `String` | **Optional**. List of output files. If length is less than that of target files and if a curation suffix is provided, then the difference is populated with the corresponding target files postfixed with the curation suffix. Alternatively, if list contains a single directory,  | targets1_annotated.xml targets2_annotated.xml
`--curation_suffix` | `-s` | `String` | **Optional**. Suffix for postfixing output files. | annotated
`--log_file` | `-l` | `String` | **Optional**. Log file to populate with more verbose curation history. Default is to not generate a log file. | curation.log
`--validate` | `-v` | `Boolean` | **Optional**. If included, output files will be checked against SBOL validation rules. Default is to not validate output files. | -v
`--in_place` | `-p` | `Boolean` | **Optional**. If included, do not copy target sequences prior to curation. Default is to curate copies of target sequences. | 

### Annotation arguments for sequences_to_features.py

Argument | Short Arg | Type | Description | Example
---- | --- | --- | --- | ---
`--feature_files` | `-f` | `String` | **Optional**. List of files and/or directories containing features for annotating sequences. Default is an empty list. Accepted file format is SBOL XML. | features1.xml feature_directory
`--min_target_length` | `-m` | `Integer` | **Optional**. Minimum length that a sequence must be to curate. Default is 2000 bp. | 2000
`--min_feature_length` | `-M` | `Integer` | **Optional**. Minimum length that feature must be for annotating sequences. Default is 40 bp. | 40
`--no_annotation` | `-na` | `Boolean` | **Optional**. If included, do not annotate target sequences. Default is to annotate target sequences. | 
`--extend_features` | `-e` | `Boolean` | **Optional**. If included, attempt to extend feature library. Derive new features from annotations in target sequences only if their names align to names of library features and the fraction mismatch between their sequences is less than the extension threshold. Default is to not extend features. | 
`--extension_threshold` | `-x` | `Float` | **Optional**. Fraction mismatch permitted between sequences when extending feature library. | 0.05
`--extension_suffix` | `-xs` | `String` | **Optional**. Suffix for postfixing extended feature file. | extended

### Pruning arguments for sequences_to_features.py

Argument | Short Arg | Type | Description | Example
---- | --- | --- | --- | ---
`--cover_offset` | `-c` | `Integer` | **Optional**. Maximum distance between the start/end of one feature and the start/end of another  feature for pruning of overlapping features to initiate. Default is 14 bp. | 14
`--deletion_roles` | `-r` | `String` | **Optional**. List of URIs for Sequence Ontology roles to remove all matching features from target sequence. Default is an empty list. | http://identifiers.org/so/SO:0000167 http://identifiers.org/so/SO:0000316
`--delete_flat` | `-d` | `Boolean` | **Optional**. If included, automatically delete annotations that do not refer to a sub-part. Default is no automatic deletion. | 
`--auto_swap` | `-a` | `Boolean` | **Optional**. If included, automatically merge any overlapping flat annotation and sub-part annotation that do not overlap with any other features. Default is to ask user if merger should take place. | 
`--non_interactive` | `-ni` | `Boolean` | **Optional**. If included, do not ask user for additional input. Default is to ask user. | 
`--no_pruning` | `-np` | `Boolean` | **Optional**. If included, do not prune target sequences. Default is to prune target sequences. | 

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
