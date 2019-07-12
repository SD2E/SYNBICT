# SYNBICT
Synthetic Biology Curation Tools

![SYNBICT architecture diagram](synbict_architecture_diagram.png)

## sequences_to_features.py

Argument | Short Arg | Type | Description | Example
--- | --- | --- | --- | ---
`--namespace` | `-n` | `String` | **Required**. Annotation namespace that you own or that you are reasonably certain is only used by you. | http://mynamespace.org
`--target_files` | `-t` | `String` | **Required**. List of paths for SBOL files containing sequences to curate. | mytargets_1.xml mytargetfile_2.xml
`--feature_files` | `-f` | `String` | **Optional**. List of paths for SBOL files containing sequence features with which to annotate. Default is an empty list. | myfeatures_1.xml myfeatures2.xml
`--curation_log` | `-l` | `String` | **Optional**. Path for log file to populate with curation history. Default is to not generate a log file. | mycuration.log
`--min_target_length` | `-m` | `Integer` | **Optional**. Minimum length that a sequence must be to curate. Default is 1000 bp. | 1000
`--min_feature_length` | `-M` | `Integer` | **Optional**. Minimum length that a sequence feature must be to use. Default is 40 bp. | 40
`--cover_offset` | `-c` | `Integer` | **Optional**. Maximum distance between the start/end of one sequence feature and the start/end of another sequence feature for pruning to initiate. Default is 14 bp. | 14
`--roles` | `-r` | `String` | **Optional**. List of Sequence Ontology terms to filter annotations and remove any without at least one such term as a role. Default is an empty list. | http://identifiers.org/so/SO:0000167 http://identifiers.org/so/SO:0000316
`--validate` | `-v` | `Boolean` | **Optional**. If included, output SBOL file will be validate. Default is no validation. | -v
