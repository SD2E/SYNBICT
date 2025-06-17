import json
import requests
import sbol2
import logging

# for only 1 sequence
def sbol_sequence(doc):
    seqURI = doc.componentDefinitions[0].sequences[0]
    seq = doc.getSequence(seqURI)
    return seq.elements
    
def load_sbol(sbol_file):

    doc = sbol2.Document()
    doc.read(sbol_file)

    doc.name = sbol_file

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')
    doc.addNamespace('http://sbolstandard.org/gff3#', 'gff3')
    doc.addNamespace('http://cellocad.org/Terms/cello#', 'cello')
    return doc

def load_non_sbol(non_sbol_file):

    conversion_request = {
        'options': {
            'language' : 'SBOL2',
            'test_equality': False,
            'check_uri_compliance': False,
            'check_completeness': False,
            'check_best_practices': False,
            'fail_on_first_error': False,
            'provide_detailed_stack_trace': False,
            'subset_uri': '',
            'uri_prefix': sbol2.getHomespace(),
            'version': '1',
            'insert_type': False,
            'main_file_name': 'main file',
            'diff_file_name': 'comparison file'
        },
        'return_file': True,
        'main_file': open(non_sbol_file).read()
    }

    conversion_response = requests.post("https://validator.sbolstandard.org/validate/", json=conversion_request)

    response_dict = json.loads(conversion_response.content.decode('utf-8'))

    doc = sbol2.Document()
    doc.readString(response_dict['result'])

    doc.name = non_sbol_file

    doc.addNamespace('http://purl.org/dc/elements/1.1/', 'dc')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/igem#', 'igem')
    doc.addNamespace('http://wiki.synbiohub.org/wiki/Terms/synbiohub#', 'sbh')

    return doc

def load_target_file(target_file):
    if target_file.endswith('.xml') or target_file.endswith('.sbol'):
        return load_sbol(target_file)
    elif (target_file.endswith('.gb')
            or target_file.endswith('.genbank')
            or target_file.endswith('.fasta')
            or target_file.endswith('.faa')
            or target_file.endswith('.fa')
            or target_file.endswith('.fas')
            or target_file.endswith('.fsa')):
        return load_non_sbol(target_file)
    else:
        logger.error('Extension of target file %s is unrecognized.', target_file)

        return None
