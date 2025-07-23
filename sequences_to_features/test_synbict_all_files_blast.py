#from memory_profiler import profile
from sequences_to_features import load_sbol, FeatureLibrary
import os, json, tracemalloc, time
from collections import Counter
from rdflib import URIRef
from Annotator import TableFeatureMapper
from FeatureAnnotatorBase import FeatureAnnotatorSimple
from FeatureExtractor import FeatureExtractor
from BlastAligner import BlastAligner
DCTERMS_TITLE = URIRef('http://purl.org/dc/terms/title')

dir = "/home/sophia/git_repo/Addgene-Annotation"
addgene_sbol_list = os.listdir(os.path.join(dir, 'addgene_sbol'))[0:100]

def get_feature_libraries_paths(feature_libraries_dir) -> str:
    # Get the paths of all feature library files in the feature_libraries_dir
    library_files = os.listdir(feature_libraries_dir)
    feature_library_paths = [os.path.join(feature_libraries_dir, f) for f in library_files]
    return feature_library_paths

#@profile
def get_feature_libraries():
    feature_libraries_dir = "/home/sophia/git_repo/SYNBICT/example/jet_libs"
    feature_libraries_paths = get_feature_libraries_paths(feature_libraries_dir)
    feature_docs = []
    for feature_file in feature_libraries_paths:
        feature_docs.append(load_sbol(feature_file))
    feature_library = FeatureLibrary(feature_docs) # load feature_docs to feature_library
    return feature_library, feature_docs

def create_index(feature_docs):
    tmp = FeatureExtractor(feature_docs)
    fasta_path = 'test.fasta' 
    index_prefix = 'test'
    tmp.write_fasta(fasta_path)
    #tmp.write_metadata('test_metadata.json') # '/home/sophia/git_repo/SYNBICT/example/test_metadata.json'
    tmp.build_index(fasta_path, index_prefix, tool='blast') 
    index_prefix = 'test'

    blast = BlastAligner(index_prefix)
    return blast


def test_annotate(feature_library, doc, blast):

    output_sam_path = 'aligned.txt'
    blast.align(doc, output_sam_path, True) # True for exact match

    mapper = TableFeatureMapper('aligned.txt')
    inline_matches, rc_matches = mapper.extract_matches(True) # True for exact match
    simple = FeatureAnnotatorSimple(feature_library, inline_matches, rc_matches)

    #doc = load_sbol("/home/sophia/git_repo/SYNBICT/example/add_gene/100005_addgene_out.xml")
    # doc is target_doc
    target_library = FeatureLibrary([doc], False)
    output_docs = []
    output_library = FeatureLibrary(output_docs, False) 
    simple.annotate(inline_matches, rc_matches, target_library, 40, in_place=True, output_library=output_library, output_matches=False)
    
def get_title_from_sbol(f2):
    titles = []
    for cd in f2.componentDefinitions:
        for sa in cd.sequenceAnnotations:
            if DCTERMS_TITLE in sa.properties:
                if(len(sa.properties[DCTERMS_TITLE]) > 0):
                    title = sa.properties[DCTERMS_TITLE][0]
                    titles.append(str(title))
            else:
                print("dcterms:title not found")
    return titles

feature_library, feature_docs = get_feature_libraries()
blast = create_index(feature_docs)

out_file = "combined_blast_results.json"
with open(out_file, 'w') as f:
    f.write("{\n")  # Start of JSON object
    for ind, file in enumerate(addgene_sbol_list):
        doc = load_sbol(os.path.join(dir, 'addgene_sbol', file))
        titles1 = get_title_from_sbol(doc)

        start = time.perf_counter()
        tracemalloc.start()

        test_annotate(feature_library, doc, blast)

        current, peak = tracemalloc.get_traced_memory()
        peak = peak / 1024 / 1024
        tracemalloc.stop()
        end = time.perf_counter()
        elapsed_ms = (end - start) * 1000
        titles2 = get_title_from_sbol(doc)
        c1 = Counter(titles1)
        c2 = Counter(titles2)
        # Subtract counts to get new elements in list2
        new_elements = list((c2 - c1).elements())

        output = {
            "time": round(elapsed_ms, 2), 
            "peak_memory": round(peak, 2), 
            "new_elements": new_elements,
            "number": len(new_elements)
        }
        key = os.path.basename(file)
        json_str = json.dumps(output, indent=2)

        # Add comma except for the last element
        comma = ',' if ind < len(addgene_sbol_list) - 1 else ''
        f.write(f'  "{key}": {json_str}{comma}\n')
        f.flush()
        os.fsync(f.fileno())

    f.write("}\n")  # End of JSON object
