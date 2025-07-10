#from memory_profiler import profile
from sequences_to_features import load_sbol, FeatureLibrary, FeatureAnnotater
import os, json, tracemalloc, sys, time
from collections import Counter
from rdflib import URIRef
DCTERMS_TITLE = URIRef('http://purl.org/dc/terms/title')

input_file = sys.argv[1]

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
    return feature_library

#@profile
def test_annotate(feature_library, doc):
    output_docs = []
    output_library = FeatureLibrary(output_docs, False)
    # doc is target_doc
    target_library = FeatureLibrary([doc], False)
    feature_annotater = FeatureAnnotater(feature_library, min_feature_length=40)
    feature_annotater.annotate(target_library, min_target_length=40, in_place=True, output_library=output_library, complete_matches=False, strip_prefixes=[])

def get_title_from_sbol(f2):
    titles = []
    for cd in f2.componentDefinitions:
        for sa in cd.sequenceAnnotations:
            if DCTERMS_TITLE in sa.properties:
                title = sa.properties[DCTERMS_TITLE][0]
                titles.append(str(title))
            else:
                print("  dcterms:title not found")
    return titles

doc = load_sbol(input_file) 
feature_library = get_feature_libraries()
titles1 = get_title_from_sbol(doc)

start = time.perf_counter()
tracemalloc.start()

test_annotate(feature_library, doc)

current, peak = tracemalloc.get_traced_memory()
tracemalloc.stop()
end = time.perf_counter()
elapsed_ms = (end - start) * 1000
titles2 = get_title_from_sbol(doc)
c1 = Counter(titles1)
c2 = Counter(titles2)
# Subtract counts to get new elements in list2
new_elements = list((c2 - c1).elements())

output = {
    "time": elapsed_ms,
    "peak_memory": peak / 1024 / 1024,
    "new_elements": new_elements,
    "number": len(new_elements)
}

out_file = name_only = os.path.splitext(os.path.basename(input_file))[0]+".json"
with open(out_file, 'w') as f:
    json.dump(output, f)
