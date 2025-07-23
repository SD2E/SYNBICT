from memory_profiler import profile
from sequences_to_features import load_sbol, FeatureLibrary
import os
from Annotator import TableFeatureMapper
from FeatureAnnotatorBase import FeatureAnnotatorSimple
from FeatureExtractor import FeatureExtractor
from BlastAligner import BlastAligner

def get_feature_libraries_paths(feature_libraries_dir) -> str:
    # Get the paths of all feature library files in the feature_libraries_dir
    library_files = os.listdir(feature_libraries_dir)
    feature_library_paths = [os.path.join(feature_libraries_dir, f) for f in library_files]
    return feature_library_paths

def get_feature_libraries():
    feature_libraries_dir = "/home/sophia/git_repo/SYNBICT/example/jet_libs"
    feature_libraries_paths = get_feature_libraries_paths(feature_libraries_dir)
    feature_docs = []
    for feature_file in feature_libraries_paths:
        feature_docs.append(load_sbol(feature_file))
    feature_library = FeatureLibrary(feature_docs) # load feature_docs to feature_library
    return feature_library, feature_docs

@profile
def test_annotate(feature_docs, feature_library, doc):
    tmp = FeatureExtractor(feature_docs)
    fasta_path = 'test.fasta' #'/home/sophia/git_repo/SYNBICT/example/test.fasta'
    index_prefix = 'test'
    tmp.write_fasta(fasta_path)
    #tmp.write_metadata('test_metadata.json') # '/home/sophia/git_repo/SYNBICT/example/test_metadata.json'
    tmp.build_index(fasta_path, index_prefix, tool='blast')

    
    index_prefix = 'test'
    blast = BlastAligner(index_prefix)
    output_sam_path = 'aligned.txt'
    blast.align(doc, output_sam_path, exact_match=False)

    mapper = TableFeatureMapper(output_sam_path)
    inline_matches, rc_matches = mapper.extract_matches(False)
    simple = FeatureAnnotatorSimple(feature_library, inline_matches, rc_matches)

    #doc = load_sbol("/home/sophia/git_repo/SYNBICT/example/add_gene/100005_addgene_out.xml")
    # doc is target_doc
    target_library = FeatureLibrary([doc], False)
    output_docs = []
    output_library = FeatureLibrary(output_docs, False) # 
    simple.annotate(inline_matches, rc_matches, target_library, 40, in_place=True, output_library=output_library, output_matches=False)

doc = load_sbol("/home/sophia/git_repo/SYNBICT/example/add_gene/100005_addgene_out.xml") 
feature_library, feature_docs = get_feature_libraries()
test_annotate(feature_docs, feature_library, doc)
print(doc)
