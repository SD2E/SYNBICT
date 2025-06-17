
class Aligner():
    """
    Abstract base class for sequence alignment tools.
    Subclasses must implement the `align` method.
    """

    def __init__(self, index_prefix):
        """
        Initialize the aligner with a given index prefix.

        Args:
            index_prefix (str): The prefix used for indexing the reference sequence.
        """
        self.index_prefix = index_prefix

    def align(self, query_sbol, output_sam_path):
        """
        Align the input SBOL query sequence and produce a SAM output.

        Args:
            query_sbol: SBOL Document or ComponentDefinition to align.
            output_sam_path (str): Path to output SAM file.
        """
        pass
