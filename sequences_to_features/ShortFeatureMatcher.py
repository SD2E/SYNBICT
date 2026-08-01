"""Exhaustive exact-substring matcher for short features.

BLASTN (and seed-based aligners in general) cannot report matches shorter than
~14 bp against a long (plasmid-length) query: the significance threshold
suppresses them regardless of word size. Features in [min_length, max_length]
(default 9-13 bp) are therefore handled here by direct forward + reverse-
complement substring search, and the results merge into the same annotation
pass as the aligner's matches for the longer features.

The output format -- (feature_list, start, end) tuples split into
inline_matches / rc_matches with 0-based, half-open [start, end) query
coordinates -- is identical to TableFeatureMapper / SAMFeatureMapper, so the
two match sets can simply be concatenated before FeatureAnnotatorSimple.annotate().
"""
from .Feature import Feature

_COMP = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N",
         "a": "t", "t": "a", "g": "c", "c": "g"}


def merge_matches(aligner_inline, aligner_rc, short_inline, short_rc):
    """Concatenate an aligner's matches (for long features) with a
    ShortFeatureMatcher's matches (for short features) so a single
    FeatureAnnotatorSimple.annotate() pass covers both length regimes.

    Returns (inline_matches, rc_matches). The two match sets use the same
    (feature_list, start, end) tuple format, so this is a plain concatenation.
    """
    return list(aligner_inline) + list(short_inline), list(aligner_rc) + list(short_rc)


def _reverse_complement(seq):
    return "".join(_COMP.get(c, "N") for c in reversed(seq))


def _find_all(haystack, needle):
    positions = []
    i = haystack.find(needle)
    while i != -1:
        positions.append(i)
        i = haystack.find(needle, i + 1)
    return positions


class ShortFeatureMatcher:
    """Index the short features of a library once, then match a target by exact
    substring search.

    Parameters
    ----------
    feature_library : FeatureLibrary
    min_length, max_length : int
        Inclusive length window handled here (default 9-13 bp). Set max_length to
        one below the aligner's minimum (e.g. aligner handles >=14 bp) so the two
        passes partition the parts without overlap.
    """

    def __init__(self, feature_library, min_length=9, max_length=13):
        self.feature_library = feature_library
        self.min_length = min_length
        self.max_length = max_length
        # Group identical sequences so a single search annotates every synonymous
        # feature ID at that locus (mirrors FlashText keyword accumulation).
        self._by_sequence = {}
        for feature in feature_library.features:
            seq = str(feature.nucleotides).upper()
            if min_length <= len(seq) <= max_length:
                self._by_sequence.setdefault(seq, []).append(feature.identity)
        self.inline_matches = []
        self.rc_matches = []

    def get_short_feature_count(self):
        return sum(len(ids) for ids in self._by_sequence.values())

    def extract_matches(self, target_nucleotides, min_feature_length=None):
        """Return (inline_matches, rc_matches) for the target sequence.

        min_feature_length, if given, raises the lower length bound (kept so this
        matcher is call-compatible with the aligner mappers). Every exact
        forward occurrence of a short feature is an inline match; every exact
        reverse-complement occurrence is an rc match.
        """
        lower = self.min_length if min_feature_length is None \
            else max(self.min_length, min_feature_length)
        self.inline_matches = []
        self.rc_matches = []
        sequence = str(target_nucleotides).upper()

        for feat_seq, identities in self._by_sequence.items():
            length = len(feat_seq)
            if length < lower or length > self.max_length:
                continue
            features = [Feature(nucleotides='', identity=i, roles='') for i in identities]

            for pos in _find_all(sequence, feat_seq):
                self.inline_matches.append((list(features), pos, pos + length))

            rc_seq = _reverse_complement(feat_seq)
            if rc_seq != feat_seq:  # a palindrome is already covered by the forward hit
                for pos in _find_all(sequence, rc_seq):
                    self.rc_matches.append((list(features), pos, pos + length))

        return self.inline_matches, self.rc_matches
