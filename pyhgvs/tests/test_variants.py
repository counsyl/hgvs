
from unittest import TestCase

from ..variants import normalize_variant
from .genome import MockGenome


_genome_seq = dict([
    (('chr17', 41246230, 41246271), 'AGCCTCATGAGGATCACTGGCCAGTAAGTCTATTTTCTCTG'),  # nopep8
    (('chr17', 41246218, 41246248), 'TTTACATATTAAAGCCTCATGAGGATCACT'),
    (('chr17', 41246249, 41246279), 'GCCAGTAAGTCTATTTTCTCTGAAGAACCA'),
])


_normalize_tests = [
    # Simple SNP.
    (('chr17', 41246250, 'G', ['C']),
     ('chr17', 41246250, 'G', ['C'])),

    # Left-align and 1bp pad.
    (('chr17', 41246251, '', ['G']),
     ('chr17', 41246248, 'T', ['TG'])),

    # Trim common prefix, left-align, and 1bp pad.
    (('chr17', 41246250, 'G', ['GG']),
     ('chr17', 41246248, 'T', ['TG'])),

    # Trim common prefix.
    (('chr17', 41246248, 'TGGC', ['TGGA']),
     ('chr17', 41246251, 'C', ['A'])),

    # Trim common prefix and suffix.
    (('chr17', 41246248, 'TGGC', ['TGAC']),
     ('chr17', 41246250, 'G', ['A'])),

    # Trim common prefix, triallelic
    (('chr17', 41246248, 'TGGC', ['TGGA', 'TGAC']),
     ('chr17', 41246249, 'GGC', ['GAC', 'GGA'])),
]


class TestVariant(TestCase):
    def test_normalize_variant(self):
        """
        Test normalize_variant against known cases.
        """
        genome = MockGenome(_genome_seq)

        for variant, true_variant in _normalize_tests:
            chrom, offset, ref, alts = variant
            norm_variant = normalize_variant(
                chrom, offset, ref, alts, genome).variant
            self.assertEqual(norm_variant, true_variant,
                             'Variant failed to normalize: %s' %
                             repr(variant))
