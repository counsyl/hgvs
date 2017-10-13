from __future__ import unicode_literals

from .. import parse_hgvs_name
from ..utils import read_transcripts
from ..variants import normalize_variant
from .genome import MockGenomeTestFile


def test_name_to_variant_long():
    """
    Convert HGVS names to variant coordinates.

    Test a large number of HGVS names from the wild.
    """
    genome = MockGenomeTestFile(
        db_filename='hg19.fa',
        filename='pyhgvs/tests/data/test_hgvs.genome',
        create_data=False)

    # Read transcripts.
    with open('pyhgvs/data/genes.refGene', 'r') as infile:
        transcripts = read_transcripts(infile)

    class NoTranscriptError(Exception):
        pass

    def get_transcript_long(name):
        """Return a transcript name for the long test."""
        transcript = transcripts.get(name)
        if not transcript:
            raise NoTranscriptError(name)
        chrom = transcript.tx_position.chrom

        # Skip alternative haplotypes.
        if '_' in chrom:
            raise NoTranscriptError(name)

        # Skip sex chromosomes.
        if chrom in ('', 'chrX', 'chrY'):
            raise NoTranscriptError(name)

        return transcript

    errors = []
    with open('pyhgvs/tests/data/test_hgvs.txt', 'r') as infile:
        for i, line in enumerate(infile):
            row = line.rstrip().split('\t')
            chrom, offset, ref, alt, hgvs_name = row[:5]
            offset = int(offset)

            try:
                hgvs_variant = parse_hgvs_name(
                    hgvs_name, genome, get_transcript=get_transcript_long)
            except NoTranscriptError:
                continue

            unnorm_variant = (chrom, offset, ref, alt)
            chrom, offset, ref, alts = normalize_variant(
                chrom, offset, ref, [alt], genome).variant
            variant = (chrom, offset, ref, alts[0])

            if hgvs_variant != variant:
                errors.append(repr([hgvs_variant, variant,
                                    unnorm_variant, hgvs_name]))

    assert not errors, '\n'.join(errors)
