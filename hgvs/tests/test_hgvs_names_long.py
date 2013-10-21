
import nose

from .. import parse_hgvs_name
from ..utils import read_transcripts
from ..variants import normalize_variant
from .genome import MockGenome


def test_name_to_variant_long():
    """
    Convert HGVS names to variant coordinates.

    Test a large number of HGVS names from the wild.
    """
    genome = MockGenome(filename='hgvs/data/test_hgvs.genome')

    #genome = MockGenome(
    #    db_filename='/seq-data/seq/misc/genomics/genomes/hg19/hg19.fa')

    # Read transcripts.
    with open('hgvs/data/genes.refGene') as infile:
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

    with open('hgvs/data/test_hgvs.txt') as infile:
        for i, line in enumerate(infile):
            row = line.rstrip().split('\t')
            chrom, offset, ref, alt, hgvs_name = row[:5]
            offset = int(offset)

            try:
                hgvs_variant = parse_hgvs_name(
                    hgvs_name, genome, get_transcript=get_transcript_long)
            except NoTranscriptError:
                continue
            except ValueError as error:
                yield (nose.tools._ok, False, repr([error, hgvs_name]))
                continue

            try:
                unnorm_variant = (chrom, offset, ref, alt)
                chrom, offset, ref, alts = normalize_variant(
                    chrom, offset, ref, [alt], genome).variant
                variant = (chrom, offset, ref, alts[0])
            except Exception as error:
                yield (nose.tools._ok, False,
                       repr([error, hgvs_variant, hgvs_name]))
                raise

            yield (nose.tools.assert_equal, hgvs_variant, variant,
                   repr([hgvs_variant, variant, unnorm_variant, hgvs_name]))

    #genome.write('hgvs/data/test_hgvs.genome')
