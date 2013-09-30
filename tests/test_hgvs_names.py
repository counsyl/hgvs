
from StringIO import StringIO

import nose

from ..import CDNACoord
from ..import CDNA_STOP_CODON
from ..import HGVSName
from ..import format_hgvs_name
from ..import parse_hgvs_name
from ..import genomic_to_cdna_coord
from ..import cdna_to_genomic_coord
from ..utils import read_refgene
from ..utils import make_transcript
from .genome import MockGenome


def get_transcript(name):
    """
    Fetch a transcript from its name or gene name.
    """
    return _transcripts.get(name)


def test_parse_cdna_coord():
    """
    Parse cDNA coordinates.
    """
    for text, expected in _parse_cdna_coords:
        nose.tools.assert_equal(CDNACoord(string=text), expected)


def test_fromat_cdna_coord():
    """
    Format cDNA coordinates.
    """
    for expected_text, coord in _parse_cdna_coords:
        nose.tools.assert_equal(str(coord), expected_text)


def test_genomic_to_cdna_coord():
    """
    Convert genomic to cDNA coordinates.
    """
    for transcript_name, genomic_coord, cdna_coord_expected in _convert_coords:
        transcript = get_transcript(transcript_name)
        cdna_coord = genomic_to_cdna_coord(transcript, genomic_coord[1])
        nose.tools.assert_equal(
            cdna_coord, cdna_coord_expected,
            repr((cdna_coord, cdna_coord_expected,
                  transcript_name, genomic_coord)))


def test_cdna_to_genomic_coord():
    """
    Convert cDNA to genomic coordinates.
    """
    for transcript_name, genomic_coord_expected, cdna_coord in _convert_coords:
        transcript = get_transcript(transcript_name)
        genomic_coord = cdna_to_genomic_coord(transcript, cdna_coord)
        nose.tools.assert_equal(
            genomic_coord, genomic_coord_expected[1],
            repr((genomic_coord, genomic_coord_expected[1],
                  transcript_name, cdna_coord)))


def test_parse_name():
    """
    Parsing HGVS names.
    """
    for name, formatable, expected in _parse_names:
        hgvs_parsed = HGVSName(name)
        for key, value in expected.items():
            nose.tools.assert_equal(
                getattr(hgvs_parsed, key), value,
                (getattr(hgvs_parsed, key), value, name, expected))


def test_format_name():
    """
    Format HGVS names.
    """
    for expected_name, formatable, attrs in _parse_names:
        if formatable:
            name = HGVSName(**attrs).format()
            nose.tools.assert_equal(name, expected_name,
                                    (name, expected_name, attrs))


def test_name_to_variant():
    """
    Convert HGVS names to variant coordinates.
    """
    genome = MockGenome(_genome_seq)

    for hgvs_name, variant, canonical in _name_variants:
        hgvs_variant = parse_hgvs_name(hgvs_name, genome,
                                       get_transcript=get_transcript)
        nose.tools.assert_equal(
            hgvs_variant, variant,
            repr([hgvs_name, variant, hgvs_variant]))


def test_variant_to_name():
    """
    Convert variant coordinates to HGVS names.
    """
    genome = MockGenome(_genome_seq)

    for expected_hgvs_name, variant, canonical in _name_variants:
        if canonical:
            transcript_name = HGVSName(expected_hgvs_name).transcript
            transcript = get_transcript(transcript_name)
            assert transcript, transcript_name
            chrom, offset, ref, alt = variant
            hgvs_name = format_hgvs_name(
                chrom, offset, ref, alt, genome, transcript,
                use_gene=False)
            nose.tools.assert_equal(
                hgvs_name, expected_hgvs_name,
                repr([hgvs_name, expected_hgvs_name, variant]))


# Test examples of cDNA coordinates.
_parse_cdna_coords = [
    ('1001', CDNACoord(1001)),
    ('-1001', CDNACoord(-1001)),
    ('*1001', CDNACoord(1001, landmark=CDNA_STOP_CODON)),
    ('1001+5', CDNACoord(1001, 5)),
    ('1001-5', CDNACoord(1001, -5)),
    ('-1001+5', CDNACoord(-1001, 5)),
    ('-1001-5', CDNACoord(-1001, -5)),
    ('*1001+5', CDNACoord(1001, 5, CDNA_STOP_CODON)),
    ('*1001-5', CDNACoord(1001, -5, CDNA_STOP_CODON)),
]


# Test examples of coverting coordinates.
_convert_coords = [
    # Positions near start codon.
    ('NM_000016.4', ('chr1', 76190473), CDNACoord(1)),
    ('NM_000016.4', ('chr1', 76190472), CDNACoord(-1)),
    ('NM_000016.4', ('chr1', 76190043), CDNACoord(-430)),
    ('NM_007294.3', ('chr17', 41276112), CDNACoord(2)),
    ('NM_007294.3', ('chr17', 41276113), CDNACoord(1)),
    ('NM_007294.3', ('chr17', 41276114), CDNACoord(-1)),

    # Positions near introns.
    ('NM_000016.4', ('chr1', 76190502), CDNACoord(30)),
    ('NM_000016.4', ('chr1', 76190503), CDNACoord(30, 1)),
    ('NM_000016.4', ('chr1', 76194085), CDNACoord(31, -1)),
    ('NM_000016.4', ('chr1', 76194086), CDNACoord(31)),
    ('NM_007294.3', ('chr17', 41276034), CDNACoord(80)),
    ('NM_007294.3', ('chr17', 41276033), CDNACoord(80, 1)),
    ('NM_007294.3', ('chr17', 41267797), CDNACoord(81, -1)),
    ('NM_007294.3', ('chr17', 41267796), CDNACoord(81)),

    # Positions near stop codon.
    ('NM_000016.4', ('chr1', 76228448), CDNACoord(1266)),
    ('NM_000016.4', ('chr1', 76228449), CDNACoord(1, 0, CDNA_STOP_CODON)),
    ('NM_000016.4', ('chr1', 76228450), CDNACoord(2, 0, CDNA_STOP_CODON)),
    ('NM_007294.3', ('chr17', 41197695), CDNACoord(5592)),
    ('NM_007294.3', ('chr17', 41197694), CDNACoord(1, 0, CDNA_STOP_CODON)),
    ('NM_007294.3', ('chr17', 41197693), CDNACoord(2, 0, CDNA_STOP_CODON)),

    # Positions near UTR introns.
    ('NM_007294.3', ('chr17', 41276142), CDNACoord(-19, -10)),
    ('NM_000038.5', ('chr5', 112090570), CDNACoord(-18)),
    ('NM_000038.5', ('chr5', 112090569), CDNACoord(-18, -1)),
    ('NM_000038.5', ('chr5', 112073622), CDNACoord(-19)),
    ('NM_000023.2', ('chr17', 48252799), CDNACoord(1, 0, CDNA_STOP_CODON)),
    ('NM_000023.2', ('chr17', 48252800), CDNACoord(2, 0, CDNA_STOP_CODON)),
    ('NM_000023.2', ('chr17', 48252810), CDNACoord(12, 0, CDNA_STOP_CODON)),
    ('NM_000023.2', ('chr17', 48252811), CDNACoord(12, 1, CDNA_STOP_CODON)),
    ('NM_000023.2', ('chr17', 48253073), CDNACoord(13, 0, CDNA_STOP_CODON)),
    ('NM_000023.2', ('chr17', 48253072), CDNACoord(13, -1, CDNA_STOP_CODON)),

    # Positions flanking the transcript.
    ('NM_007294.3', ('chr17', 41196313), CDNACoord(1382, 0, CDNA_STOP_CODON)),
    ('NM_007294.3', ('chr17', 41196312), CDNACoord(1383, 0, CDNA_STOP_CODON)),
    ('NM_007294.3', ('chr17', 41196311), CDNACoord(1384, 0, CDNA_STOP_CODON)),
    ('NM_007294.3', ('chr17', 41277500), CDNACoord(-232)),
    ('NM_007294.3', ('chr17', 41277501), CDNACoord(-233)),
    ('NM_000016.4', ('chr1', 76190042), CDNACoord(-431)),
    ('NM_000016.4', ('chr1', 76190043), CDNACoord(-430)),
    ('NM_000016.4', ('chr1', 76229354), CDNACoord(906, 0, CDNA_STOP_CODON)),
    ('NM_000016.4', ('chr1', 76229355), CDNACoord(907, 0, CDNA_STOP_CODON)),
    ('NM_000016.4', ('chr1', 76229356), CDNACoord(908, 0, CDNA_STOP_CODON)),
]


# Test examples of HGVS names.
# (name, formatable, hgvs_attrs)
_parse_names = [
    # Names with prefixes.
    ('NM_007294.3:c.2207A>C', True,
     {
         'transcript': 'NM_007294.3',
         'gene': '',
         'kind': 'c',
         'cdna_start': CDNACoord(2207),
         'cdna_end': CDNACoord(2207),
         'ref_allele': 'A',
         'alt_allele': 'C',
         'mutation_type': '>',
     }),
    ('NM_007294.3(BRCA1):c.2207A>C', False,
     {
         'transcript': 'NM_007294.3',
         'gene': 'BRCA1',
     }),
    ('BRCA1{NM_007294.3}:c.2207A>C', True,
     {
         'transcript': 'NM_007294.3',
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(2207),
         'cdna_end': CDNACoord(2207),
         'ref_allele': 'A',
         'alt_allele': 'C',
         'mutation_type': '>',

     }),
    ('BRCA1:c.2207A>C', True,
     {
         'transcript': '',
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(2207),
         'cdna_end': CDNACoord(2207),
         'ref_allele': 'A',
         'alt_allele': 'C',
         'mutation_type': '>',
     }),
    ('NC_000001:g.1000000A>C', False,
     {
         'chrom': 'NC_000001',
     }),
    ('chr7:g.1000000A>C', True,
     {
         'chrom': 'chr7',
         'kind': 'g',
         'start': 1000000,
         'end': 1000000,
         'ref_allele': 'A',
         'alt_allele': 'C',
         'mutation_type': '>',
     }),
    ('BRCA1{NM_007294.3}:g.2207A>C', True,
     {
         'transcript': 'NM_007294.3',
         'gene': 'BRCA1',
         'chrom': '',
         'kind': 'g',
         'start': 2207,
         'end': 2207,
         'ref_allele': 'A',
         'alt_allele': 'C',
         'mutation_type': '>',
     }),

    # cDNA 1bp mutations.
    ('BRCA1:c.101A>C', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(101),
         'cdna_end': CDNACoord(101),
         'ref_allele': 'A',
         'alt_allele': 'C',
         'mutation_type': '>',
     }),
    ('BRCA1:c.101insA', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(101),
         'cdna_end': CDNACoord(101),
         'ref_allele': '',
         'alt_allele': 'A',
         'mutation_type': 'ins',
     }),
    ('BRCA1:c.101delA', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(101),
         'cdna_end': CDNACoord(101),
         'ref_allele': 'A',
         'alt_allele': '',
         'mutation_type': 'del',
     }),
    ('BRCA1:c.101dupA', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(101),
         'cdna_end': CDNACoord(101),
         'ref_allele': 'A',
         'alt_allele': 'AA',
         'mutation_type': 'dup',
     }),

    # cDNA range mutations.
    ('BRCA1:c.1000_1001insATG', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(1000),
         'cdna_end': CDNACoord(1001),
         'ref_allele': '',
         'alt_allele': 'ATG',
         'mutation_type': 'ins',
     }),
    ('BRCA1:c.1000_1002delATG', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(1000),
         'cdna_end': CDNACoord(1002),
         'ref_allele': 'ATG',
         'alt_allele': '',
         'mutation_type': 'del',
     }),
    ('BRCA1:c.1000_1002del', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(1000),
         'cdna_end': CDNACoord(1002),
         'ref_allele': '',
         'alt_allele': '',
         'mutation_type': 'del',
     }),
    ('BRCA1:c.1000_1002dupATG', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(1000),
         'cdna_end': CDNACoord(1002),
         'ref_allele': 'ATG',
         'alt_allele': 'ATGATG',
         'mutation_type': 'dup',
     }),
    ('BRCA1:c.1000+5_1000+6insATG', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(1000, 5),
         'cdna_end': CDNACoord(1000, 6),
         'ref_allele': '',
         'alt_allele': 'ATG',
         'mutation_type': 'ins',
     }),

    # cDNA Indels.
    ('BRCA1:c.3428delCinsTA', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(3428),
         'cdna_end': CDNACoord(3428),
         'ref_allele': 'C',
         'alt_allele': 'TA',
         'mutation_type': 'delins',
     }),
    ('BRCA1:c.3428_3429delCAinsTA', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(3428),
         'cdna_end': CDNACoord(3429),
         'ref_allele': 'CA',
         'alt_allele': 'TA',
         'mutation_type': 'delins',
     }),
    ('BRCA1:c.3428_3429delinsTA', True,
     {
         'gene': 'BRCA1',
         'kind': 'c',
         'cdna_start': CDNACoord(3428),
         'cdna_end': CDNACoord(3429),
         'ref_allele': '',
         'alt_allele': 'TA',
         'mutation_type': 'delins',
     }),

    # Protein 1bp change.
    ('p.Glu1000=', True,
     {
         'kind': 'p',
         'start': 1000,
         'end': 1000,
         'ref_allele': 'Glu',
         'ref2_allele': 'Glu',
         'alt_allele': 'Glu',
         'pep_extra': '=',
         'mutation_type': '>',
     }),
    ('p.Glu1000Ser', True,
     {
         'kind': 'p',
         'start': 1000,
         'end': 1000,
         'ref_allele': 'Glu',
         'ref2_allele': 'Glu',
         'alt_allele': 'Ser',
         'pep_extra': '',
         'mutation_type': '>',
     }),
    ('p.Glu1000_Ser1003?fs', True,
     {
         'kind': 'p',
         'start': 1000,
         'end': 1003,
         'ref_allele': 'Glu',
         'ref2_allele': 'Ser',
         'alt_allele': '',
         'pep_extra': '?fs',
         'mutation_type': 'delins',
     }),
    ('p.Glu1000_Ser1003Aln?fs', False,
     {
         'kind': 'p',
         'start': 1000,
         'end': 1003,
         'ref_allele': 'Glu',
         'ref2_allele': 'Ser',
         'alt_allele': 'Aln',
         'pep_extra': '?fs',
         'mutation_type': 'delins',
     }),

    # Genomic 1bp mutations.
    ('BRCA1:g.101A>C', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 101,
         'end': 101,
         'ref_allele': 'A',
         'alt_allele': 'C',
         'mutation_type': '>',
     }),
    ('BRCA1:g.101insA', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 101,
         'end': 101,
         'ref_allele': '',
         'alt_allele': 'A',
         'mutation_type': 'ins',
     }),
    ('BRCA1:g.101delA', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 101,
         'end': 101,
         'ref_allele': 'A',
         'alt_allele': '',
         'mutation_type': 'del',
     }),
    ('BRCA1:g.101dupA', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 101,
         'end': 101,
         'ref_allele': 'A',
         'alt_allele': 'AA',
         'mutation_type': 'dup',
     }),

    # Genomic range mutations.
    ('BRCA1:g.1000_1001insATG', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 1000,
         'end': 1001,
         'ref_allele': '',
         'alt_allele': 'ATG',
         'mutation_type': 'ins',
     }),
    ('BRCA1:g.1000_1002delATG', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 1000,
         'end': 1002,
         'ref_allele': 'ATG',
         'alt_allele': '',
         'mutation_type': 'del',
     }),
    ('BRCA1:g.1000_1002del', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 1000,
         'end': 1002,
         'ref_allele': '',
         'alt_allele': '',
         'mutation_type': 'del',
     }),
    ('BRCA1:g.1000_1002dupATG', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 1000,
         'end': 1002,
         'ref_allele': 'ATG',
         'alt_allele': 'ATGATG',
         'mutation_type': 'dup',
     }),

    # Genomic indels.
    ('BRCA1:g.3428delCinsTA', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 3428,
         'end': 3428,
         'ref_allele': 'C',
         'alt_allele': 'TA',
         'mutation_type': 'delins',
     }),
    ('BRCA1:g.3428_3429delCAinsTA', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 3428,
         'end': 3429,
         'ref_allele': 'CA',
         'alt_allele': 'TA',
         'mutation_type': 'delins',
     }),
    ('BRCA1:g.3428_3429delinsTA', True,
     {
         'gene': 'BRCA1',
         'kind': 'g',
         'start': 3428,
         'end': 3429,
         'ref_allele': '',
         'alt_allele': 'TA',
         'mutation_type': 'delins',
     }),
]


# Example HGVS names and variants.
# format: (name, variant, is_canonical)
_name_variants = [
    # Simple SNPs.
    ('NM_000352.3:c.215A>G', ('chr11', 17496508, 'T', 'C'), True),
    ('NM_000352.3:c.72C>A', ('chr11', 17498252, 'G', 'T'),  True),
    ('NM_000352.3:c.3885C>G', ('chr11', 17418843, 'G', 'C'), True),

    # SNPs within introns.
    ('NM_000352.3:c.1630+1G>A', ('chr11', 17464266, 'C', 'T'), True),
    ('NM_000352.3:c.1630+1G>C', ('chr11', 17464266, 'C', 'G'), True),
    ('NM_000352.3:c.1630+1G>T', ('chr11', 17464266, 'C', 'A'), True),
    ('NM_000352.3:c.1672-20A>G', ('chr11', 17452526, 'T', 'C'), True),
    ('NM_000352.3:c.1672-74G>A', ('chr11', 17452580, 'C', 'T'), True),
    ('NM_000352.3:c.1923+5G>T', ('chr11', 17450107, 'C', 'A'), True),
    ('NM_000352.3:c.2041-21G>A', ('chr11', 17449510, 'C', 'T'), True),
    ('NM_000352.3:c.2116+3A>G', ('chr11', 17449411, 'T', 'C'), True),
    ('NM_000352.3:c.2116+1G>T', ('chr11', 17449413, 'C', 'A'), True),
    ('NM_000352.3:c.2116+2T>C', ('chr11', 17449412, 'A', 'G'), True),

    # Indels.
    ('NM_000018.3:c.922_930delGCAGAGGTGinsTCAAAGCAC',
     ('chr17', 7126028, 'AGCAGAGGTG', 'ATCAAAGCAC'), False),
    ('NM_000018.3:c.1077_1077+1delGGinsCAC',
     ('chr17', 7126183, 'CGG', 'CCAC'), True),
    ('NM_000019.3:c.163_167delTTTTTinsAA',
     ('chr11', 108004588, 'TTTTTT', 'TAA'), False),
    ('NM_000022.2:c.781-3_781delCAGAinsTGGAAGAGCAGATCTGG',
     ('chr20', 43251292, 'ATCTG', 'ACCAGATCTGCTCTTCCA'), False),
    ('NM_000023.2:c.585-2_585-1delAGinsT',
     ('chr17', 48246450, 'CAG', 'CT'), True),
    ('NM_000030.2:c.2_3delTGinsAT',
     ('chr2', 241808283, 'ATG', 'AAT'), True),
    ('NM_000030.2:c.2_3del2ins2',
     ('chr2', 241808283, 'ATG', 'ANN'), False),
    ('NM_007294.3:c.4185+2_4185+22del21insA',
     ('chr17', 41242938, 'GCACACACACACACGCTTTTTA', 'GT'), False),

    # Single letter del and insert.
    ('NM_000016.4:c.945+4delAinsGC',
     ('chr1', 76216234, u'AA', 'AGC'), True),
    ('NM_000030.2:c.308delGinsTCCTGGTTGA',
     ('chr2', 241808728, u'GG', 'GTCCTGGTTGA'), False),
    ('NM_000038.5:c.1617delCinsGAA',
     ('chr5', 112163693, u'AC', 'AGAA'), True),
    ('NM_000038.5:c.4256delGinsCC',
     ('chr5', 112175546, u'AG', 'ACC'), True),
    ('NM_000038.5:c.4256delGins5',
     ('chr5', 112175546, u'AG', 'ANNNNN'), False),

    # Delete region.
    ('NM_000016.4:c.291_296delTCTTGG',
     ('chr1', 76199214, u'AGGTCTT', 'A'),  False),
    ('NM_000016.4:c.306_307insG',
     ('chr1', 76199232, u'T', 'TG'), False),
    ('NM_000016.4:c.343_348delGGATGT',
     ('chr1', 76199267, u'ATGGATG', 'A'), False),
    ('NM_000016.4:c.430_432delAAG', ('chr1', 76200511, u'AAAG', 'A'), True),
    ('NM_000016.4:c.430_432del3', ('chr1', 76200511, u'AAAG', 'A'), False),

    # Single letter insert, delete, duplication.
    ('NM_000016.4:c.203delA', ('chr1', 76198412, 'GA', 'G'), True),
    ('NM_000016.4:c.244dupT', ('chr1', 76198564, 'C', 'CT'), True),
    ('NM_000016.4:c.387+1delG', ('chr1', 76199309, 'TG', 'T'), True),
    ('NM_000016.4:c.475delT', ('chr1', 76205669, 'AT', 'A'), True),
    ('NM_000016.4:c.1189dupT', ('chr1', 76227049, 'C', 'CT'), True),
    ('NM_000016.4:c.1191delT', ('chr1', 76227051, 'AT', 'A'), True),
    ('NM_000016.4:c.307insG', ('chr1', 76199232, 'T', 'TG'), True),

    # Alignment tests for HGVS 3' and VCF left-alignment.
    ('NM_000492.3:c.935_937delTCT', ('chr7', 117180210, 'CCTT', 'C'), True),
    ('NM_000492.3:c.442delA', ('chr7', 117171120, 'CA', 'C'), True),
    ('NM_000492.3:c.805_806delAT', ('chr7', 117176660, 'AAT', 'A'), True),
    ('NM_000492.3:c.1155_1156dupTA', ('chr7', 117182104, 'A', 'AAT'), True),
    ('NM_000492.3:c.3889dupT', ('chr7', 117292905, 'A', 'AT'), True),

    # Transcript prefix.
    ('NM_007294.3:c.2207A>C', ('chr17', 41245341, 'T', 'G'), False),
    ('NM_007294.3:c.2207A>C', ('chr17', 41245341, 'T', 'G'), False),
    ('NM_007294.3(BRCA1):c.2207A>C', ('chr17', 41245341, 'T', 'G'), False),
    ('BRCA1{NM_007294.3}:c.2207A>C', ('chr17', 41245341, 'T', 'G'), False),

    # After stop codon.
    ('NM_000492.3:c.*3A>C', ('chr7', 117307165, 'A', 'C'), True),

    # Genomic simple SNPs.
    ('chr11:g.17496508T>C', ('chr11', 17496508, 'T', 'C'), False),

    # Genomic indels.
    ('chr17:g.7126029_7126037delGCAGAGGTGinsTCAAAGCAC',
     ('chr17', 7126028, 'AGCAGAGGTG', 'ATCAAAGCAC'), False),

    # Genomic single letter del and insert.
    ('chr1:g.76216235delAinsGC', ('chr1', 76216234, u'AA', 'AGC'), False),

    # Genomic delete region.
    ('chr1:g.76199215_76199220delGGTCTT',
     ('chr1', 76199214, u'AGGTCTT', 'A'), False),
]


# Mock refGene transcripts.
_refgene = '\n'.join([
    '1166	NM_000016.4	chr1	+	76190042	76229355	76190472	76228448	12	76190042,76194085,76198328,76198537,76199212,76200475,76205664,76211490,76215103,76216135,76226806,76228376,	76190502,76194173,76198426,76198607,76199313,76200556,76205795,76211599,76215244,76216231,76227055,76229355,	0	ACADM	cmpl	cmpl	0,0,1,0,1,0,0,2,0,0,0,0,',  # nopep8
    '89	NM_000352.3	chr11	-	17414431	17498449	17414537	17498323	39	17414431,17415243,17415812,17416718,17417156,17417398,17418462,17418739,17419230,17419885,17424207,17426058,17427040,17428168,17428434,17428900,17429938,17432062,17434212,17434940,17436058,17436850,17438476,17448595,17449413,17449835,17450111,17452360,17453750,17464266,17464724,17470062,17474665,17482034,17483129,17484984,17491647,17496432,17498175,	17414675,17415306,17415946,17416822,17417265,17417477,17418593,17418860,17419344,17419988,17424300,17426216,17427110,17428335,17428676,17429000,17430064,17432200,17434293,17435025,17436157,17436886,17438509,17448701,17449489,17449952,17450217,17452506,17453791,17464429,17464859,17470218,17474830,17482223,17483372,17485151,17491769,17496574,17498449,	0	ABCC8	cmpl	cmpl	0,0,1,2,1,0,1,0,0,2,2,0,2,0,1,0,0,0,0,2,2,2,2,1,0,0,2,0,1,0,0,0,0,0,0,1,2,1,0,',  # nopep8
    '21	NM_000019.3	chr11	+	107992257	108018891	107992333	108018117	12	107992257,108002633,108004546,108004947,108005868,108009624,108010791,108012331,108013163,108014709,108016928,108017996,	107992405,108002681,108004664,108005043,108005969,108009768,108010942,108012427,108013277,108014774,108017086,108018891,	0	ACAT1	cmpl	cmpl	0,0,0,1,1,0,0,1,1,1,0,2,',  # nopep8
    '639	NM_000018.3	chr17	+	7123149	7128586	7123303	7128416	20	7123149,7123440,7123782,7123922,7124084,7124242,7124856,7125270,7125495,7125985,7126451,7126962,7127131,7127286,7127464,7127639,7127798,7127960,7128127,7128275,	7123365,7123516,7123848,7123995,7124149,7124377,7125001,7125400,7125621,7126184,7126556,7127049,7127194,7127388,7127562,7127712,7127871,7128033,7128203,7128586,	0	ACADVL	cmpl	cmpl	0,2,0,0,1,0,0,1,2,2,0,0,0,0,0,2,0,1,2,0,',  # nopep8
    '899	NM_007294.3	chr17	-	41196311	41277500	41197694	41276113	23	41196311,41199659,41201137,41203079,41209068,41215349,41215890,41219624,41222944,41226347,41228504,41234420,41242960,41243451,41247862,41249260,41251791,41256138,41256884,41258472,41267742,41276033,41277287,	41197819,41199720,41201211,41203134,41209152,41215390,41215968,41219712,41223255,41226538,41228631,41234592,41243049,41246877,41247939,41249306,41251897,41256278,41256973,41258550,41267796,41276132,41277500,	0	BRCA1	cmpl	cmpl	1,0,1,0,0,1,1,0,1,2,1,0,1,1,2,1,0,1,2,2,2,0,-1,',  # nopep8
    '953	NM_000023.2	chr17	+	48243365	48253293	48243401	48252798	10	48243365,48244728,48244942,48245307,48245734,48246452,48247503,48248000,48252617,48253072,	48243438,48244848,48245097,48245380,48245933,48246615,48247712,48248027,48252810,48253293,	0	SGCA	cmpl	cmpl	0,1,1,0,1,2,0,2,2,-1,',  # nopep8
    '2429	NM_000030.2	chr2	+	241808161	241818536	241808282	241818238	11	241808161,241808586,241810060,241810765,241812395,241813394,241814525,241815351,241816953,241817438,241818130,	241808447,241808779,241810125,241810866,241812466,241813479,241814621,241815421,241817049,241817567,241818536,	0	AGXT	cmpl	cmpl	0,0,1,0,2,1,2,2,0,0,0,',  # nopep8
    '114	NM_000022.2	chr20	-	43248162	43280376	43248474	43280248	12	43248162,43248939,43249658,43251228,43251469,43251647,43252842,43254209,43255096,43257687,43264867,43280215,	43248488,43249042,43249788,43251293,43251571,43251719,43252970,43254325,43255240,43257810,43264929,43280376,	0	ADA	cmpl	cmpl	1,0,2,0,0,0,1,2,2,2,0,0,',  # nopep8
    '1440	NM_000038.5	chr5	+	112073555	112181936	112090587	112179823	16	112073555,112090569,112102022,112102885,112111325,112116486,112128142,112136975,112151191,112154662,112157592,112162804,112163625,112164552,112170647,112173249,	112073622,112090722,112102107,112103087,112111434,112116600,112128226,112137080,112151290,112155041,112157688,112162944,112163703,112164669,112170862,112181936,	0	APC	cmpl	cmpl	-1,0,0,1,2,0,0,0,0,0,1,1,0,0,0,2,',  # nopep8
    '184	NM_000492.3	chr7	+	117120016	117308718	117120148	117307162	27	117120016,117144306,117149087,117170952,117174329,117175301,117176601,117180153,117182069,117188694,117199517,117227792,117230406,117231987,117234983,117242879,117243585,117246727,117250572,117251634,117254666,117267575,117282491,117292895,117304741,117305512,117306961,	117120201,117144417,117149196,117171168,117174419,117175465,117176727,117180400,117182162,117188877,117199709,117227887,117230493,117232711,117235112,117242917,117243836,117246807,117250723,117251862,117254767,117267824,117282647,117292985,117304914,117305618,117308718,	0	CFTR	cmpl	cmpl	0,2,2,0,0,0,2,2,0,0,0,0,2,2,0,0,2,1,0,1,1,0,0,0,0,2,0,',  # nopep8
])


# Mock transcripts.
_transcripts = {}
for trans in map(make_transcript,
                 read_refgene(StringIO(_refgene))):
    _transcripts[trans.name] = trans
    _transcripts[trans.name + '.' + str(trans.version)] = trans


# Mock genome sequence.
_genome_seq = dict([
    (('chr11', 17496507, 17496508), 'T'),
    (('chr11', 17496477, 17496507), 'AGCAGCATGAAGGTCAGGATCCACCGCAGG'),
    (('chr11', 17498251, 17498252), 'G'),
    (('chr11', 17498221, 17498251), 'CACCACGTTGAGCGCGTCCACAAAGCAGCC'),
    (('chr11', 17418842, 17418843), 'G'),
    (('chr11', 17418812, 17418842), 'CTCCATGTCTGCCAGGTTCCTCACCATCCA'),
    (('chr11', 17464265, 17464266), 'C'),
    (('chr11', 17464235, 17464265), 'ACCCAGGGCTGGCTGTGTGGGGTGAACTCA'),
    (('chr11', 17464265, 17464266), 'C'),
    (('chr11', 17464235, 17464265), 'ACCCAGGGCTGGCTGTGTGGGGTGAACTCA'),
    (('chr11', 17464265, 17464266), 'C'),
    (('chr11', 17464235, 17464265), 'ACCCAGGGCTGGCTGTGTGGGGTGAACTCA'),
    (('chr11', 17452525, 17452526), 'T'),
    (('chr11', 17452495, 17452525), 'CCCACGAAAGTCTGTGGACAGAGGCACAAG'),
    (('chr11', 17452579, 17452580), 'C'),
    (('chr11', 17452549, 17452579), 'GGGTCCCTCCCACACTGGAAACGCTCAGCA'),
    (('chr11', 17450106, 17450107), 'C'),
    (('chr11', 17450076, 17450106), 'ATGCTGGGAGTAGCAAGGGGAGGCCGGGCA'),
    (('chr11', 17449509, 17449510), 'C'),
    (('chr11', 17449479, 17449509), 'CTCCCATGATCTTCATTAGGCGTGTCCCAC'),
    (('chr11', 17449410, 17449411), 'T'),
    (('chr11', 17449380, 17449410), 'AGGGGCATGCTGGAGGGGTGGACTGGGCCA'),
    (('chr11', 17449412, 17449413), 'C'),
    (('chr11', 17449382, 17449412), 'GGGCATGCTGGAGGGGTGGACTGGGCCATA'),
    (('chr11', 17449411, 17449412), 'A'),
    (('chr11', 17449381, 17449411), 'GGGGCATGCTGGAGGGGTGGACTGGGCCAT'),
    (('chr17', 7126027, 7126037), 'AGCAGAGGTG'),
    (('chr17', 7125997, 7126027), 'GAAGAAGATGGGCATCAAGGCTTCAAACAC'),
    (('chr17', 7126182, 7126185), 'CGG'),
    (('chr17', 7126152, 7126182), 'CAGGTACCATGAGAGGCATCATTGCTAAGG'),
    (('chr11', 108004587, 108004593), 'TTTTTT'),
    (('chr11', 108004557, 108004587), 'AGTAAGTGCTACAAGAACACCCATTGGATC'),
    (('chr20', 43251291, 43251296), 'ATCTG'),
    (('chr20', 43251261, 43251291), 'GCACCAGTGAGGTAGCTGGACCAGGGGCAG'),
    (('chr17', 48246449, 48246452), 'CAG'),
    (('chr17', 48246419, 48246449), 'ACTCTGCTGACAGTGACTTCTATCTGGTCC'),
    (('chr2', 241808282, 241808285), 'ATG'),
    (('chr2', 241808252, 241808282), 'GGTTCCCGAGCGGCAGGTTGGGTGCGGACC'),
    (('chr2', 241808282, 241808285), 'ATG'),
    (('chr2', 241808252, 241808282), 'GGTTCCCGAGCGGCAGGTTGGGTGCGGACC'),
    (('chr17', 41242937, 41242959), 'GCACACACACACACGCTTTTTA'),
    (('chr17', 41242907, 41242937), 'GAATGCAAAGGACACCACACACACGCATGT'),
    (('chr1', 76216233, 76216235), 'AA'),
    (('chr1', 76216203, 76216233), 'GAAAACTTTCGGAAAGCTACTTGTAGAGGT'),
    (('chr2', 241808727, 241808729), 'GG'),
    (('chr2', 241808697, 241808727), 'GTGCTGGAGCCTGGGGACTCCTTCCTGGTT'),
    (('chr5', 112163692, 112163694), 'AC'),
    (('chr5', 112163662, 112163692), 'TTGTGGCCCAACTAAAATCTGAAAGTGAAG'),
    (('chr5', 112175545, 112175547), 'AG'),
    (('chr5', 112175515, 112175545), 'CCATGCAGTGGAATGGTAAGTGGCATTATA'),
    (('chr5', 112175545, 112175547), 'AG'),
    (('chr5', 112175515, 112175545), 'CCATGCAGTGGAATGGTAAGTGGCATTATA'),
    (('chr1', 76199215, 76199222), 'GTCTTGG'),
    (('chr1', 76199185, 76199215), 'ATGTGTTGAAACATTTTGATACTGTAGGAG'),
    (('chr1', 76199231, 76199232), 'T'),
    (('chr1', 76199201, 76199231), 'TGATACTGTAGGAGGTCTTGGACTTGGAAC'),
    (('chr1', 76199267, 76199274), 'TGGATGT'),
    (('chr1', 76199237, 76199267), 'TGCTTGTTTAATTAGTGAAGAATTGGCTTA'),
    (('chr1', 76200516, 76200520), 'GAAG'),
    (('chr1', 76200486, 76200516), 'TATTATTGCTGGAAATGATCAACAAAAGAA'),
    (('chr1', 76200516, 76200520), 'GAAG'),
    (('chr1', 76200486, 76200516), 'TATTATTGCTGGAAATGATCAACAAAAGAA'),
    (('chr1', 76198411, 76198413), 'GA'),
    (('chr1', 76198381, 76198411), 'GAGGAAATCATCCCAGTGGCTGCAGAATAT'),
    (('chr1', 76198563, 76198564), 'C'),
    (('chr1', 76198533, 76198563), 'CTAGTATCCAGTCCCCCTAATTAGAAGAGC'),
    (('chr1', 76199312, 76199314), 'GG'),
    (('chr1', 76199282, 76199312), 'TCAGACTGCTATTGAAGGAAATTCTTTGGG'),
    (('chr1', 76205669, 76205671), 'TT'),
    (('chr1', 76205639, 76205669), 'TTTCTCTTGTTTTTATATATTCAAGGCTTA'),
    (('chr1', 76227048, 76227049), 'C'),
    (('chr1', 76227018, 76227048), 'TGTAGAAAAACTAATGAGGGATGCCAAAAT'),
    (('chr1', 76227050, 76227052), 'AT'),
    (('chr1', 76227020, 76227050), 'TAGAAAAACTAATGAGGGATGCCAAAATCT'),
    (('chr1', 76199231, 76199232), 'T'),
    (('chr1', 76199201, 76199231), 'TGATACTGTAGGAGGTCTTGGACTTGGAAC'),
    (('chr7', 117180217, 117180221), 'TTCT'),
    (('chr7', 117180187, 117180217), 'GTGAGATACTTCAATAGCTCAGCCTTCTTC'),
    (('chr7', 117171119, 117171121), 'CA'),
    (('chr7', 117171089, 117171119), 'CCTACACCCAGCCATTTTTGGCCTTCATCA'),
    (('chr7', 117176661, 117176664), 'TAT'),
    (('chr7', 117176631, 117176661), 'ACTTGTGATTACCTCAGAAATGATTGAAAA'),
    (('chr7', 117182106, 117182107), 'A'),
    (('chr7', 117182076, 117182106), 'TACAAAAGCAAGAATATAAGACATTGGAAT'),
    (('chr7', 117292909, 117292910), 'T'),
    (('chr7', 117292879, 117292909), 'TTCTTTTTTGCTATAGAAAGTATTTATTTT'),
    (('chr17', 41245340, 41245341), 'T'),
    (('chr17', 41245310, 41245340), 'TCTTCAGCATTATTAGACACTTTAACTGTT'),
    (('chr17', 41245340, 41245341), 'T'),
    (('chr17', 41245310, 41245340), 'TCTTCAGCATTATTAGACACTTTAACTGTT'),
    (('chr17', 41245340, 41245341), 'T'),
    (('chr17', 41245310, 41245340), 'TCTTCAGCATTATTAGACACTTTAACTGTT'),
    (('chr17', 41245340, 41245341), 'T'),
    (('chr17', 41245310, 41245340), 'TCTTCAGCATTATTAGACACTTTAACTGTT'),
    (('chr7', 117307164, 117307165), 'A'),
    (('chr7', 117307134, 117307164), 'AGAAGAGGTGCAAGATACAAGGCTTTAGAG'),
    (('chr11', 17496507, 17496508), 'T'),
    (('chr11', 17496477, 17496507), 'AGCAGCATGAAGGTCAGGATCCACCGCAGG'),
    (('chr17', 7126027, 7126037), 'AGCAGAGGTG'),
    (('chr17', 7125997, 7126027), 'GAAGAAGATGGGCATCAAGGCTTCAAACAC'),
    (('chr1', 76216233, 76216235), 'AA'),
    (('chr1', 76216203, 76216233), 'GAAAACTTTCGGAAAGCTACTTGTAGAGGT'),
    (('chr1', 76199213, 76199220), 'AGGTCTT'),
    (('chr1', 76199183, 76199213), 'CAATGTGTTGAAACATTTTGATACTGTAGG'),
    (('chr1', 76200411, 76200612), 'CATCTCTGAATTTACATATCCAATAAAAATGACTTGATTTTTTAATGTCAATTTTCTTCGGTAGCAAATGCCTATTATTATTGCTGGAAATGATCAACAAAAGAAGAAGTATTTGGGGAGAATGACTGAGGAGCCATTGATGTGTGTGAGTATGTGTAACTGCCGCTTTATTTCACACTTAAGAAGGGAACAAAGGTGCTA'),  # nopep8
    (('chr1', 76198312, 76198513), 'AATAATTTTCCCTTAGAGTTCACCGAACAGCAGAAAGAATTTCAAGCTACTGCTCGTAAATTTGCCAGAGAGGAAATCATCCCAGTGGCTGCAGAATATGATAAAACTGGTGAAGTAGGTATATACATTTTAAAGAGGGAAAAATCTTTTACATTTTTTACAAGATTATGTAATCAAACTATCTGGATTTCAAAATATATT'),  # nopep8
    (('chr1', 76198464, 76198665), 'ATTTTTTACAAGATTATGTAATCAAACTATCTGGATTTCAAAATATATTTTAACTCAGTTCTTTTTCTTCTAGTATCCAGTCCCCCTAATTAGAAGAGCCTGGGAACTTGGTTTAATGAACACACACATTCCAGAGAACTGTGGTAAGCTTTCTTTATATTTTTAATACTGGAATGCATATGAGTAAGAAAAATTGTGGAA'),  # nopep8
    (('chr1', 76198564, 76198565), 'T'),
    (('chr1', 76198565, 76198566), 'G'),
    (('chr1', 76199209, 76199410), 'TAGGAGGTCTTGGACTTGGAACTTTTGATGCTTGTTTAATTAGTGAAGAATTGGCTTATGGATGTACAGGGGTTCAGACTGCTATTGAAGGAAATTCTTTGGGGGTAAGTGACTTAGAAAATTAACTACCTAACTCAGCTCTTGTTAATGAGATAGTTACTCCTGAAGAAGTTGGATTTTTAGAAGAAAAAAAAAGGAAAG'),  # nopep8
    (('chr1', 76205569, 76205770), 'aataCAGCTTAAAAAATAAAACAATCCTGTTTCCAAACAGTCAAAATTTAATCACTAACATTTAATTTCATTTCTCTTGTTTTTATATATTCAAGGCTTATTGTGTAACAGAACCTGGAGCAGGCTCTGATGTAGCTGGTATAAAGACCAAAGCAGAAAAGAAAGGAGATGAGTATATTATTAATGGTCAGAAGATGTGGA'),  # nopep8
    (('chr1', 76226949, 76227150), 'TATTGCAAATCAGTTAGCTACTGATGCTGTGCAGATACTTGGAGGCAATGGATTTAATACAGAATATCCTGTAGAAAAACTAATGAGGGATGCCAAAATCTATCAGGTAAGGTTAAAGATGATTTTTTTGGTTTGCAAGGAGAGAGAATATTCATAAATGACAACGTGGATTTCTGATTATTTAGAAATTTTATGTCCCTA'),  # nopep8
    (('chr1', 76227049, 76227050), 'T'),
    (('chr1', 76227050, 76227051), 'A'),
    (('chr1', 76226951, 76227152), 'TTGCAAATCAGTTAGCTACTGATGCTGTGCAGATACTTGGAGGCAATGGATTTAATACAGAATATCCTGTAGAAAAACTAATGAGGGATGCCAAAATCTATCAGGTAAGGTTAAAGATGATTTTTTTGGTTTGCAAGGAGAGAGAATATTCATAAATGACAACGTGGATTTCTGATTATTTAGAAATTTTATGTCCCTAGT'),  # nopep8
    (('chr1', 76199132, 76199333), 'AAGCAATTAAAATAGTTTACCTTTATTTCTATTGTGATGTACTACATATTACAATGTGTTGAAACATTTTGATACTGTAGGAGGTCTTGGACTTGGAACTTTTGATGCTTGTTTAATTAGTGAAGAATTGGCTTATGGATGTACAGGGGTTCAGACTGCTATTGAAGGAAATTCTTTGGGGGTAAGTGACTTAGAAAATTA'),  # nopep8
    (('chr1', 76199231, 76199232), 'T'),
    (('chr1', 76199232, 76199233), 'T'),
    (('chr7', 117180110, 117180311), 'AATAAAATAACATCCTGAATTTTATTGTTATTGTTTTTTATAGAACAGAACTGAAACTGACTCGGAAGGCAGCCTATGTGAGATACTTCAATAGCTCAGCCTTCTTCTTCTCAGGGTTCTTTGTGGTGTTTTTATCTGTGCTTCCCTATGCACTAATCAAAGGAATCATCCTCCGGAAAATATTCACCACCATCTCATTCT'),  # nopep8
    (('chr7', 117171020, 117171221), 'GGAGGAACGCTCTATCGCGATTTATCTAGGCATAGGCTTATGCCTTCTCTTTATTGTGAGGACACTGCTCCTACACCCAGCCATTTTTGGCCTTCATCACATTGGAATGCAGATGAGAATAGCTATGTTTAGTTTGATTTATAAGAAGGTAATACTTCCTTGCACAGGCCCCATGGCACATATATTCTGTATCGTACATGT'),  # nopep8
    (('chr7', 117176560, 117176761), 'TACTATTAgattgattgattgattgattgattgattTACAGAGATCAGAGAGCTGGGAAGATCAGTGAAAGACTTGTGATTACCTCAGAAATGATTGAAAATATCCAATCTGTTAAGGCATACTGCTGGGAAGAAGCAATGGAAAAAATGATTGAAAACTTAAGACAGTAAGTTGTTCCAATAATTTCAATATTGTTAGTA'),  # nopep8
    (('chr7', 117182004, 117182205), 'AGATGTAATAATGCATTAATGCTATTCTGATTCTATAATATGTTTTTGCTCTCTTTTATAAATAGGATTTCTTACAAAAGCAAGAATATAAGACATTGGAATATAACTTAACGACTACAGAAGTAGTGATGGAGAATGTAACAGCCTTCTGGGAGGAGGTCAGAATTTTTAAAAAATTGTTTGCTCTAAACACCTAACTGT'),  # nopep8
    (('chr7', 117182107, 117182109), 'TA'),
    (('chr7', 117182109, 117182111), 'AC'),
    (('chr7', 117292805, 117293006), 'CTACTTGAAATATTTTACAATACAATAAGGGAAAAATAAAAAGTTATTTAAGTTATTCATACTTTCTTCTTCTTTTCTTTTTTGCTATAGAAAGTATTTATTTTTTCTGGAACATTTAGAAAAAACTTGGATCCCTATGAACAGTGGAGTGATCAAGAAATATGGAAAGTTGCAGATGAGGTAAGGCTGCTAACTGAAATG'),  # nopep8
    (('chr7', 117292910, 117292911), 'T'),
    (('chr7', 117292911, 117292912), 'C'),
])
