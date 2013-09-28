
import nose

from ..import CDNACoord
from ..import CDNA_STOP_CODON
from ..import HGVSName


def test_parse_cdna_coord():
    """
    Parse cDNA coordinates.
    """
    for text, expected in _parse_cdna_coords:
        nose.tools.assert_equal(CDNACoord(string=text), expected)


def test_parse_name():
    """
    Test parsing of HGVS names.
    """
    for name, formatable, expected in _parse_names:
        hgvs_parsed = HGVSName(name)
        for key, value in expected.items():
            nose.tools.assert_equal(
                getattr(hgvs_parsed, key), value,
                (getattr(hgvs_parsed, key), value, name, expected))


def test_format_name():
    """
    Test parsing of HGVS names.
    """
    for expected_name, formatable, attrs in _parse_names:
        if formatable:
            name = HGVSName(**attrs).format()
            nose.tools.assert_equal(name, expected_name,
                                    (name, expected_name, attrs))


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
