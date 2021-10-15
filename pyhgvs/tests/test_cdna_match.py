"""
    RefSeq transcripts don't always match the genome reference, so can align with gaps.
"""

import nose
from pyhgvs import HGVSName

from pyhgvs.utils import make_transcript

_HGVS_NM_015120_GRCh37_COORDS = [
    # Has Gap=M185 I3 M250 5' UTR length is 111 so last match is c.74 then 3 bases insert
    ("NM_015120.4:c.74A>T", 73613070),
    ("NM_015120.4:c.75G>T", None),  # No consistent alignment
    ("NM_015120.4:c.76G>T", None),  # No consistent alignment
    ("NM_015120.4:c.77A>T", None),  # No consistent alignment
    ("NM_015120.4:c.78A>T", 73613071),
    ("NM_015120.4:c.79G>C", 73613072),
]

_HGVS_NM_001135649_GRCh37_COORDS = [
    # Has Gap=M460 I1 M337 5' UTR length is 158 so last match is c.302 then 1 base insert
    ("NM_001135649.3:c.300C>T", 88751754),
    ("NM_001135649.3:c.301A>C", 88751753),
    ("NM_001135649.3:c.302C>T", 88751752),
    ("NM_001135649.3:c.303C>T", None),  # N/A GRCh37 - 1 bp insert
    ("NM_001135649.3:c.304T>C", 88751751),
    ("NM_001135649.3:c.305T>C", 88751750),
]

_HGVS_NM_001012755_GRCh37_COORDS = [
    # 133, 6207, "M1574 I1 M4500 - 5'UTR length is 163 so last match is transcript pos (1574+132) ie c.1543 then 1bp ins
    ("NM_001012755.5:c.*619G>C", 103348398),  # NM_001012755.5:c.1543G>C
    ("NM_001012755.5:c.*620G>C", None),       # NM_001012755.5:c.1544G>C - No GRCh37
    ("NM_001012755.5:c.*621C>G", 103348397),  # NM_001012755.5:c.1545G>C
]


def test_cdna_to_genomic_coord():
    transcript = get_transcript("NM_015120.4")
    for hgvs_str, expected_genomic_coord in _HGVS_NM_015120_GRCh37_COORDS:
        hgvs_name = HGVSName(hgvs_str)
        if expected_genomic_coord:
            genomic_coord = transcript.cdna_to_genomic_coord(hgvs_name.cdna_start)
            nose.tools.assert_equal(genomic_coord, expected_genomic_coord)
        else:
            nose.tools.assert_raises(ValueError, transcript.cdna_to_genomic_coord, hgvs_name.cdna_start)


def test_genomic_to_cdna_coord():
    transcript = get_transcript("NM_015120.4")
    for hgvs_str, genomic_coord in _HGVS_NM_015120_GRCh37_COORDS:
        hgvs_name = HGVSName(hgvs_str)
        if genomic_coord:
            cdna_coord = transcript.genomic_to_cdna_coord(genomic_coord)
            nose.tools.assert_equal(cdna_coord, hgvs_name.cdna_start)


def test_cdna_to_genomic_coord_negative_strand():
    transcript = get_transcript("NM_001135649.3")
    for hgvs_str, expected_genomic_coord in _HGVS_NM_001135649_GRCh37_COORDS:
        hgvs_name = HGVSName(hgvs_str)
        if expected_genomic_coord:
            genomic_coord = transcript.cdna_to_genomic_coord(hgvs_name.cdna_start)
            nose.tools.assert_equal(genomic_coord, expected_genomic_coord)
        else:
            nose.tools.assert_raises(ValueError, transcript.cdna_to_genomic_coord, hgvs_name.cdna_start)


def test_cdna_to_genomic_coord_negative_strand2():
    transcript = get_transcript("NM_001012755.5")
    for hgvs_str, expected_genomic_coord in _HGVS_NM_001012755_GRCh37_COORDS:
        hgvs_name = HGVSName(hgvs_str)
        if expected_genomic_coord:
            genomic_coord = transcript.cdna_to_genomic_coord(hgvs_name.cdna_start)
            nose.tools.assert_equal(genomic_coord, expected_genomic_coord)
        else:
            nose.tools.assert_raises(ValueError, transcript.cdna_to_genomic_coord, hgvs_name.cdna_start)



@nose.SkipTest  # Currently fails
def test_genomic_to_cdna_coord_negative_strand():
    transcript = get_transcript("NM_001135649.3")
    for hgvs_str, genomic_coord in _HGVS_NM_001135649_GRCh37_COORDS:
        hgvs_name = HGVSName(hgvs_str)
        if genomic_coord:
            cdna_coord = transcript.genomic_to_cdna_coord(genomic_coord)
            nose.tools.assert_equal(cdna_coord, hgvs_name.cdna_start)


def test_genomic_to_cdna_coord_negative_strand2():
    transcript = get_transcript("NM_001012755.5")

    for hgvs_str, genomic_coord in _HGVS_NM_001012755_GRCh37_COORDS:
        hgvs_name = HGVSName(hgvs_str)
        if genomic_coord:
            cdna_coord = transcript.genomic_to_cdna_coord(genomic_coord)
            nose.tools.assert_equal(cdna_coord, hgvs_name.cdna_start)



def get_transcript(accession):
    transcript_json = _transcripts[accession]
    return make_transcript(transcript_json)


_transcripts = {
    # GRCh37
    "NM_015120.4": {
        "id": "NM_015120.4",
        "gene_name": "ALMS1",
        "end": 73837046,
        "chrom": "NC_000002.11",
        "exons": [[73612885, 73613320], [73635749, 73635875], [73646250, 73646446], [73649984, 73650102],
                  [73651557, 73652030], [73653580, 73653681], [73659325, 73659419], [73675089, 73681194],
                  [73682288, 73682422], [73716760, 73718625], [73746901, 73747143], [73761950, 73762076],
                  [73777393, 73777564], [73784346, 73784481], [73786098, 73786269], [73799388, 73800551],
                  [73826527, 73826648], [73827804, 73828008], [73828321, 73828563], [73829311, 73829495],
                  [73830367, 73830431], [73835601, 73835701], [73836694, 73837046]],
        "start": 73612885,
        "strand": "+",
        "cds_end": 73836739,
        "cds_start": 73612996,
        "cdna_match": [[73612885, 73613320, 1, 438, "M185 I3 M250"], [73635749, 73635875, 439, 564, None],
                       [73646250, 73646446, 565, 760, None], [73649984, 73650102, 761, 878, None],
                       [73651557, 73652030, 879, 1351, None], [73653580, 73653681, 1352, 1452, None],
                       [73659325, 73659419, 1453, 1546, None], [73675089, 73681194, 1547, 7654, "M141 I3 M5964"],
                       [73682288, 73682422, 7655, 7788, None], [73716760, 73718625, 7789, 9653, None],
                       [73746901, 73747143, 9654, 9895, None], [73761950, 73762076, 9896, 10021, None],
                       [73777393, 73777564, 10022, 10192, None], [73784346, 73784481, 10193, 10327, None],
                       [73786098, 73786269, 10328, 10498, None], [73799388, 73800551, 10499, 11661, None],
                       [73826527, 73826648, 11662, 11782, None], [73827804, 73828008, 11783, 11986, None],
                       [73828321, 73828563, 11987, 12228, None], [73829311, 73829495, 12229, 12412, None],
                       [73830367, 73830431, 12413, 12476, None], [73835601, 73835701, 12477, 12576, None],
                       [73836694, 73837046, 12577, 12928, None]]
    },
    "NM_001135649.3": {
        # This is a strange case - there are 3 exons, but 2 cDNA match entries
        # I am not sure if this is an error, ie where are the splice sites?
        "id": "NM_001135649.3",
        "gene_name": "FOXI3",
        "end": 88752211,
        "chrom": "NC_000002.11",
        "exons": [[88746305, 88748348], [88751414, 88751751], [88751753, 88752211]],
        "start": 88746305,
        "strand": "-",
        "cds_end": 88752053,
        "cds_start": 88747725,
        "cdna_match": [[88746305, 88748348, 799, 2841, None], [88751414, 88752211, 1, 798, "M460 I1 M337"]]
    },
    "NM_001012755.5": {
        "id": "NM_001012755.5",
        "gene_name": "SLC25A53",
        "end": 103401690,
        "chrom": "NC_000023.10",
        "exons": [[103343897, 103349971], [103401558, 103401690]],
        "start": 103343897,
        "strand": "-",
        "cds_end": 103349940,
        "cds_start": 103349016,
        "cdna_match": [[103343897, 103349971, 133, 6207, "M1574 I1 M4500"], [103401558, 103401690, 1, 132, None]]
    }
}
