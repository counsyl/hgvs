"""
Helper functions.
"""

from __future__ import absolute_import
from __future__ import unicode_literals

import operator
from .models.variants import Position
from .models.transcript import Transcript, CDNA_Match, Exon


def read_refgene(infile):
    """ refGene = genePred with extra column at front (and ignored ones after) """
    return read_genepred(infile, skip_first_column=True)


def read_genepred(infile, skip_first_column=False):
    """
    GenePred extension format:
    http://genome.ucsc.edu/FAQ/FAQformat.html#GenePredExt

    Column definitions:
    0. string name;                 "Name of gene (usually transcript_id from GTF)"
    1. string chrom;                "Chromosome name"
    2. char[1] strand;              "+ or - for strand"
    3. uint txStart;                "Transcription start position"
    4. uint txEnd;                  "Transcription end position"
    5. uint cdsStart;               "Coding region start"
    6. uint cdsEnd;                 "Coding region end"
    7. uint exonCount;              "Number of exons"
    8. uint[exonCount] exonStarts;  "Exon start positions"
    9. uint[exonCount] exonEnds;    "Exon end positions"
    10. uint id;                    "Unique identifier"
    11. string name2;               "Alternate name (e.g. gene_id from GTF)"
    """
    for line in infile:
        # Skip comments.
        if line.startswith('#'):
            continue
        row = line.rstrip('\n').split('\t')
        if skip_first_column:
            row = row[1:]

        # Skip trailing ,
        exon_starts = list(map(int, row[8].split(',')[:-1]))
        exon_ends = list(map(int, row[9].split(',')[:-1]))
        exons = list(zip(exon_starts, exon_ends))

        yield {
            'chrom': row[1],
            'start': int(row[3]),
            'end': int(row[4]),
            'id': row[0],
            'strand': row[2],
            'cds_start': int(row[5]),
            'cds_end': int(row[6]),
            'gene_name': row[11],
            'exons': exons,
        }


def make_transcript(transcript_json):
    """
    Make a Transcript form a JSON object.
    """

    transcript_name = transcript_json['id']
    if '.' in transcript_name:
        name, version = transcript_name.split('.')
    else:
        name, version = transcript_name, None

    transcript = Transcript(
        name=name,
        version=int(version) if version is not None else None,
        gene=transcript_json['gene_name'],
        tx_position=Position(
            transcript_json['chrom'],
            transcript_json['start'],
            transcript_json['end'],
            transcript_json['strand'] == '+'),
        cds_position=Position(
            transcript_json['chrom'],
            transcript_json['cds_start'],
            transcript_json['cds_end'],
            transcript_json['strand'] == '+'),
        start_codon_transcript_pos=transcript_json.get("start_codon_transcript_pos"),
        stop_codon_transcript_pos=transcript_json.get("stop_codon_transcript_pos"),
    )

    exons = transcript_json['exons']
    exons.sort(key=operator.itemgetter(0))
    cdna_match = transcript_json.get('cdna_match', [])
    cdna_match.sort(key=operator.itemgetter(0))

    if not transcript.tx_position.is_forward_strand:
        exons.reverse()
        cdna_match.reverse()

    # We don't use exons, but run everything through cDNA match so there's just 1 path
    # exons are treated as a perfect cDNA match
    if not cdna_match:
        cdna_match = json_perfect_exons_to_cdna_match(exons)

    for number, (exon_start, exon_end, cdna_start, cdna_end, gap) in enumerate(cdna_match, 1):
        transcript.cdna_match.append(CDNA_Match(transcript=transcript,
                                                tx_position=Position(
                                                    transcript_json['chrom'],
                                                    exon_start,
                                                    exon_end,
                                                    transcript_json['strand'] == '+'),
                                                cdna_start=cdna_start,
                                                cdna_end=cdna_end,
                                                gap=gap,
                                                number=number))

    return transcript


def json_perfect_exons_to_cdna_match(ordered_exons, single=False):
    """ Perfectly matched exons are basically a no-gap case of cDNA match
        single - use a single cDNA match (deletions for introns) - this is currently broken do not use
    """
    cdna_match = []
    if single:
        ordered_exons = list(ordered_exons)
        start = ordered_exons[0][0]
        end = ordered_exons[-1][1]
        last_exon_end = None
        gap_list = []
        cdna_length = 0
        for (exon_start, exon_end) in ordered_exons:
            # end up looking like "M D M D (M=exon, D=intron length)"
            if last_exon_end:
                intron_length = abs(exon_start - last_exon_end)
                gap_list.append("D%d" % intron_length)
            exon_length = exon_end - exon_start
            cdna_length += exon_length
            gap_list.append("M%d" % exon_length)
            last_exon_end = exon_end
        cdna_match = [[start, end, 1, cdna_length, " ".join(gap_list)]]
    else:
        cdna_start = 1
        for (exon_start, exon_end) in ordered_exons:
            exon_length = exon_end - exon_start
            cdna_end = cdna_start + exon_length - 1
            cdna_match.append([exon_start, exon_end, cdna_start, cdna_end, None])
            cdna_start = cdna_end + 1
    return cdna_match


def read_transcripts(refgene_file):
    """
    Read all transcripts in a RefGene file.
    """
    transcripts = {}
    for trans in (make_transcript(record)
                  for record in read_refgene(refgene_file)):
        transcripts[trans.name] = trans
        transcripts[trans.full_name] = trans

    return transcripts
