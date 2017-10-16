"""
Helper functions.
"""

from __future__ import absolute_import
from __future__ import unicode_literals

from .models import Exon
from .models import Position
from .models import Transcript


def read_refgene(infile):
    """
    Iterate through a refGene file.

    GenePred extension format:
    http://genome.ucsc.edu/FAQ/FAQformat.html#GenePredExt

    Column definitions:
    0. uint undocumented id
    1. string name;             "Name of gene (usually transcript_id from GTF)"
    2. string chrom;                "Chromosome name"
    3. char[1] strand;              "+ or - for strand"
    4. uint txStart;                "Transcription start position"
    5. uint txEnd;                  "Transcription end position"
    6. uint cdsStart;               "Coding region start"
    7. uint cdsEnd;                 "Coding region end"
    8. uint exonCount;              "Number of exons"
    9. uint[exonCount] exonStarts;  "Exon start positions"
    10. uint[exonCount] exonEnds;   "Exon end positions"
    11. uint id;                    "Unique identifier"
    12. string name2;               "Alternate name (e.g. gene_id from GTF)"
    13. string cdsStartStat;        "enum('none','unk','incmpl','cmpl')"
    14. string cdsEndStat;          "enum('none','unk','incmpl','cmpl')"
    15. lstring exonFrames;         "Exon frame offsets {0,1,2}"
    """
    for line in infile:
        # Skip comments.
        if line.startswith('#'):
            continue
        row = line.rstrip('\n').split('\t')
        if len(row) != 16:
            raise ValueError(
                'File has incorrect number of columns '
                'in at least one line.')

        # Skip trailing ,
        exon_starts = list(map(int, row[9].split(',')[:-1]))
        exon_ends = list(map(int, row[10].split(',')[:-1]))
        exon_frames = list(map(int, row[15].split(',')[:-1]))
        exons = list(zip(exon_starts, exon_ends))

        yield {
            'chrom': row[2],
            'start': int(row[4]),
            'end': int(row[5]),
            'id': row[1],
            'strand': row[3],
            'cds_start': int(row[6]),
            'cds_end': int(row[7]),
            'gene_name': row[12],
            'exons': exons,
            'exon_frames': exon_frames
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
            transcript_json['strand'] == '+'))

    exons = transcript_json['exons']
    if not transcript.tx_position.is_forward_strand:
        exons = reversed(exons)

    for exon_number, (exon_start, exon_end) in enumerate(exons, 1):
        transcript.exons.append(
            Exon(transcript=transcript,
                 tx_position=Position(
                     transcript_json['chrom'],
                     exon_start,
                     exon_end,
                     transcript_json['strand'] == '+'),
                 exon_number=exon_number))

    return transcript


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
