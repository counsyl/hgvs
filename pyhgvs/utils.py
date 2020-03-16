"""
Helper functions.
"""

from __future__ import absolute_import
from __future__ import unicode_literals

from .models import Exon
from .models import Position
from .models import Transcript


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
