"""
Models for representing genomic elements.
"""

from __future__ import unicode_literals

from collections import namedtuple


class Position(object):
    """A position in the genome."""

    def __init__(self, chrom, chrom_start, chrom_stop, is_forward_strand):
        self.chrom = chrom
        self.chrom_start = chrom_start
        self.chrom_stop = chrom_stop
        self.is_forward_strand = is_forward_strand

    def __repr__(self):
        return "<Position %s[%d:%d]>" % (
            self.chrom, self.chrom_start, self.chrom_stop)


class Gene(object):
    def __init__(self, name):
        self.name = name


class Transcript(object):
    """RefGene Transcripts for hg19

    A gene may have multiple transcripts with different combinations of exons.
    """

    def __init__(self, name, version, gene, tx_position, cds_position,
                 is_default=False, exons=None):
        self.name = name
        self.version = version
        self.gene = Gene(gene)
        self.tx_position = tx_position
        self.cds_position = cds_position
        self.is_default = is_default
        self.exons = exons if exons else []

    @property
    def full_name(self):
        if self.version is not None:
            return '%s.%d' % (self.name, self.version)
        else:
            return self.name

    @property
    def is_coding(self):
        # Coding transcripts have CDS with non-zero length.
        return (self.cds_position.chrom_stop -
                self.cds_position.chrom_start > 0)

    @property
    def strand(self):
        return ('+' if self.tx_position.is_forward_strand else '-')

    @property
    def coding_exons(self):
        return [exon.get_as_interval(coding_only=True)
                for exon in self.exons]


BED6Interval_base = namedtuple(
    "BED6Interval_base", (
        "chrom",
        "chrom_start",
        "chrom_end",
        "name",
        "score",
        "strand"))


class BED6Interval(BED6Interval_base):
    def distance(self, offset):
        """Return the distance to the interval.

        if offset is inside the exon, distance is zero.
        otherwise, distance is the distance to the nearest edge.

        distance is positive if the exon comes after the offset.
        distance is negative if the exon comes before the offset.
        """

        start = self.chrom_start + 1
        end = self.chrom_end

        if start <= offset <= end:
            return 0

        start_distance = start - offset
        end_distance = offset - end

        if abs(start_distance) < abs(end_distance):
            return start_distance
        else:
            return -end_distance


class Exon(object):
    def __init__(self, transcript, tx_position, exon_number):
        self.transcript = transcript
        self.tx_position = tx_position
        self.exon_number = exon_number

    @property
    def get_exon_name(self):
        return "%s.%d" % (self.transcript.name, self.exon_number)

    def get_as_interval(self, coding_only=False):
        """Returns the coding region for this exon as a BED6Interval.

        This function returns a BED6Interval objects containing  position
        information for this exon. This may be used as input for
        pybedtools.create_interval_from_list() after casting chrom_start
        and chrom_end as strings.

        coding_only: only include exons in the coding region

        """

        exon_start = self.tx_position.chrom_start
        exon_stop = self.tx_position.chrom_stop

        # Get only exon coding region if requested
        if coding_only:
            if (exon_stop <= self.transcript.cds_position.chrom_start or
                    exon_start >= self.transcript.cds_position.chrom_stop):
                return None
            exon_start = max(exon_start,
                             self.transcript.cds_position.chrom_start)
            exon_stop = min(
                max(exon_stop, self.transcript.cds_position.chrom_start),
                self.transcript.cds_position.chrom_stop)

        return BED6Interval(
            self.tx_position.chrom,
            exon_start,
            exon_stop,
            self.get_exon_name,
            '.',
            self.strand,
        )

    @property
    def strand(self):
        strand = '+' if self.tx_position.is_forward_strand else '-'
        return strand
