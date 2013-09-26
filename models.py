
from collections import namedtuple

from . import get_refseq_type


class Position(object):
    def __init__(self, chrom, chrom_start, chrom_stop, is_forward_strand):
        self.chrom = chrom
        self.chrom_start = chrom_start
        self.chrom_stop = chrom_stop
        self.is_forward_strand = is_forward_strand


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
    def is_coding(self):
        return get_refseq_type(self.name) == 'mRNA'

    @property
    def strand(self):
        return ('+' if self.tx_position.is_forward_strand else '-')

    @property
    def coding_exons(self):
        return [exon.get_as_interval(coding_only=True)
                for exon in self.exons]

    def genomic_offset_to_cds_coord(self, genomic_offset):
        """Convert a genomic coordinate to a cDNA coordinate and offset.
        """
        if self.is_coding:
            exons = [exon for exon in self.coding_exons if exon is not None]
        else:
            exons = [exon.get_as_interval()
                     for exon in self.exons]

        if len(exons) == 0:
            return None

        strand = self.strand

        if strand == "+":
            exons.sort(key=lambda exon: exon.chrom_start)
        else:
            exons.sort(key=lambda exon: -exon.chrom_end)

        distances = [exon.distance(genomic_offset)
                     for exon in exons]
        min_distance_to_exon = min(map(abs, distances))

        coding_offset = 0
        for exon in exons:
            exon_length = exon.chrom_end - exon.chrom_start
            distance = exon.distance(genomic_offset)
            if abs(distance) == min_distance_to_exon:
                if strand == "+":
                    exon_start_cds_offset = coding_offset + 1
                    exon_end_cds_offset = coding_offset + exon_length
                else:
                    exon_start_cds_offset = coding_offset + exon_length
                    exon_end_cds_offset = coding_offset + 1
                # this is the exon we want to annotate against
                if distance == 0:
                    # inside the exon
                    if strand == "+":
                        coord = (exon_start_cds_offset +
                                 (genomic_offset -
                                  (exon.chrom_start + 1)))
                    else:
                        coord = (exon_end_cds_offset +
                                 (exon.chrom_end -
                                  genomic_offset))
                    cds_coord = (coord, 0)
                else:
                    # outside the exon
                    if distance > 0:
                        nearest_coding = exon_start_cds_offset
                    else:
                        nearest_coding = exon_end_cds_offset
                    if strand == "+":
                        distance *= -1
                    cds_coord = (nearest_coding, distance)
                break
            coding_offset += exon_length

        return cds_coord

    def genomic_offset_to_cds_offset(self, genomic_offset):
        """Convert a genomic offset to a cDNA offset as in HGVS nomenclature.
        """
        coord, offset = self.genomic_offset_to_cds_coord(genomic_offset)
        if offset == 0:
            return "%d" % coord
        elif offset < 0:
            return "%d%d" % (coord, offset)
        else:
            return "%d+%d" % (coord, offset)


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
