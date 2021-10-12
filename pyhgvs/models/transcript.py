"""
Models for representing genomic elements.
"""

from __future__ import unicode_literals

from collections import namedtuple

from lazy import lazy

from pyhgvs.models.cdna import CDNA_START_CODON, CDNA_STOP_CODON, CDNACoord


class Gene(object):
    def __init__(self, name):
        self.name = name


class Transcript(object):
    """

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

    @lazy
    def ordered_exons(self):
        """Yield exons in coding order."""
        transcript_strand = self.tx_position.is_forward_strand
        if hasattr(self.exons, 'select_related'):
            exons = list(self.exons.select_related('tx_position'))
        else:
            exons = list(self.exons)
        exons.sort(key=lambda exon: exon.tx_position.chrom_start)
        if not transcript_strand:
            exons.reverse()
        return exons

    def get_coding_exons(self):
        """Yield non-empty coding exonic regions in coding order."""
        for exon in self.ordered_exons:
            region = exon.get_as_interval(coding_only=True)
            if region:
                yield region

    def get_utr5p_size(self):
        """Return the size of the 5prime UTR of a transcript."""

        transcript_strand = self.tx_position.is_forward_strand

        # Find the exon containing the start codon.
        start_codon = (self.cds_position.chrom_start if transcript_strand
                       else self.cds_position.chrom_stop - 1)
        cdna_len = 0
        for exon in self.ordered_exons:
            exon_start = exon.tx_position.chrom_start
            exon_end = exon.tx_position.chrom_stop
            if exon_start <= start_codon < exon_end:
                break
            cdna_len += exon_end - exon_start
        else:
            raise ValueError("transcript contains no exons")

        if transcript_strand:
            return cdna_len + (start_codon - exon_start)
        else:
            return cdna_len + (exon_end - start_codon - 1)

    def find_stop_codon(self, cds_position):
        """Return the position along the cDNA of the base after the stop codon."""
        if cds_position.is_forward_strand:
            stop_pos = cds_position.chrom_stop
        else:
            stop_pos = cds_position.chrom_start
        cdna_pos = 0
        for exon in self.ordered_exons:
            exon_start = exon.tx_position.chrom_start
            exon_stop = exon.tx_position.chrom_stop

            if exon_start <= stop_pos <= exon_stop:
                if cds_position.is_forward_strand:
                    return cdna_pos + stop_pos - exon_start
                else:
                    return cdna_pos + exon_stop - stop_pos
            else:
                cdna_pos += exon_stop - exon_start
        raise ValueError('Stop codon is not in any of the exons')

    def cdna_to_genomic_coord(self, coord):
        """Convert a HGVS cDNA coordinate to a genomic coordinate."""
        transcript_strand = self.tx_position.is_forward_strand
        utr5p = (self.get_utr5p_size()
                 if self.is_coding else 0)

        # compute starting position along spliced transcript.
        if coord.landmark == CDNA_START_CODON:
            if coord.coord > 0:
                pos = utr5p + coord.coord
            else:
                pos = utr5p + coord.coord + 1
        elif coord.landmark == CDNA_STOP_CODON:
            if coord.coord < 0:
                raise ValueError('CDNACoord cannot have a negative coord and '
                                 'landmark CDNA_STOP_CODON')
            pos = self.find_stop_codon(self.cds_position) + coord.coord
        else:
            raise ValueError('unknown CDNACoord landmark "%s"' % coord.landmark)

        # 5' flanking sequence.
        if pos < 1:
            if transcript_strand:
                return self.tx_position.chrom_start + pos
            else:
                return self.tx_position.chrom_stop - pos + 1

        # Walk along transcript until we find an exon that contains pos.
        cdna_start = 1
        cdna_end = 1
        for exon in self.ordered_exons:
            exon_start = exon.tx_position.chrom_start + 1
            exon_end = exon.tx_position.chrom_stop
            cdna_end = cdna_start + (exon_end - exon_start)
            if cdna_start <= pos <= cdna_end:
                break
            cdna_start = cdna_end + 1
        else:
            # 3' flanking sequence
            if transcript_strand:
                return self.cds_position.chrom_stop + coord.coord
            else:
                return self.cds_position.chrom_start + 1 - coord.coord

        # Compute genomic coordinate using offset.
        if transcript_strand:
            # Plus strand.
            return exon_start + (pos - cdna_start) + coord.offset
        else:
            # Minus strand.
            return exon_end - (pos - cdna_start) - coord.offset

    def genomic_to_cdna_coord(self, genomic_coord):
        """Convert a genomic coordinate to a cDNA coordinate and offset.
        """
        exons = [exon.get_as_interval()
                 for exon in self.ordered_exons]

        if len(exons) == 0:
            return None

        strand = self.strand

        if strand == "+":
            exons.sort(key=lambda exon: exon.chrom_start)
        else:
            exons.sort(key=lambda exon: -exon.chrom_end)

        distances = [exon.distance(genomic_coord)
                     for exon in exons]
        min_distance_to_exon = min(map(abs, distances))

        coding_offset = 0
        for exon in exons:
            exon_length = exon.chrom_end - exon.chrom_start
            distance = exon.distance(genomic_coord)
            if abs(distance) == min_distance_to_exon:
                if strand == "+":
                    exon_start_cds_offset = coding_offset + 1
                    exon_end_cds_offset = coding_offset + exon_length
                else:
                    exon_start_cds_offset = coding_offset + exon_length
                    exon_end_cds_offset = coding_offset + 1
                # This is the exon we want to annotate against.
                if distance == 0:
                    # Inside the exon.
                    if strand == "+":
                        coord = (exon_start_cds_offset +
                                 (genomic_coord -
                                  (exon.chrom_start + 1)))
                    else:
                        coord = (exon_end_cds_offset +
                                 (exon.chrom_end -
                                  genomic_coord))
                    cdna_coord = CDNACoord(coord, 0)
                else:
                    # Outside the exon.
                    if distance > 0:
                        nearest_exonic = exon_start_cds_offset
                    else:
                        nearest_exonic = exon_end_cds_offset
                    if strand == "+":
                        distance *= -1

                    # If outside transcript, don't use offset.
                    if (genomic_coord < self.tx_position.chrom_start + 1 or
                            genomic_coord > self.tx_position.chrom_stop):
                        nearest_exonic += distance
                        distance = 0
                    cdna_coord = CDNACoord(nearest_exonic, distance)
                break
            coding_offset += exon_length

        # Adjust coordinates for coding transcript.
        if self.is_coding:
            # Detect if position before start codon.
            utr5p = self.get_utr5p_size() if self.is_coding else 0
            cdna_coord.coord -= utr5p
            if cdna_coord.coord <= 0:
                cdna_coord.coord -= 1
            else:
                # Detect if position is after stop_codon.
                stop_codon = self.find_stop_codon(self.cds_position)
                stop_codon -= utr5p
                if (cdna_coord.coord > stop_codon or
                        cdna_coord.coord == stop_codon and cdna_coord.offset > 0):
                    cdna_coord.coord -= stop_codon
                    cdna_coord.landmark = CDNA_STOP_CODON

        return cdna_coord

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
