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

    We need both exons and cdna_match as need to know exact exon boundaries to work out flanking
    """

    def __init__(self, name, version, gene, tx_position, cds_position,
                 is_default=False, exons=None, cdna_match=None):
        self.name = name
        self.version = version
        self.gene = Gene(gene)
        self.tx_position = tx_position
        self.cds_position = cds_position
        self.is_default = is_default
        self.exons = exons or []
        self.cdna_match = cdna_match or []

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

    @lazy
    def ordered_cdna_match(self):
        """ Coding order """
        transcript_strand = self.tx_position.is_forward_strand
        cdna_match = list(self.cdna_match)
        cdna_match.sort(key=lambda cm: cm.tx_position.chrom_start)
        if not transcript_strand:
            cdna_match.reverse()
        return cdna_match

    def get_utr5p_size(self):
        """Return the size of the 5prime UTR of a transcript."""

        transcript_strand = self.tx_position.is_forward_strand

        # Find the exon containing the start codon.
        start_codon = (self.cds_position.chrom_start if transcript_strand
                       else self.cds_position.chrom_stop - 1)
        cdna_offset = 0
        for cdna_match in self.ordered_cdna_match:
            start = cdna_match.tx_position.chrom_start
            end = cdna_match.tx_position.chrom_stop
            if start <= start_codon < end:
                # We're inside this match
                if transcript_strand:
                    position = start_codon - start
                else:
                    position = end - start_codon - 1
                return cdna_offset + position + cdna_match.get_offset(position)
            cdna_offset += cdna_match.length

        # Couldn't find it
        cdna_matches = []
        for cdna_match in self.ordered_cdna_match:
            start = cdna_match.tx_position.chrom_start
            end = cdna_match.tx_position.chrom_stop
            cdna_matches.append(f"{start}-{end}")
        cdna_match_summary = ", ".join(cdna_matches)
        raise ValueError("Couldn't find start_codon (%d) in cdna_match: {%s" % (start_codon, cdna_match_summary))

    def find_stop_codon(self, cds_position):
        """Return the position along the cDNA of the base after the stop codon."""
        if cds_position.is_forward_strand:
            stop_pos = cds_position.chrom_stop
        else:
            stop_pos = cds_position.chrom_start
        cdna_offset = 0
        for cdna_match in self.ordered_cdna_match:
            start = cdna_match.tx_position.chrom_start
            stop = cdna_match.tx_position.chrom_stop

            if start <= stop_pos <= stop:
                # We're inside this match
                if cds_position.is_forward_strand:
                    position = stop_pos - start
                else:
                    position = stop - stop_pos
                return cdna_offset + position + cdna_match.get_offset(position)
            else:
                cdna_offset += cdna_match.length
        raise ValueError('Stop codon is not in any of the exons')

    def cdna_to_genomic_coord(self, coord):
        """Convert a HGVS cDNA coordinate to a genomic coordinate."""
        transcript_strand = self.tx_position.is_forward_strand

        # compute starting position along spliced transcript.
        if coord.landmark == CDNA_START_CODON:
            utr5p = (self.get_utr5p_size()
                     if self.is_coding else 0)

            if coord.coord > 0:
                cdna_pos = utr5p + coord.coord
            else:
                cdna_pos = utr5p + coord.coord + 1
        elif coord.landmark == CDNA_STOP_CODON:
            if coord.coord < 0:
                raise ValueError('CDNACoord cannot have a negative coord and '
                                 'landmark CDNA_STOP_CODON')
            cdna_pos = self.find_stop_codon(self.cds_position) + coord.coord
        else:
            raise ValueError('unknown CDNACoord landmark "%s"' % coord.landmark)

        # 5' flanking sequence (no need to account for gaps)
        if cdna_pos < 1:
            if transcript_strand:
                return self.tx_position.chrom_start + cdna_pos
            else:
                return self.tx_position.chrom_stop - cdna_pos + 1

        # Walk along transcript until we find an exon that contains cdna_pos.
        for cdna_match in self.ordered_cdna_match:
            if cdna_match.cdna_start <= cdna_pos <= cdna_match.cdna_end:
                match_pos = cdna_pos - cdna_match.cdna_start
                match_pos -= cdna_match.get_offset(match_pos)
                # Compute genomic coordinate using offset.
                if transcript_strand:
                    # Plus strand.
                    start = cdna_match.tx_position.chrom_start + 1
                    return start + match_pos + coord.offset
                else:
                    # Minus strand.
                    end = cdna_match.tx_position.chrom_stop
                    return end - match_pos - coord.offset
        else:
            # 3' flanking sequence (no need to account for gaps)
            if transcript_strand:
                return self.cds_position.chrom_stop + coord.coord
            else:
                return self.cds_position.chrom_start + 1 - coord.coord

    def genomic_to_cdna_coord(self, genomic_coord):
        """ Convert a genomic coordinate to a cDNA coordinate and offset. """
        exons = [exon.get_as_interval()
                 for exon in self.ordered_cdna_match]

        if len(exons) == 0:
            return None

        strand = self.strand
        distances = [exon.distance(genomic_coord)
                     for exon in exons]
        min_distance_to_exon = min(map(abs, distances))
        if min_distance_to_exon:
            # We're outside of exon - so need to find closest point
            for exon in exons:
                distance = exon.distance(genomic_coord)
                if abs(distance) == min_distance_to_exon:

                    # Outside the exon.
                    if distance > 0:
                        genomic_nearest_exon = exon.chrom_start + 1
                    else:
                        genomic_nearest_exon = exon.chrom_end

                    if strand == "+":
                        distance *= -1

                    nearest_exon_coord = self._exon_genomic_to_cdna_coord(genomic_nearest_exon)

                    # If outside transcript, don't use offset.
                    if (genomic_coord < self.tx_position.chrom_start + 1 or
                            genomic_coord > self.tx_position.chrom_stop):
                        nearest_exon_coord += distance
                        distance = 0
                    cdna_coord = CDNACoord(nearest_exon_coord, distance)
                    break
            else:
                raise ValueError("Could not find closest exon!")  # Should never happen
        else:
            coord = self._exon_genomic_to_cdna_coord(genomic_coord)
            cdna_coord = CDNACoord(coord, 0)

        # Adjust coordinates for coding transcript.
        if self.is_coding:
            # Detect if position before start codon.
            utr5p = self.get_utr5p_size()
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

    def _exon_genomic_to_cdna_coord(self, genomic_coord):
        cdna_offset = 0
        for cdna_match in self.ordered_cdna_match:
            # Inside the exon.
            if cdna_match.tx_position.chrom_start <= genomic_coord <= cdna_match.tx_position.chrom_stop:
                if self.strand == "+":
                    position = genomic_coord - (cdna_match.tx_position.chrom_start + 1)
                else:
                    position = cdna_match.tx_position.chrom_stop - genomic_coord
                offset = cdna_match.get_offset(position, validate=False)
                return cdna_offset + position + offset + 1
            cdna_offset += cdna_match.length

        raise ValueError(f"Couldn't find {genomic_coord=}!")


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
    """ We still need exons to work out the flanking boundaries """
    def __init__(self, transcript, tx_position, number):
        self.transcript = transcript
        self.tx_position = tx_position
        self.number = number

    @property
    def name(self):
        return "%s.%d" % (self.transcript.name, self.number)

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
            self.name,
            '.',
            self.strand,
        )

    @property
    def strand(self):
        strand = '+' if self.tx_position.is_forward_strand else '-'
        return strand


class CDNA_Match(Exon):
    def __init__(self, transcript, tx_position, cdna_start, cdna_end, gap, number):
        super(CDNA_Match, self).__init__(transcript, tx_position, number)
        self.cdna_start = cdna_start
        self.cdna_end = cdna_end
        self.gap = gap

    @property
    def length(self):
        return self.cdna_end - self.cdna_start + 1

    def get_offset(self, position: int, validate=True):
        """ cdna_match GAP attribute looks like: 'M185 I3 M250' which is code/length
            @see https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md#the-gap-attribute
            codes operation
            M 	match
            I 	insert a gap into the reference sequence
            D 	insert a gap into the target (delete from reference)

            If you want the whole exon, then pass the end
        """

        if not self.gap:
            return 0

        position_1_based = position + 1
        cdna_match_index = 1
        offset = 0
        for gap_op in self.gap.split():
            code = gap_op[0]
            length = int(gap_op[1:])
            if code == "M":
                cdna_match_index += length
            elif code == "I":
                if validate and position_1_based < cdna_match_index + length:
                    raise ValueError("Coordinate (%d) inside insertion (%s) - no mapping possible!" % (position_1_based, gap_op))
                offset += length
            elif code == "D":
                if validate and position < cdna_match_index + length:
                    raise ValueError("Coordinate (%d) inside deletion (%s) - no mapping possible!" % (position_1_based, gap_op))
                offset -= length
            else:
                raise ValueError(f"Unknown code in cDNA GAP: {gap_op}")

            if cdna_match_index > position_1_based:
                break

        return offset
