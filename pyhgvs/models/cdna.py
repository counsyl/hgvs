import re

CDNA_START_CODON = 'cdna_start'
CDNA_STOP_CODON = 'cdna_stop'


class CDNACoord(object):
    """
    A HGVS cDNA-based coordinate.

    A cDNA coordinate can take one of these forms:

    N = nucleotide N in protein coding sequence (e.g. 11A>G)

    -N = nucleotide N 5' of the ATG translation initiation codon (e.g. -4A>G)
         NOTE: so located in the 5'UTR or 5' of the transcription initiation
         site (upstream of the gene, incl. promoter)

    *N = nucleotide N 3' of the translation stop codon (e.g. *6A>G)
         NOTE: so located in the 3'UTR or 3' of the polyA-addition site
         (including downstream of the gene)

    N+M = nucleotide M in the intron after (3' of) position N in the coding DNA
          reference sequence (e.g. 30+4A>G)

    N-M = nucleotide M in the intron before (5' of) position N in the coding
          DNA reference sequence (e.g. 301-2A>G)

    -N+M / -N-M = nucleotide in an intron in the 5'UTR (e.g. -45+4A>G)

    *N+M / *N-M = nucleotide in an intron in the 3'UTR (e.g. *212-2A>G)
    """

    def __init__(self, coord=0, offset=0, landmark=CDNA_START_CODON,
                 string=''):
        """
        coord: main coordinate along cDNA on the same strand as the transcript

        offset: an additional genomic offset from the main coordinate.  This
                allows referencing non-coding (e.g. intronic) positions.
                Offset is also interpreted on the coding strand.

        landmark: ('cdna_start', 'cdna_stop') indicating that 'coord'
                  is relative to one of these landmarks.

        string: a coordinate from an HGVS name.  If given coord, offset, and
                landmark should not be specified.
        """

        if string:
            if coord != 0 or offset != 0 or landmark != CDNA_START_CODON:
                raise ValueError("coord, offset, and landmark should not "
                                 "be given with string argument")

            self.parse(string)
        else:
            self.coord = coord
            self.offset = offset
            self.landmark = landmark

    def parse(self, coord_text):
        """
        Parse a HGVS formatted cDNA coordinate.
        """

        match = re.match(r"(|-|\*)(\d+)((-|\+)(\d+))?", coord_text)
        if not match:
            raise ValueError("unknown coordinate format '%s'" % coord_text)
        coord_prefix, coord, _, offset_prefix, offset = match.groups()

        self.coord = int(coord)
        self.offset = int(offset) if offset else 0

        if offset_prefix == '-':
            self.offset *= -1
        elif offset_prefix == '+' or offset is None:
            pass
        else:
            raise ValueError("unknown offset_prefix '%s'" % offset_prefix)

        if coord_prefix == '':
            self.landmark = CDNA_START_CODON
        elif coord_prefix == "-":
            self.coord *= -1
            self.landmark = CDNA_START_CODON
        elif coord_prefix == '*':
            self.landmark = CDNA_STOP_CODON
        else:
            raise ValueError("unknown coord_prefix '%s'" % coord_prefix)
        return self

    def __str__(self):
        """
        Return a formatted cDNA coordinate
        """
        if self.landmark == CDNA_STOP_CODON:
            coord_prefix = '*'
        else:
            coord_prefix = ''

        if self.offset < 0:
            offset = '%d' % self.offset
        elif self.offset > 0:
            offset = '+%d' % self.offset
        else:
            offset = ''

        return '%s%d%s' % (coord_prefix, self.coord, offset)

    def __eq__(self, other):
        """Equality operator."""
        return ((self.coord, self.offset, self.landmark) ==
                (other.coord, other.offset, other.landmark))

    def __repr__(self):
        """
        Returns a string representation of a cDNA coordinate.
        """
        if self.landmark != CDNA_START_CODON:
            return "CDNACoord(%d, %d, '%s')" % (
                self.coord, self.offset, self.landmark)
        else:
            return "CDNACoord(%d, %d)" % (self.coord, self.offset)
