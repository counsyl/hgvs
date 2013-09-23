"""
Methods for manipulating HGVS names

HGVS language:

HGVS = ALLELE
     | PREFIX_NAME : ALLELE

PREFIX_NAME = TRANSCRIPT
            | TRANSCRIPT '(' GENE ')'
            | GENE '{' TRANSCRIPT '}' # Counsyl style

TRANSCRIPT = TRANSCRIPT_NAME
           | TRANSCRIPT_NAME '.' TRANSCRIPT_VERSION

TRANSCRIPT_VERSION = NUMBER

ALLELE = 'c.' CDNA_ALLELE    # cDNA
       | 'g.' GENOMIC_ALLELE # genomic
       | 'm.' MIT_ALLELE     # mitochondrial sequence
       | 'n.' NC_ALLELE      # non-coding RNA reference sequence
       | 'r.' RNA_ALLELE     # RNA sequence (like r.76a>u)
       | 'p.' PROTEIN_ALLELE # protein sequence (like  p.Lys76Asn)

NC_ALLELE =
RNA_ALLELE =
CDNA_ALLELE = CDNA_COORD BASE '>' BASE                 # substitution
            | CDNA_COORD 'ins' BASE                    # 1bp insertion
            | CDNA_COORD 'del' BASE                    # 1bp deletion
            | CDNA_COORD 'dup' BASE                    # 1bp duplication
            | CDNA_COORD_RANGE 'del' BASES             # deletion
            | CDNA_COORD_RANGE 'ins' BASES             # insertion
            | CDNA_COORD_RANGE 'dup' BASES             # duplication
            | CDNA_COORD 'del' BASE 'ins' BASE         # 1bp indel
            | CDNA_COORD_RANGE 'del' BASES 'ins' BASES # indel
            | CDNA_COORD_RANGE 'delins' BASES          # indel

GENOMIC_ALLELE =
MIT_ALLELE = COORD BASE '>' BASE                 # substitution
           | COORD 'ins' BASE                    # 1bp insertion
           | COORD 'del' BASE                    # 1bp deletion
           | COORD 'dup' BASE                    # 1bp duplication
           | COORD_RANGE 'del' BASES             # deletion
           | COORD_RANGE 'ins' BASES             # insertion
           | COORD_RANGE 'dup' BASES             # duplication
           | COORD 'del' BASE 'ins' BASE         # 1bp indel
           | COORD_RANGE 'del' BASES 'ins' BASES # indel
           | COORD_RANGE 'delins' BASES          # indel

PROTEIN_ALLELE = TODO...

# A genomic range:
COORD_RANGE = COORD '_' COORD

# A cDNA range:
CDNA_COORD_RANGE = CDNA_COORD '_' CDNA_COORD

# A cDNA coordinate:
CDNA_COORD = COORD_PREFIX COORD
           | COORD_PREFIX COORD OFFSET_PREFIX OFFSET
COORD_PREFIX = '' | '-' | '*'
COORD = NUMBER
OFFSET_PREFIX = '-' | '+'
OFFSET = NUMBER

# Primatives:
NUMBER = \d+
BASE = [ACGT]
BASES = BASE+

"""

import re

from .variants import justify_indel
from .variants import normalize_variant
from .variants import revcomp


CDNA_START_CODON = 'cdna_start'
CDNA_STOP_CODON = 'cdna_stop'


class HGVSRegex(object):
    """
    All regular expression for HGVS names.
    """

    # DNA syntax
    BASE = "[acgtACGT]|\d+"
    BASES = "[acgtACGT]+|\d+"
    DNA_REF = "(?P<ref>" + BASES + ")"
    DNA_ALT = "(?P<alt>" + BASES + ")"

    # Mutation types
    SUB = "(?P<mutation_type>>)"
    INS = "(?P<mutation_type>ins)"
    DEL = "(?P<mutation_type>del)"
    DUP = "(?P<mutation_type>dup)"

    # Simple coordinate syntax
    COORD_START = "(?P<start>\d+)"
    COORD_END = "(?P<end>\d+)"
    COORD_RANGE = COORD_START + "_" + COORD_END

    # cDNA coordinate syntax
    CDNA_COORD = ("(?P<coord_prefix>|-|\*)(?P<coord>\d+)"
                  "((?P<offset_prefix>-|\+)(?P<offset>\d+))?")
    CDNA_START = ("(?P<start>(?P<start_coord_prefix>|-|\*)(?P<start_coord>\d+)"
                  "((?P<start_offset_prefix>-|\+)(?P<start_offset>\d+))?)")
    CDNA_END = (r"(?P<end>(?P<end_coord_prefix>|-|\*)(?P<end_coord>\d+)"
                "((?P<end_offset_prefix>-|\+)(?P<end_offset>\d+))?)")
    CDNA_RANGE = CDNA_START + "_" + CDNA_END

    # cDNA allele syntax
    CDNA_ALLELE = [
        # Substitution
        CDNA_START + DNA_REF + SUB + DNA_ALT,

        # 1bp insertion, deletion, duplication
        CDNA_START + INS + DNA_ALT,
        CDNA_START + DEL + DNA_REF,
        CDNA_START + DUP + DNA_REF,

        # Insertion, deletion, duplication
        CDNA_RANGE + INS + DNA_ALT,
        CDNA_RANGE + DEL + DNA_REF,
        CDNA_RANGE + DUP + DNA_REF,
        CDNA_RANGE + DEL,

        # Indels
        "(?P<delins>" + CDNA_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + CDNA_RANGE + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + CDNA_RANGE + 'delins' + DNA_ALT + ")",
    ]

    CDNA_ALLELE_REGEXES = [re.compile("^" + regex + "$")
                           for regex in CDNA_ALLELE]

    # Peptide syntax
    PEP = "([A-Z]([a-z]{2}))+"
    PEP_REF = "(?P<ref>" + PEP + ")"
    PEP_REF2 = "(?P<ref2>" + PEP + ")"
    PEP_ALT = "(?P<alt>" + PEP + ")"

    PEP_EXTRA = "(?P<extra>(|=|\?)(|fs))"

    # Peptide allele syntax
    PEP_ALLELE = [
        # No peptide change
        # Example: Glu1161=
        PEP_REF + COORD_START + PEP_EXTRA,

        # Peptide change
        # Example: Glu1161Ser
        PEP_REF + COORD_START + PEP_ALT + PEP_EXTRA,

        # Peptide indel
        # Example: Glu1161_Ser1164?fs
        PEP_REF + COORD_START + "_" + PEP_REF2 + COORD_END + PEP_EXTRA,
        PEP_REF + COORD_START + "_" + PEP_REF2 + COORD_END + PEP_ALT +
        PEP_EXTRA,
    ]

    PEP_ALLELE_REGEXES = [re.compile("^" + regex + "$")
                          for regex in PEP_ALLELE]

    # Genomic allele syntax
    GENOMIC_ALLELE = [
        # Substitution
        COORD_START + DNA_REF + SUB + DNA_ALT,

        # 1bp insertion, deletion, duplication
        COORD_START + INS + DNA_ALT,
        COORD_START + DEL + DNA_REF,
        COORD_START + DUP + DNA_REF,

        # Insertion, deletion, duplication
        COORD_RANGE + INS + DNA_ALT,
        COORD_RANGE + DEL + DNA_REF,
        COORD_RANGE + DUP + DNA_REF,
        COORD_RANGE + DEL,

        # Indels
        "(?P<delins>" + COORD_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + COORD_RANGE + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + COORD_RANGE + 'delins' + DNA_ALT + ")",
    ]

    GENOMIC_ALLELE_REGEXES = [re.compile("^" + regex + "$")
                              for regex in GENOMIC_ALLELE]


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
        elif self.coord < 0:
            coord_prefix = '-'
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
                (other.coord, other.offset, other.lankmark))

    def __repr__(self):
        """
        Returns a string representation of a cDNA coordinate.
        """
        if self.landmark != CDNA_START_CODON:
            return "CDNACoord(%d, %d, '%s')" % (
                self.coord, self.offset, self.landmark)
        else:
            return "CDNACoord(%d, %d)" % (self.coord, self.offset)


def get_exons(transcript):
    """Yield exons in coding order."""
    transcript_strand = transcript.tx_position.is_forward_strand
    exons = list(transcript.exons)
    exons.sort(key=lambda exon: exon.exon_number)
    #if transcript_strand:
    #    exons = (transcript.exons.select_related('tx_position')
    #             .order_by('tx_position__chrom_start'))
    #else:
    #    exons = (transcript.exons.select_related('tx_position')
    #             .order_by('-tx_position__chrom_start'))
    return exons


def get_coding_exons(transcript):
    """Yield non-empty coding exonic regions in coding order."""
    for exon in get_exons(transcript):
        region = exon.get_as_interval(coding_only=True)
        if region:
            yield region


def get_utr5p_size(transcript):
    """Return the size of the 5prime UTR of a transcript."""

    transcript_strand = transcript.tx_position.is_forward_strand
    exons = get_exons(transcript)

    # Find the exon containing the start codon.
    start_codon = (transcript.cds_position.chrom_start if transcript_strand
                   else transcript.cds_position.chrom_stop)
    cdna_len = 0
    for exon in exons:
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
        return cdna_len + (exon_end - start_codon)


def get_genomic_sequence(genome, chrom, start, end):
    """
    Return a sequence for the genomic region

    start, end: 0-based, end-exclusive coordinates of the sequence.
    """
    if start > end:
        return ''
    else:
        return str(genome[str(chrom)][start - 1:end]).upper()


def get_cdna_genomic_coordinate(transcript, coord):
    """Convert a HGVS cDNA coordinate to a genomic coordinate."""

    # TODO: still need to implement landmark handling (after stop codon)

    transcript_strand = transcript.tx_position.is_forward_strand
    exons = get_exons(transcript)
    utr5p = get_utr5p_size(transcript)

    # compute position along spliced transcript
    if coord.coord > 0:
        pos = utr5p + coord.coord
    else:
        pos = utr5p + coord.coord + 1

    # Walk along transcript until we find an exon that contains pos.
    cdna_start = 1
    cdna_end = 1
    for exon in exons:
        exon_start = exon.tx_position.chrom_start + 1
        exon_end = exon.tx_position.chrom_stop
        cdna_end = cdna_start + (exon_end - exon_start)
        if cdna_start <= pos <= cdna_end:
            break
        cdna_start = cdna_end + 1
    else:
        raise ValueError("coordinate not within an exon")

    # Compute genomic coordinate.
    if transcript_strand:
        # Plus strand.
        return exon_start + (pos - cdna_start) + coord.offset
    else:
        # Minus strand.
        return exon_end - (pos - cdna_start) - coord.offset


def get_allele(hgvs, genome, transcript=None):
    """Get an allele from a HGVSName, a genome, and a transcript."""
    chrom, start, end = hgvs.get_coords(transcript)
    _, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True)
    ref = get_genomic_sequence(genome, chrom, start, end)
    return chrom, start, end, ref, alt


def get_vcf_allele(hgvs, genome, transcript=None):
    """Get an VCF-style allele from a HGVSName, a genome, and a transcript."""
    chrom, start, end = hgvs.get_vcf_coords(transcript)
    _, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True)
    ref = get_genomic_sequence(genome, chrom, start, end)

    if hgvs.mutation_type in {'ins', 'del', 'dup', 'delins'}:
        # Left-pad alternate allele.
        alt = ref[0] + alt

    return chrom, start, end, ref, alt


def matches_ref_allele(hgvs, genome, transcript=None):
    """Return True if reference allele matches genomic sequence."""
    ref, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True)
    chrom, start, end = hgvs.get_coords(transcript)
    genome_ref = get_genomic_sequence(genome, chrom, start, end)
    return genome_ref == ref


class InvalidHGVSName(ValueError):
    def __init__(self, name='', part='name'):
        if name:
            message = 'Invalid HGVS %s "%s"' % (part, name)
        else:
            message = 'Invalid HGVS %s' % part
        super(InvalidHGVSName, self).__init__(message)

        self.name = name
        self.part = part


class HGVSName(object):
    """
    Represents a HGVS variant name.

    http://www.hgvs.org/mutnomen/standards.html
    """

    def __init__(self, name='', chrom='', transcript='', gene='', kind='',
                 mutation_type=None, start=0, end=0, ref_allele='',
                 ref2_allele='', alt_allele='',
                 cdna_start=None, cdna_end=None, pep_extra=''):

        # Full HGVS name.
        self.name = name

        # Name parts.
        self.chrom = chrom
        self.transcript = transcript
        self.gene = gene
        self.kind = kind
        self.mutation_type = mutation_type
        self.start = start
        self.end = end
        self.ref_allele = ref_allele    # reference allele
        self.ref2_allele = ref2_allele  # reference allele at end of pep indel
        self.alt_allele = alt_allele    # alternate allele

        # cDNA-specific fields
        self.cdna_start = cdna_start if cdna_start else CDNACoord()
        self.cdna_end = cdna_start if cdna_end else CDNACoord()

        # Protein-specific fields
        self.pep_extra = pep_extra

        if name:
            self.parse(name)

    def parse(self, name):
        """Parse a HGVS name."""
        try:
            # Does HGVS name have transcript/gene prefix?
            if ':' in name:
                prefix, allele = name.split(':')
            else:
                prefix = ''
                allele = name
        except:
            raise InvalidHGVSName(name)

        self.name = name

        # Parse prefix and allele.
        self.parse_allele(allele)
        if self.kind == 'g':
            self.chrom = prefix
        else:
            self.parse_transcript(prefix)

    def parse_transcript(self, prefix):
        """
        Parse a HGVS trancript/gene prefix.

        Some examples of full hgvs names with transcript include:
          NM_007294.3:c.2207A>C
          NM_007294.3(BRCA1):c.2207A>C
          BRCA1{NM_007294.3}:c.2207A>C
        """

        # No prefix.
        if prefix == '':
            self.transcript = ''
            self.gene = ''
            return

        # Transcript and gene given with parens:
        # example: NM_007294.3(BRCA1):c.2207A>C
        match = re.match("^(?P<transcript>[^(]+)\((?P<gene>[^)]+)\)$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = match.group('gene')
            return

        # Transcript and gene given with braces:
        # example: BRCA1{NM_007294.3}:c.2207A>C
        match = re.match("^(?P<gene>[^{]+)\{(?P<transcript>[^}]+)\}$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = match.group('gene')
            return

        # Only transcript:
        match = re.match("^(?P<transcript>.+)$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = ''
            return

        raise InvalidHGVSName(prefix, 'prefix')

    def parse_allele(self, allele):
        """
        Parse a HGVS allele description.

        Some examples include:
          cDNA substitution: c.101A>C,
          cDNA indel: c.3428delCinsTA, c.1000_1003delATG, c.1000_1001insATG
          No protein change: p.Glu1161=
          Protein change: p.Glu1161Ser
          Protein frameshift: p.Glu1161_Ser1164?fs
          Genomic substitution: g.1000100A>T
          Genomic indel: g.1000100_1000102delATG
        """
        if '.' not in allele:
            InvalidHGVSName(allele, 'allele')

        # Determine HGVS name kind.
        kind, details = allele.split('.', 1)
        self.kind = kind
        self.mutation_type = None

        if kind == "c":
            self.parse_cdna(details)
        elif kind == "p":
            self.parse_protein(details)
        elif kind == "g":
            self.parse_genome(details)
        else:
            raise NotImplementedError("not implemented: %s" % allele)

    def parse_cdna(self, details):
        """
        Parse a HGVS cDNA name.

        Some examples include:
          Substitution: 101A>C,
          Indel: 3428delCinsTA, 1000_1003delATG, 1000_1001insATG
        """

        for regex in HGVSRegex.CDNA_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups = match.groupdict()

                # Parse mutation type.
                if groups.get('delins'):
                    self.mutation_type = 'delins'
                else:
                    self.mutation_type = groups['mutation_type']

                # Parse coordinates.
                self.cdna_start = CDNACoord(string=groups.get('start'))
                if groups.get('end'):
                    self.cdna_end = CDNACoord(string=groups.get('end'))
                else:
                    self.cdna_end = CDNACoord(string=groups.get('start'))

                # Parse alleles.
                self.ref_allele = groups.get('ref', '')
                self.alt_allele = groups.get('alt', '')

                # Convert numerical allelles.
                if self.ref_allele.isdigit():
                    self.ref_allele = "N" * int(self.ref_allele)
                if self.alt_allele.isdigit():
                    self.alt_allele = "N" * int(self.alt_allele)

                # Convert duplication alleles.
                if self.mutation_type == "dup":
                    self.alt_allele = self.ref_allele * 2
                return

        raise InvalidHGVSName(details, 'cDNA allele')

    def parse_protein(self, details):
        """
        Parse a HGVS protein name.

        Some examples include:
          No change: Glu1161=
          Change: Glu1161Ser
          Frameshift: Glu1161_Ser1164?fs
        """
        for regex in HGVSRegex.PEP_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups = match.groupdict()

                # Parse mutation type.
                if groups.get('delins'):
                    self.mutation_type = 'delins'
                else:
                    self.mutation_type = '>'

                # Parse coordinates.
                self.start = int(groups.get('start'))
                if groups.get('end'):
                    self.end = int(groups.get('end'))
                else:
                    self.end = self.start

                # Parse alleles.
                self.ref_allele = groups.get('ref', '')
                if groups.get('ref2'):
                    self.ref2_allele = groups.get('ref2')
                    self.alt_allele = groups.get('alt', '')
                else:
                    # If alt is not given, assume matching with ref
                    self.ref2_allele = self.ref_allele
                    self.alt_allele = groups.get(
                        'alt', self.ref_allele)

                self.pep_extra = groups.get('extra')
                return

        raise InvalidHGVSName(details, 'protein allele')

    def parse_genome(self, details):
        """
        Parse a HGVS genomic name.

        Som examples include:
          Substitution: 1000100A>T
          Indel: 1000100_1000102delATG
        """

        for regex in HGVSRegex.GENOMIC_ALLELE_REGEXES:
            match = re.match(regex, details)
            if match:
                groups = match.groupdict()

                # Parse mutation type.
                if groups.get('delins'):
                    self.mutation_type = 'delins'
                else:
                    self.mutation_type = groups['mutation_type']

                # Parse coordinates.
                self.start = int(groups.get('start'))
                if groups.get('end'):
                    self.end = int(groups.get('end'))
                else:
                    self.end = self.start

                # Parse alleles.
                self.ref_allele = groups.get('ref', '')
                self.alt_allele = groups.get('alt', '')

                # Convert numerical allelles.
                if self.ref_allele.isdigit():
                    self.ref_allele = "N" * int(self.ref_allele)
                if self.alt_allele.isdigit():
                    self.alt_allele = "N" * int(self.alt_allele)

                # Convert duplication alleles.
                if self.mutation_type == "dup":
                    self.alt_allele = self.ref_allele * 2
                return

        raise InvalidHGVSName(details, 'genomic allele')

    def __repr__(self):
        try:
            return "HGVSName(%s)" % self.format()
        except NotImplementedError:
            return "HGVSName(%s)" % self.name

    def __unicode__(self):
        return self.format()

    def format(self, use_gene=True, use_counsyl=True):
        """Generate a HGVS name as a string."""

        if self.kind == 'c':
            allele = 'c.' + self.format_cdna()
        elif self.kind == "p":
            allele = 'p.' + self.format_protein()
        elif self.kind == "g":
            allele = 'g.' + self.format_genome()
        else:
            raise NotImplementedError("not implemented: '%s'" % self.kind)

        return self.format_transcript(
            use_gene=use_gene, use_counsyl=use_counsyl) + allele

    def format_transcript(self, use_gene=True, use_counsyl=True):
        """
        Generate HGVS trancript/gene prefix.

        Some examples of full hgvs names with transcript include:
          NM_007294.3:c.2207A>C
          NM_007294.3(BRCA1):c.2207A>C
          BRCA1{NM_007294.3}:c.2207A>C
        """
        if not self.transcript:
            return ''

        elif use_gene and self.gene:
            if use_counsyl:
                return '%s{%s}' % (self.gene, self.transcript)
            else:
                return '%s(%s)' % (self.transcript, self.gene)
        else:
            return self.transcript

    def format_cdna_coords(self):
        """
        Generate HGVS cDNA coordinates string.
        """
        # Format coordinates.
        if self.cdna_start == self.cdna_end:
            return str(self.cdna_start)
        else:
            return "%s_%s" % (self.cdna_start, self.cdna_end)

    def format_dna_allele(self):
        """
        Generate HGVS DNA allele.
        """
        if self.mutation_type == '>':
            # SNP.
            # example: 101A>C
            return self.ref_allele + '>' + self.alt_allele

        elif self.mutation_type == 'delins':
            # Indel.
            # example: 112_117delAGGTCAinsTG, 112_117delinsTG
            return 'del' + self.ref_allele + 'ins' + self.alt_allele

        elif self.mutation_type in ('del', 'dup'):
            # Delete, duplication.
            # example: 1000_1003delATG, 1000_1003dupATG
            return self.mutation_type + self.ref_allele

        elif self.mutation_type == 'ins':
            # Insert.
            # example: 1000_1001insATG
            return self.mutation_type + self.alt_allele

        else:
            raise AssertionError(
                "unknown mutation type: '%s'" % self.mutation_type)

    def format_cdna(self):
        """
        Generate HGVS cDNA allele.

        Some examples include:
          Substitution: 101A>C,
          Indel: 3428delCinsTA, 1000_1003delATG, 1000_1001insATG
        """
        return self.format_cdna_coords() + self.format_dna_allele()

    def format_protein(self):
        """
        Generate HGVS protein name.

        Some examples include:
          No change: Glu1161=
          Change: Glu1161Ser
          Frameshift: Glu1161_Ser1164?fs
        """
        if (self.start == self.end and
                self.ref_allele == self.ref2_allele ==
                self.alt_allele):
            # Match.
            # Example: Glu1161=
            return (self.ref_allele + str(self.start) + '=' +
                    self.pep_extra)

        elif (self.start == self.end and
                self.ref_allele == self.ref2_allele and
                self.ref_allele != self.alt_allele):
            # Change.
            # Example: Glu1161Ser
            return (self.ref_allele + str(self.start) +
                    self.pep_start_alt_allele + self.pep_extra)

        elif self.start != self.end:
            # Range change.
            # Example: Glu1161_Ser1164?fs
            return (self.ref_allele + str(self.start) + '_' +
                    self.pep_start_alt_allele + str(self.end) +
                    self.pep_extra)

        else:
            raise NotImplementedError('protein name formatting.')

    def format_coords(self):
        """
        Generate HGVS cDNA coordinates string.
        """
        # Format coordinates.
        if self.start == self.end:
            return str(self.start)
        else:
            return "%s_%s" % (self.start, self.end)

    def format_genome(self):
        """
        Generate HGVS genomic allele.

        Som examples include:
          Substitution: 1000100A>T
          Indel: 1000100_1000102delATG
        """
        return self.format_coords() + self.format_dna_allele()

    def get_coords(self, transcript=None):
        """Return genomic coordinates of reference allele."""
        if self.kind == 'c':
            chrom = transcript.tx_position.chrom
            start = get_cdna_genomic_coordinate(transcript, self.cdna_start)
            end = get_cdna_genomic_coordinate(transcript, self.cdna_end)

            if not transcript.tx_position.is_forward_strand:
                if end > start:
                    raise AssertionError(
                        "cdna_start cannot be greater than cdna_end")
                start, end = end, start

            if self.mutation_type == "ins":
                # Inserts have empty interval.
                if start < end:
                    start += 1
                    end -= 1
                else:
                    end = start - 1

            elif self.mutation_type == "dup":
                end = start - 1

        elif self.kind == 'g':
            chrom = self.chrom
            start = self.start
            end = self.end

        else:
            raise NotImplementedError(
                'Coordinates are available for this kind of HGVS name "%s"'
                % self.kind)

        return chrom, start, end

    def get_vcf_coords(self, transcript=None):
        """Return genomic coordinates of reference allele in VCF-style."""
        chrom, start, end = self.get_coords(transcript)

        # Inserts and deletes require left-padding by 1 base
        if self.mutation_type == ">":
            pass
        elif self.mutation_type in ("del", "ins", "dup", "delins"):
            # Indels have left-padding.
            start -= 1
        else:
            raise NotImplementedError("unknown mutation_type '%s'" %
                                      self.mutation_type)
        return chrom, start, end

    def get_ref_alt(self, is_forward_strand=True):
        """Return reference and alternate alleles."""
        if self.kind == 'p':
            raise NotImplementedError(
                'get_ref_alt is not implemented for protein HGVS names')
        alleles = [self.ref_allele, self.alt_allele]

        # Represent duplications are inserts.
        if self.mutation_type == "dup":
            alleles[0] = ""
            alleles[1] = alleles[1][:len(alleles[1]) / 2]

        if is_forward_strand:
            return alleles
        else:
            return tuple(map(revcomp, alleles))


def is_duplication(chrom, offset, ref, alt, genome):
    """
    Determines if allele is a duplication.

    chrom: Chromosome name.
    offset: 1-index genomic coordinate.
    ref: Reference allele (no padding).
    alt: Alternate allele (no padding).
    genome: pygr compatible genome object.

    Returns duplicated region [start, end] if allele is an insert that
    could be represented as a duplication. Otherwise, returns None.
    """

    if len(ref) == len(alt) == 0:
        # it's a SNP, just return.
        return None

    if len(ref) > 0 and len(alt) > 0:
        # complex indel, don't know how to dup check
        return None

    if len(ref) > len(alt):
        # deletion -- don't dup check
        return None

    indel_seq = alt
    indel_length = len(indel_seq)

    # Convert offset to 0-index.
    offset -= 1

    # Get genomic sequence around the lesion.
    prev_seq = unicode(
        genome[str(chrom)][offset - indel_length:offset]).upper()
    next_seq = unicode(
        genome[str(chrom)][offset:offset + indel_length]).upper()

    if prev_seq == indel_seq:
        return (offset - indel_length + 1, offset)
    if next_seq == indel_seq:
        return (offset + 1, offset + indel_length)
    else:
        return None


def hgvs_justify_indel(chrom, offset, ref, alt, strand, genome):
    """
    3' justify an indel according to the HGVS standard.

    Returns (offset, ref, alt).
    """
    if len(ref) == len(alt) == 0:
        # It's a SNP, just return.
        return offset, ref, alt

    if len(ref) > 0 and len(alt) > 0:
        # Complex indel, don't know how to justify.
        return offset, ref, alt

    # Get genomic sequence around the lesion.
    start = max(offset - 1000, 0)
    end = offset + 1000
    seq = unicode(genome[str(chrom)][start-1:end]).upper()
    cds_offset = offset - start

    # indel -- strip off the ref base to get the actual lesion sequence
    is_insert = len(alt) > 0
    if is_insert:
        indel_seq = alt
        cds_offset_end = cds_offset
    else:
        indel_seq = ref
        cds_offset_end = cds_offset + len(indel_seq)

    # now 3' justify (vs. cDNA not genome) the offset
    justify = 'right' if strand == '+' else 'left'
    offset, _, indel_seq = justify_indel(
        cds_offset, cds_offset_end, indel_seq, seq, justify)
    offset += start

    if is_insert:
        alt = indel_seq
    else:
        ref = indel_seq

    return offset, ref, alt


def parse_hgvs(
        hgvs_name, genome, transcript=None,
        get_transcript=lambda name: None,
        flank_length=30, normalize=True):
    """
    Parse an HGVS name into (chrom, start, end, ref, alt)

    hgvs_name: HGVS name to parse.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to HGVS name.
    normalize: If True, normalize allele according to VCF standard.
    """
    hgvs = HGVSName(hgvs_name)

    # Determine transcript.
    if hgvs.kind == 'c' and not transcript:
        if hgvs.transcript and get_transcript:
            transcript = get_transcript(hgvs.transcript)
        if not transcript:
            raise ValueError('transcript is required')

    chrom, start, end, ref, alt = get_vcf_allele(hgvs, genome, transcript)
    if normalize:
        chrom, start, ref, [alt] = normalize_variant(
            chrom, start, ref, [alt], genome,
            flank_length=flank_length).variant
    return (chrom, start, ref, alt)


def format_hgvs_name(chrom, offset, ref, alt, genome, transcript,
                     use_gene=True, max_allele_length=4):
    """
    Generate a HGVS name from a genomic coordinate.

    chrom: Chromosome name.
    offset: Genomic offset of allele.
    ref: Reference allele.
    alt: Alternate allele.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to allele.
    use_gene: Include transcript and gene in HGVS name.
    max_allele_length: If allele is greater than this use allele length.
    """

    def format_hgvs_name_prefix(transcript):
        return "%s{%s.%d}" % (
            transcript.gene.name, transcript.name, transcript.version)

    def format_coords(coords):
        if len(coords) == 1:
            return str(coords[0])
        elif len(coords) == 2:
            return str(coords[0]) + '_' + str(coords[1])
        else:
            raise ValueError('coords should be length 1 or 2')

    def format_allele(chrom, offset, ref, alt, strand, genome):

        # Find cDNA allele sequence.
        genomic_ref, genomic_alt = ref, alt
        if strand == "-":
            ref = revcomp(ref)
            alt = revcomp(alt)

        ref_len = len(ref)
        alt_len = len(alt)

        # Format allele.
        if len(ref) == 1 and len(alt) == 1:
            # SNP.
            mutation_type = '>'
            delta_seq = '%s>%s' % (ref, alt)
            delta_len = 1
        elif len(ref) > 0 and len(alt) == 0:
            # Deletion.
            mutation_type = 'del'
            if ref_len <= max_allele_length:
                delta_seq = "del" + ref
            else:
                delta_seq = "del" + str(ref_len)
            delta_len = ref_len
        elif len(ref) == 0 and len(alt) > 0:
            # Insertion.
            mutation_type = 'ins'
            if alt_len <= max_allele_length:
                delta_seq = "ins" + alt
            else:
                delta_seq = "ins" + str(alt_len)
            delta_len = alt_len

            # Check for possible duplication.
            dup_region = is_duplication(
                chrom, offset, genomic_ref, genomic_alt, genome)
            if dup_region:
                offset = dup_region[0]
                delta_seq = re.sub("ins", "dup", delta_seq)
                mutation_type = 'dup'
        else:
            # Indel.
            mutation_type = 'indel'
            if (ref_len > max_allele_length or alt_len > max_allele_length):
                delta_seq = "del" + str(ref_len) + "ins" + str(alt_len)
            else:
                delta_seq = "del" + ref + "ins" + alt
            delta_len = ref_len

        return mutation_type, offset, delta_seq, delta_len

    def compute_allele_coords(mutation_type, offset, delta_len, transcript):
        if delta_len == 1:
            cds_coords = [CDNACoord(
                *transcript.genomic_offset_to_cds_coord(offset))]
        else:
            if mutation_type == 'ins':
                # Insert uses coordinates around the insert point.
                offset_start = offset - 1
                offset_end = offset
            else:
                offset_start = offset
                offset_end = offset + delta_len - 1
            if transcript.strand == '-':
                offset_start, offset_end = offset_end, offset_start
            cds_coords = [
                CDNACoord(
                    *transcript.genomic_offset_to_cds_coord(offset_start)),
                CDNACoord(
                    *transcript.genomic_offset_to_cds_coord(offset_end))]

        return cds_coords

    def format_hgvs_name_snp(chrom, offset, ref, alt, genome, transcript):
        if transcript:
            coords = [CDNACoord(
                *transcript.genomic_offset_to_cds_coord(offset))]
            mutation_type, offset, delta_seq, delta_len = format_allele(
                chrom, offset, ref, alt, transcript.strand, genome)
            return 'c.' + format_coords(coords) + delta_seq

        else:
            return 'g.' + offset + format_allele(
                chrom, offset, ref, alt, '+', genome)

    def format_hgvs_name_indel(chrom, offset, ref, alt, genome, transcript):
        # remove 1bp padding
        offset += 1
        ref = ref[1:]
        alt = alt[1:]

        strand = transcript.strand if transcript else '+'

        # 3' justify allele.
        shifted_offset, ref, alt = hgvs_justify_indel(
            chrom, offset, ref, alt, strand, genome)

        # Format the allele part of the name.
        mutation_type, shifted_offset, delta_seq, delta_len = format_allele(
            chrom, shifted_offset, ref, alt, strand, genome)

        # If there is a transcript, compute the cDNA coordinate.
        if transcript:
            cds_coords = compute_allele_coords(
                mutation_type, shifted_offset, delta_len, transcript)
            c_str = 'c.' + format_coords(cds_coords)
            return c_str + delta_seq
        else:
            # No transcript in the DB for some reason, fall back to "g." coord.
            g_start = offset
            g_end = g_start + delta_len - 1
            g_str = "g.%d_%d" % (g_start, g_end)
            return g_str + delta_seq

    if len(ref) == len(alt) == 1:
        name = format_hgvs_name_snp(
            chrom, offset, ref, alt, genome, transcript)
    else:
        name = format_hgvs_name_indel(
            chrom, offset, ref, alt, genome, transcript)

    # Add gene and transcript prefix:
    if use_gene:
        prefix_name = format_hgvs_name_prefix(transcript)
        name = prefix_name + ':' + name
    return name
