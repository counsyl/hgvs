"""
Methods for manipulating HGVS names

Recommendations for the HGVS naming standard:
http://www.hgvs.org/mutnomen/standards.html

Definition of which transcript to use coding variants:
ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene


HGVS language currently implemented.

HGVS = ALLELE
     | PREFIX_NAME : ALLELE

PREFIX_NAME = TRANSCRIPT
            | TRANSCRIPT '(' GENE ')'

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
CDNA_ALLELE = CDNA_COORD SINGLE_BASE_CHANGE
            | CDNA_COORD_RANGE MULTI_BASE_CHANGE

GENOMIC_ALLELE =
MIT_ALLELE = COORD SINGLE_BASE_CHANGE
           | COORD_RANGE MULTI_BASE_CHANGE

SINGLE_BASE_CHANGE = CDNA_ALLELE = CDNA_COORD BASE '='        # no change
                   | CDNA_COORD BASE '>' BASE                 # substitution
                   | CDNA_COORD 'ins' BASE                    # 1bp insertion
                   | CDNA_COORD 'del' BASE                    # 1bp deletion
                   | CDNA_COORD 'dup' BASE                    # 1bp duplication
                   | CDNA_COORD 'ins'                         # 1bp insertion
                   | CDNA_COORD 'del'                         # 1bp deletion
                   | CDNA_COORD 'dup'                         # 1bp duplication
                   | CDNA_COORD 'del' BASE 'ins' BASE         # 1bp indel
                   | CDNA_COORD 'delins' BASE                 # 1bp indel

MULTI_BASE_CHANGE = COORD_RANGE 'del' BASES             # deletion
                  | COORD_RANGE 'ins' BASES             # insertion
                  | COORD_RANGE 'dup' BASES             # duplication
                  | COORD_RANGE 'del'                   # deletion
                  | COORD_RANGE 'dup'                   # duplication
                  | COORD_RANGE 'del' BASES 'ins' BASES # indel
                  | COORD_RANGE 'delins' BASES          # indel


AMINO1 = [GAVLIMFWPSTCYNQDEKRH]

AMINO3 = 'Gly' | 'Ala' | 'Val' | 'Leu' | 'Ile' | 'Met' | 'Phe' | 'Trp' | 'Pro'
       | 'Ser' | 'Thr' | 'Cys' | 'Tyr' | 'Asn' | 'Gln' | 'Asp' | 'Glu' | 'Lys'
       | 'Arg' | 'His'

PROTEIN_ALLELE = AMINO3 COORD '='               # no peptide change
               | AMINO1 COORD '='               # no peptide change
               | AMINO3 COORD AMINO3 PEP_EXTRA  # peptide change
               | AMINO1 COORD AMINO1 PEP_EXTRA  # peptide change
               | AMINO3 COORD '_' AMINO3 COORD PEP_EXTRA        # indel
               | AMINO1 COORD '_' AMINO1 COORD PEP_EXTRA        # indel
               | AMINO3 COORD '_' AMINO3 COORD PEP_EXTRA AMINO3 # indel
               | AMINO1 COORD '_' AMINO1 COORD PEP_EXTRA AMINO1 # indel

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
from __future__ import absolute_import
from __future__ import unicode_literals

import re

from .variants import justify_indel
from .variants import normalize_variant
from .variants import revcomp


CHROM_PREFIX = 'chr'
CDNA_START_CODON = 'cdna_start'
CDNA_STOP_CODON = 'cdna_stop'


class HGVSRegex(object):
    """
    All regular expression for HGVS names.
    """

    # DNA syntax
    # http://www.hgvs.org/mutnomen/standards.html#nucleotide
    BASE = "[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]|\d+"
    BASES = "[acgtbdhkmnrsvwyACGTBDHKMNRSVWY]+|\d+"
    DNA_REF = "(?P<ref>" + BASES + ")"
    DNA_ALT = "(?P<alt>" + BASES + ")"

    # Mutation types
    EQUAL = "(?P<mutation_type>=)"
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
        # No change
        CDNA_START + DNA_REF + EQUAL,

        # Substitution
        CDNA_START + DNA_REF + SUB + DNA_ALT,

        # 1bp insertion, deletion, duplication
        CDNA_START + INS + DNA_ALT,
        CDNA_START + DEL + DNA_REF,
        CDNA_START + DUP + DNA_REF,
        CDNA_START + DEL,
        CDNA_START + DUP,

        # Insertion, deletion, duplication
        CDNA_RANGE + INS + DNA_ALT,
        CDNA_RANGE + DEL + DNA_REF,
        CDNA_RANGE + DUP + DNA_REF,
        CDNA_RANGE + DEL,
        CDNA_RANGE + DUP,

        # Indels
        "(?P<delins>" + CDNA_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + CDNA_RANGE + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + CDNA_START + 'delins' + DNA_ALT + ")",
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
        "(?P<delins>" + PEP_REF + COORD_START + "_" + PEP_REF2 + COORD_END +
        PEP_EXTRA + ")",
        "(?P<delins>" + PEP_REF + COORD_START + "_" + PEP_REF2 + COORD_END +
        PEP_ALT + PEP_EXTRA + ")",
    ]

    PEP_ALLELE_REGEXES = [re.compile("^" + regex + "$")
                          for regex in PEP_ALLELE]

    # Genomic allele syntax
    GENOMIC_ALLELE = [
        # No change
        COORD_START + DNA_REF + EQUAL,

        # Substitution
        COORD_START + DNA_REF + SUB + DNA_ALT,

        # 1bp insertion, deletion, duplication
        COORD_START + INS + DNA_ALT,
        COORD_START + DEL + DNA_REF,
        COORD_START + DUP + DNA_REF,
        COORD_START + DEL,
        COORD_START + DUP,

        # Insertion, deletion, duplication
        COORD_RANGE + INS + DNA_ALT,
        COORD_RANGE + DEL + DNA_REF,
        COORD_RANGE + DUP + DNA_REF,
        COORD_RANGE + DEL,
        COORD_RANGE + DUP,

        # Indels
        "(?P<delins>" + COORD_START + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + COORD_RANGE + 'del' + DNA_REF + 'ins' + DNA_ALT + ")",
        "(?P<delins>" + COORD_START + 'delins' + DNA_ALT + ")",
        "(?P<delins>" + COORD_RANGE + 'delins' + DNA_ALT + ")",
    ]

    GENOMIC_ALLELE_REGEXES = [re.compile("^" + regex + "$")
                              for regex in GENOMIC_ALLELE]


class ChromosomeSubset(object):
    """
    Allow direct access to a subset of the chromosome.
    """

    def __init__(self, name, genome=None):
        self.name = name
        self.genome = genome

    def __getitem__(self, key):
        """Return sequence from region [start, end)

        Coordinates are 0-based, end-exclusive."""
        if isinstance(key, slice):
            start, end = (key.start, key.stop)
            start -= self.genome.start
            end -= self.genome.start
            return self.genome.genome[self.genome.seqid][start:end]
        else:
            raise TypeError('Expected a slice object but '
                            'received a {0}.'.format(type(key)))

    def __repr__(self):
        return 'ChromosomeSubset("%s")' % self.name


class GenomeSubset(object):
    """
    Allow the direct access of a subset of the genome.
    """

    def __init__(self, genome, chrom, start, end, seqid):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.end = end
        self.seqid = seqid
        self._chroms = {}

    def __getitem__(self, chrom):
        """Return a chromosome by its name."""
        if chrom in self._chroms:
            return self._chroms[chrom]
        else:
            chromosome = ChromosomeSubset(chrom, self)
            self._chroms[chrom] = chromosome
            return chromosome


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


# The RefSeq standard for naming contigs/transcripts/proteins:
# http://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly  # nopep8
REFSEQ_PREFIXES = [
    ('AC_', 'genomic',
     'Complete genomic molecule, usually alternate assembly'),
    ('NC_', 'genomic',
     'Complete genomic molecule, usually reference assembly'),
    ('NG_', 'genomic', 'Incomplete genomic region'),
    ('NT_', 'genomic', 'Contig or scaffold, clone-based or WGS'),
    ('NW_', 'genomic', 'Contig or scaffold, primarily WGS'),
    ('NS_', 'genomic', 'Environmental sequence'),
    ('NZ_', 'genomic', 'Unfinished WGS'),
    ('NM_', 'mRNA', ''),
    ('NR_', 'RNA', ''),
    ('XM_', 'mRNA', 'Predicted model'),
    ('XR_', 'RNA', 'Predicted model'),
    ('AP_', 'Protein', 'Annotated on AC_ alternate assembly'),
    ('NP_', 'Protein', 'Associated with an NM_ or NC_ accession'),
    ('YP_', 'Protein', ''),
    ('XP_', 'Protein', 'Predicted model, associated with an XM_ accession'),
    ('ZP_', 'Protein', 'Predicted model, annotated on NZ_ genomic records'),
]

REFSEQ_PREFIX_LOOKUP = dict(
    (prefix, (kind, description))
    for prefix, kind, description in REFSEQ_PREFIXES
)


def get_refseq_type(name):
    """
    Return the RefSeq type for a refseq name.
    """
    prefix = name[:3]
    return REFSEQ_PREFIX_LOOKUP.get(prefix, (None, ''))[0]


def get_exons(transcript):
    """Yield exons in coding order."""
    transcript_strand = transcript.tx_position.is_forward_strand
    if hasattr(transcript.exons, 'select_related'):
        exons = list(transcript.exons.select_related('tx_position'))
    else:
        exons = list(transcript.exons)
    exons.sort(key=lambda exon: exon.tx_position.chrom_start)
    if not transcript_strand:
        exons.reverse()
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
                   else transcript.cds_position.chrom_stop - 1)
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
        return cdna_len + (exon_end - start_codon - 1)


def find_stop_codon(exons, cds_position):
    """Return the position along the cDNA of the base after the stop codon."""
    if cds_position.is_forward_strand:
        stop_pos = cds_position.chrom_stop
    else:
        stop_pos = cds_position.chrom_start
    cdna_pos = 0
    for exon in exons:
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


def get_genomic_sequence(genome, chrom, start, end):
    """
    Return a sequence for the genomic region.

    start, end: 1-based, end-inclusive coordinates of the sequence.
    """
    if start > end:
        return ''
    else:
        return str(genome[str(chrom)][start - 1:end]).upper()


def cdna_to_genomic_coord(transcript, coord):
    """Convert a HGVS cDNA coordinate to a genomic coordinate."""
    transcript_strand = transcript.tx_position.is_forward_strand
    exons = get_exons(transcript)
    utr5p = (get_utr5p_size(transcript)
             if transcript.is_coding else 0)

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
        pos = find_stop_codon(exons, transcript.cds_position) + coord.coord
    else:
        raise ValueError('unknown CDNACoord landmark "%s"' % coord.landmark)

    # 5' flanking sequence.
    if pos < 1:
        if transcript_strand:
            return transcript.tx_position.chrom_start + pos
        else:
            return transcript.tx_position.chrom_stop - pos + 1

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
        # 3' flanking sequence
        if transcript_strand:
            return transcript.cds_position.chrom_stop + coord.coord
        else:
            return transcript.cds_position.chrom_start + 1 - coord.coord

    # Compute genomic coordinate using offset.
    if transcript_strand:
        # Plus strand.
        return exon_start + (pos - cdna_start) + coord.offset
    else:
        # Minus strand.
        return exon_end - (pos - cdna_start) - coord.offset


def genomic_to_cdna_coord(transcript, genomic_coord):
    """Convert a genomic coordinate to a cDNA coordinate and offset.
    """
    exons = [exon.get_as_interval()
             for exon in get_exons(transcript)]

    if len(exons) == 0:
        return None

    strand = transcript.strand

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
                if (genomic_coord < transcript.tx_position.chrom_start + 1 or
                        genomic_coord > transcript.tx_position.chrom_stop):
                    nearest_exonic += distance
                    distance = 0
                cdna_coord = CDNACoord(nearest_exonic, distance)
            break
        coding_offset += exon_length

    # Adjust coordinates for coding transcript.
    if transcript.is_coding:
        # Detect if position before start codon.
        utr5p = get_utr5p_size(transcript) if transcript.is_coding else 0
        cdna_coord.coord -= utr5p
        if cdna_coord.coord <= 0:
            cdna_coord.coord -= 1
        else:
            # Detect if position is after stop_codon.
            exons = get_exons(transcript)
            stop_codon = find_stop_codon(exons, transcript.cds_position)
            stop_codon -= utr5p
            if (cdna_coord.coord > stop_codon or
                    cdna_coord.coord == stop_codon and cdna_coord.offset > 0):
                cdna_coord.coord -= stop_codon
                cdna_coord.landmark = CDNA_STOP_CODON

    return cdna_coord


def get_allele(hgvs, genome, transcript=None):
    """Get an allele from a HGVSName, a genome, and a transcript."""
    chrom, start, end = hgvs.get_coords(transcript)
    _, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True)
    ref = get_genomic_sequence(genome, chrom, start, end)
    return chrom, start, end, ref, alt


_indel_mutation_types = set(['ins', 'del', 'dup', 'delins'])


def get_vcf_allele(hgvs, genome, transcript=None):
    """Get an VCF-style allele from a HGVSName, a genome, and a transcript."""
    chrom, start, end = hgvs.get_vcf_coords(transcript)
    _, alt = hgvs.get_ref_alt(
        transcript.tx_position.is_forward_strand if transcript else True)
    ref = get_genomic_sequence(genome, chrom, start, end)

    if hgvs.mutation_type in _indel_mutation_types:
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
    def __init__(self, name='', part='name', reason=''):
        if name:
            message = 'Invalid HGVS %s "%s"' % (part, name)
        else:
            message = 'Invalid HGVS %s' % part
        if reason:
            message += ': ' + reason
        super(InvalidHGVSName, self).__init__(message)

        self.name = name
        self.part = part
        self.reason = reason


class HGVSName(object):
    """
    Represents a HGVS variant name.
    """

    def __init__(self, name='', prefix='', chrom='', transcript='', gene='',
                 kind='', mutation_type=None, start=0, end=0, ref_allele='',
                 ref2_allele='', alt_allele='',
                 cdna_start=None, cdna_end=None, pep_extra=''):

        # Full HGVS name.
        self.name = name

        # Name parts.
        self.prefix = prefix
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
        self.cdna_end = cdna_end if cdna_end else CDNACoord()

        # Protein-specific fields
        self.pep_extra = pep_extra

        if name:
            self.parse(name)

    def parse(self, name):
        """Parse a HGVS name."""
        # Does HGVS name have transcript/gene prefix?
        if ':' in name:
            prefix, allele = name.split(':', 1)
        else:
            prefix = ''
            allele = name

        self.name = name

        # Parse prefix and allele.
        self.parse_allele(allele)
        self.parse_prefix(prefix, self.kind)
        self._validate()

    def parse_prefix(self, prefix, kind):
        """
        Parse a HGVS prefix (gene/transcript/chromosome).

        Some examples of full hgvs names with transcript include:
          NM_007294.3:c.2207A>C
          NM_007294.3(BRCA1):c.2207A>C
          BRCA1{NM_007294.3}:c.2207A>C
        """

        self.prefix = prefix

        # No prefix.
        if prefix == '':
            self.chrom = ''
            self.transcript = ''
            self.gene = ''
            return

        # Transcript and gene given with parens.
        # example: NM_007294.3(BRCA1):c.2207A>C
        match = re.match("^(?P<transcript>[^(]+)\((?P<gene>[^)]+)\)$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = match.group('gene')
            return

        # Transcript and gene given with braces.
        # example: BRCA1{NM_007294.3}:c.2207A>C
        match = re.match("^(?P<gene>[^{]+)\{(?P<transcript>[^}]+)\}$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = match.group('gene')
            return

        # Determine using Ensembl type.
        if prefix.startswith('ENST'):
            self.transcript = prefix
            return

        # Determine using refseq type.
        refseq_type = get_refseq_type(prefix)
        if refseq_type in ('mRNA', 'RNA'):
            self.transcript = prefix
            return

        # Determine using refseq type.
        if prefix.startswith(CHROM_PREFIX) or refseq_type == 'genomic':
            self.chrom = prefix
            return

        # Assume gene name.
        self.gene = prefix

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
            raise InvalidHGVSName(allele, 'allele',
                                  'expected kind "c.", "p.", "g.", etc')

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
            raise NotImplementedError("unknown kind: %s" % allele)

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

                # Convert no match alleles.
                if self.mutation_type == "=":
                    self.alt_allele = self.ref_allele
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

                # Convert no match alleles.
                if self.mutation_type == "=":
                    self.alt_allele = self.ref_allele
                return

        raise InvalidHGVSName(details, 'genomic allele')

    def _validate(self):
        """
        Check for internal inconsistencies in representation
        """
        if self.start > self.end:
            raise InvalidHGVSName(reason="Coordinates are nonincreasing")

    def __repr__(self):
        try:
            return "HGVSName('%s')" % self.format()
        except NotImplementedError:
            return "HGVSName('%s')" % self.name

    def __unicode__(self):
        return self.format()

    def format(self, use_prefix=True, use_gene=True, use_counsyl=False):
        """Generate a HGVS name as a string."""

        if self.kind == 'c':
            allele = 'c.' + self.format_cdna()
        elif self.kind == 'p':
            allele = 'p.' + self.format_protein()
        elif self.kind == 'g':
            allele = 'g.' + self.format_genome()
        else:
            raise NotImplementedError("not implemented: '%s'" % self.kind)

        prefix = self.format_prefix(use_gene=use_gene) if use_prefix else ''

        if prefix:
            return prefix + ':' + allele
        else:
            return allele

    def format_prefix(self, use_gene=True):
        """
        Generate HGVS trancript/gene prefix.

        Some examples of full hgvs names with transcript include:
          NM_007294.3:c.2207A>C
          NM_007294.3(BRCA1):c.2207A>C
        """

        if self.kind == 'g':
            if self.chrom:
                return self.chrom

        if self.transcript:
            if use_gene and self.gene:
                return '%s(%s)' % (self.transcript, self.gene)
            else:
                return self.transcript
        else:
            if use_gene:
                return self.gene
            else:
                return ''

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
        if self.mutation_type == '=':
            # No change.
            # example: 101A=
            return self.ref_allele + '='

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

        elif self.mutation_type == 'inv':
            return self.mutation_type

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
                self.ref_allele == self.ref2_allele == self.alt_allele):
            # Match.
            # Example: Glu1161=
            pep_extra = self.pep_extra if self.pep_extra else '='
            return self.ref_allele + str(self.start) + pep_extra

        elif (self.start == self.end and
                self.ref_allele == self.ref2_allele and
                self.ref_allele != self.alt_allele):
            # Change.
            # Example: Glu1161Ser
            return (self.ref_allele + str(self.start) +
                    self.alt_allele + self.pep_extra)

        elif self.start != self.end:
            # Range change.
            # Example: Glu1161_Ser1164?fs
            return (self.ref_allele + str(self.start) + '_' +
                    self.ref2_allele + str(self.end) +
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
            start = cdna_to_genomic_coord(transcript, self.cdna_start)
            end = cdna_to_genomic_coord(transcript, self.cdna_end)

            if not transcript.tx_position.is_forward_strand:
                if end > start:
                    raise AssertionError(
                        "cdna_start cannot be greater than cdna_end")
                start, end = end, start
            else:
                if start > end:
                    raise AssertionError(
                        "cdna_start cannot be greater than cdna_end")

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
                'Coordinates are not available for this kind of HGVS name "%s"'
                % self.kind)

        return chrom, start, end

    def get_vcf_coords(self, transcript=None):
        """Return genomic coordinates of reference allele in VCF-style."""
        chrom, start, end = self.get_coords(transcript)

        # Inserts and deletes require left-padding by 1 base
        if self.mutation_type in ("=", ">"):
            pass
        elif self.mutation_type in ("del", "ins", "dup", "delins"):
            # Indels have left-padding.
            start -= 1
        else:
            raise NotImplementedError("Unknown mutation_type '%s'" %
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
            alleles[1] = alleles[1][:len(alleles[1]) // 2]

        if is_forward_strand:
            return alleles
        else:
            return tuple(map(revcomp, alleles))


def hgvs_justify_dup(chrom, offset, ref, alt, genome):
    """
    Determines if allele is a duplication and justifies.

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
        return chrom, offset, ref, alt, '>'

    if len(ref) > 0 and len(alt) > 0:
        # complex indel, don't know how to dup check
        return chrom, offset, ref, alt, 'delins'

    if len(ref) > len(alt):
        # deletion -- don't dup check
        return chrom, offset, ref, alt, 'del'

    indel_seq = alt
    indel_length = len(indel_seq)

    # Convert offset to 0-index.
    offset -= 1

    # Get genomic sequence around the lesion.
    prev_seq = str(
        genome[str(chrom)][offset - indel_length:offset]).upper()
    next_seq = str(
        genome[str(chrom)][offset:offset + indel_length]).upper()

    # Convert offset back to 1-index.
    offset += 1

    if prev_seq == indel_seq:
        offset = offset - indel_length
        mutation_type = 'dup'
        ref = indel_seq
        alt = indel_seq * 2
    elif next_seq == indel_seq:
        mutation_type = 'dup'
        ref = indel_seq
        alt = indel_seq * 2
    else:
        mutation_type = 'ins'

    return chrom, offset, ref, alt, mutation_type


def hgvs_justify_indel(chrom, offset, ref, alt, strand, genome):
    """
    3' justify an indel according to the HGVS standard.

    Returns (chrom, offset, ref, alt).
    """
    if len(ref) == len(alt) == 0:
        # It's a SNP, just return.
        return chrom, offset, ref, alt

    if len(ref) > 0 and len(alt) > 0:
        # Complex indel, don't know how to justify.
        return chrom, offset, ref, alt

    # Get genomic sequence around the lesion.
    start = max(offset - 100, 0)
    end = offset + 100
    seq = str(genome[str(chrom)][start - 1:end]).upper()
    cds_offset = offset - start

    # indel -- strip off the ref base to get the actual lesion sequence
    is_insert = len(alt) > 0
    if is_insert:
        indel_seq = alt
        cds_offset_end = cds_offset
    else:
        indel_seq = ref
        cds_offset_end = cds_offset + len(indel_seq)

    # Now 3' justify (vs. cDNA not genome) the offset
    justify = 'right' if strand == '+' else 'left'
    offset, _, indel_seq = justify_indel(
        cds_offset, cds_offset_end, indel_seq, seq, justify)
    offset += start

    if is_insert:
        alt = indel_seq
    else:
        ref = indel_seq

    return chrom, offset, ref, alt


def hgvs_normalize_variant(chrom, offset, ref, alt, genome, transcript=None):
    """Convert VCF-style variant to HGVS-style."""
    if len(ref) == len(alt) == 1:
        if ref == alt:
            mutation_type = '='
        else:
            mutation_type = '>'
    else:
        # Remove 1bp padding
        offset += 1
        ref = ref[1:]
        alt = alt[1:]

        # 3' justify allele.
        strand = transcript.strand if transcript else '+'
        chrom, offset, ref, alt = hgvs_justify_indel(
            chrom, offset, ref, alt, strand, genome)

        # Represent as duplication if possible.
        chrom, offset, ref, alt, mutation_type = hgvs_justify_dup(
            chrom, offset, ref, alt, genome)
    return chrom, offset, ref, alt, mutation_type


def parse_hgvs_name(hgvs_name, genome, transcript=None,
                    get_transcript=lambda name: None,
                    flank_length=30, normalize=True, lazy=False):
    """
    Parse an HGVS name into (chrom, start, end, ref, alt)

    hgvs_name: HGVS name to parse.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to HGVS name.
    normalize: If True, normalize allele according to VCF standard.
    lazy: If True, discard version information from incoming transcript/gene.
    """
    hgvs = HGVSName(hgvs_name)

    # Determine transcript.
    if hgvs.kind == 'c' and not transcript:
        if '.' in hgvs.transcript and lazy:
            hgvs.transcript, version = hgvs.transcript.split('.')
        elif '.' in hgvs.gene and lazy:
            hgvs.gene, version = hgvs.gene.split('.')
        if get_transcript:
            if hgvs.transcript:
                transcript = get_transcript(hgvs.transcript)
            elif hgvs.gene:
                transcript = get_transcript(hgvs.gene)
        if not transcript:
            raise ValueError('transcript is required')

    if transcript and hgvs.transcript in genome:
        # Reference sequence is directly known, use it.
        genome = GenomeSubset(genome, transcript.tx_position.chrom,
                              transcript.tx_position.chrom_start,
                              transcript.tx_position.chrom_stop,
                              hgvs.transcript)

    chrom, start, end, ref, alt = get_vcf_allele(hgvs, genome, transcript)
    if normalize:
        chrom, start, ref, [alt] = normalize_variant(
            chrom, start, ref, [alt], genome,
            flank_length=flank_length).variant
    return (chrom, start, ref, alt)


def variant_to_hgvs_name(chrom, offset, ref, alt, genome, transcript,
                         max_allele_length=4, use_counsyl=False):
    """
    Populate a HGVSName from a genomic coordinate.

    chrom: Chromosome name.
    offset: Genomic offset of allele.
    ref: Reference allele.
    alt: Alternate allele.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to allele.
    max_allele_length: If allele is greater than this use allele length.
    """
    # Convert VCF-style variant to HGVS-style.
    chrom, offset, ref, [alt] = normalize_variant(
        chrom, offset, ref, [alt], genome).variant
    chrom, offset, ref, alt, mutation_type = hgvs_normalize_variant(
        chrom, offset, ref, alt, genome, transcript)

    # Populate HGVSName parse tree.
    hgvs = HGVSName()

    # Populate coordinates.
    if not transcript:
        # Use genomic coordinate when no transcript is available.
        hgvs.kind = 'g'
        hgvs.start = offset
        hgvs.end = offset + len(ref) - 1
    else:
        # Use cDNA coordinates.
        hgvs.kind = 'c'
        is_single_base_indel = (
            (mutation_type == 'ins' and len(alt) == 1) or
            (mutation_type in ('del', 'delins', 'dup') and len(ref) == 1))

        if mutation_type == '>' or (use_counsyl and is_single_base_indel):
            # Use a single coordinate.
            hgvs.cdna_start = genomic_to_cdna_coord(transcript, offset)
            hgvs.cdna_end = hgvs.cdna_start
        else:
            # Use a range of coordinates.
            if mutation_type == 'ins':
                # Insert uses coordinates around the insert point.
                offset_start = offset - 1
                offset_end = offset
            else:
                offset_start = offset
                offset_end = offset + len(ref) - 1
            if transcript.strand == '-':
                offset_start, offset_end = offset_end, offset_start
            hgvs.cdna_start = genomic_to_cdna_coord(transcript, offset_start)
            hgvs.cdna_end = genomic_to_cdna_coord(transcript, offset_end)

    # Populate prefix.
    if transcript:
        hgvs.transcript = transcript.full_name
        hgvs.gene = transcript.gene.name

    # Convert alleles to transcript strand.
    if transcript and transcript.strand == '-':
        ref = revcomp(ref)
        alt = revcomp(alt)

    # Convert to allele length if alleles are long.
    ref_len = len(ref)
    alt_len = len(alt)
    if ((mutation_type == 'dup' and ref_len > max_allele_length) or
            (mutation_type != 'dup' and
             (ref_len > max_allele_length or alt_len > max_allele_length))):
        ref = str(ref_len)
        alt = str(alt_len)

    # Populate alleles.
    hgvs.mutation_type = mutation_type
    hgvs.ref_allele = ref
    hgvs.alt_allele = alt

    return hgvs


def format_hgvs_name(chrom, offset, ref, alt, genome, transcript,
                     use_prefix=True, use_gene=True, use_counsyl=False,
                     max_allele_length=4):
    """
    Generate a HGVS name from a genomic coordinate.

    chrom: Chromosome name.
    offset: Genomic offset of allele.
    ref: Reference allele.
    alt: Alternate allele.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to allele.
    use_prefix: Include a transcript/gene/chromosome prefix in HGVS name.
    use_gene: Include gene name in HGVS prefix.
    max_allele_length: If allele is greater than this use allele length.
    """
    hgvs = variant_to_hgvs_name(chrom, offset, ref, alt, genome, transcript,
                                max_allele_length=max_allele_length,
                                use_counsyl=use_counsyl)
    return hgvs.format(use_prefix=use_prefix, use_gene=use_gene,
                       use_counsyl=use_counsyl)
