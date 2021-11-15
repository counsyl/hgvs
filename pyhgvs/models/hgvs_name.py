"""
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

import re

from .cdna import CDNACoord
from .variants import revcomp

CHROM_PREFIX = 'chr'


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
        self.ref_allele = ref_allele  # reference allele
        self.ref2_allele = ref2_allele  # reference allele at end of pep indel
        self.alt_allele = alt_allele  # alternate allele

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
        match = re.match(r"^(?P<transcript>[^(]+)\((?P<gene>[^)]+)\)$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = match.group('gene')
            return

        # Transcript and gene given with braces.
        # example: BRCA1{NM_007294.3}:c.2207A>C
        match = re.match(r"^(?P<gene>[^{]+){(?P<transcript>[^}]+)}$", prefix)
        if match:
            self.transcript = match.group('transcript')
            self.gene = match.group('gene')
            return

        # Determine using Ensembl type.
        if prefix.startswith('ENST'):
            self.transcript = prefix
            return

        # Determine using LRG type.
        if prefix.startswith('LRG_'):
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

        if kind in ("c", 'n'):
            self.parse_cdna(details)
            if kind == 'n':  # Ensure no 3'UTR or 5'UTR coords in non-coding
                if self.cdna_start.coord < 0:
                    raise InvalidHGVSName(allele, 'allele',
                                          "Non-coding transcript cannot contain negative (5'UTR) coordinates")
                if self.cdna_start.landmark == 'cdna_stop' or self.cdna_end and self.cdna_end.landmark == 'cdna_stop':
                    raise InvalidHGVSName(allele, 'allele',
                                          "Non-coding transcript cannot contain '*' (3'UTR) coordinates")
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
        elif self.kind == 'n':
            allele = 'n.' + self.format_cdna()
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

    def get_raw_coords(self, transcript=None):
        """ return genomic coordinates """
        if self.kind in ('c', 'n'):
            chrom = transcript.tx_position.chrom
            start = transcript.cdna_to_genomic_coord(self.cdna_start)
            end = transcript.cdna_to_genomic_coord(self.cdna_end)

            if not transcript.tx_position.is_forward_strand:
                if end > start:
                    raise AssertionError(
                        "cdna_start cannot be greater than cdna_end")
                start, end = end, start

            if start > end:
                raise AssertionError(
                    "cdna_start cannot be greater than cdna_end")
        elif self.kind == 'g':
            chrom = self.chrom
            start = self.start
            end = self.end
        else:
            raise NotImplementedError(
                'Coordinates are not available for this kind of HGVS name "%s"'
                % self.kind)

        return chrom, start, end

    def get_ref_coords(self, transcript=None):
        """Return genomic coordinates of reference allele."""

        chrom, start, end = self.get_raw_coords(transcript)

        if self.mutation_type == "ins":
            # Inserts have empty interval.
            if start < end:
                start += 1
                end -= 1
            else:
                end = start - 1

        elif self.mutation_type == "dup":
            end = start - 1
        return chrom, start, end

    def get_vcf_coords(self, transcript=None):
        """Return genomic coordinates of reference allele in VCF-style."""
        chrom, start, end = self.get_ref_coords(transcript)

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
