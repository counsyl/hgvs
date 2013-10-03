"""

Utitlities for manipulating variants.

"""


_COMP = dict(A='T', C='G', G='C', T='A', N='N',
             a='t', c='g', g='c', t='a', n='n')


def revcomp(seq):
    """Reverse complement."""
    return ''.join(_COMP[base] for base in reversed(seq))


def get_sequence(genome, chrom, start, end, is_forward_strand=True):
    """Return a sequence for the genomic region.

    Coordinates are 0-based, end-exclusive.
    """
    if start >= end:
        return ''
    else:
        seq = genome[str(chrom)][start:end]
        if not is_forward_strand:
            seq = -seq
        return str(seq).upper()


def get_sequence_from_position(genome, position):
    """Return a sequence for the genomic region

    Position is 0-based, end-exclusive.
    """
    return get_sequence(genome, position.chrom,
                        position.chrom_start, position.chrom_stop,
                        position.is_forward_strand)


def justify_indel(start, end, indel, seq, justify):
    """
    Justify an indel to the left or right along a sequence 'seq'.

    start, end: 0-based, end-exclusive coordinates of 'indel' within the
        sequence 'seq'. Inserts denote the insertion point using start=end
        and deletions indicate the deleted region with (start,end).
    indel: indel sequence, can be insertion or deletion.
    seq: a larger sequence containing the indel. Can be a fragment from the
        genome.
    justify: Which direction to justify the indel ('left', 'right').
    """
    if justify == 'left':
        while seq[start - 1] == indel[-1]:
            seq_added = seq[start - 1]
            indel = seq_added + indel[:-1]
            start -= 1
            end -= 1
    elif justify == 'right':
        while seq[end] == indel[0]:
            seq_added = seq[end]
            indel = indel[1:] + seq_added
            start += 1
            end += 1
    else:
        raise ValueError('unknown justify "%s"' % justify)
    return start, end, indel


def normalize_variant(chrom, offset, ref_sequence, alt_sequences, genome,
                      flank_length=30):
    """
    Normalize variant according to the GATK/VCF standard.

    chrom: chromsome containing variant.
    offset: 1-based coordinate of reference allele in the genome.
    ref_sequence: reference allele.
    alt_sequences: list of all alternate sequences.
    genome: pygr-compatiable genome object.
    """
    start = offset - 1
    end = offset + len(ref_sequence)
    upstream_flank = get_sequence(
        genome, chrom, start - flank_length, start)
    position = Position(
        chrom=chrom,
        chrom_start=start,
        chrom_stop=end,
        is_forward_strand=True)
    return NormalizedVariant(position, ref_sequence, alt_sequences,
                             seq_5p=upstream_flank)


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


class NormalizedVariant(object):
    """
    Normalizes variant representation to match GATK/VCF.
    """

    def __init__(self, position, ref_allele, alt_alleles,
                 seq_5p='', seq_3p=''):
        """
        position: a 0-index genomic Position.
        ref_allele: the reference allele sequence.
        alt_alleles: a list of alternate allele sequences.
        seq_5p: 5 prime flanking sequence of variant.
        seq_3p: 3 prime flanking sequence of variant.
        """
        self.position = position
        self.alleles = [ref_allele] + list(alt_alleles)
        self.seq_5p = seq_5p
        self.seq_3p = seq_3p
        self.log = []

        self._on_forward_strand()
        self._trim_common_prefix()
        self._trim_common_suffix()
        self._left_align()
        self._1bp_pad()
        self._set_1based_position()

    def _on_forward_strand(self):
        """
        Ensure variant is on forward strand.
        """
        if not self.position.is_forward_strand:
            self.log.append('flip strand')
            seq_5p = self.seq_5p
            seq_3p = self.seq_3p
            self.seq_5p = revcomp(seq_3p)
            self.seq_3p = revcomp(seq_5p)
            self.alleles = map(revcomp, self.alleles)

    def _trim_common_prefix(self):
        """
        Trim the common prefix amongst all alleles.
        """
        minlength = min(map(len, self.alleles))
        common_prefix = 0
        for i in range(minlength):
            if len(set(allele[i] for allele in self.alleles)) > 1:
                # Not all alleles match at this site, so common prefix ends.
                break
            common_prefix = i + 1

        # Remove common prefix from all alleles.
        if common_prefix:
            self.log.append('trim common prefix')
            self.position.chrom_start += common_prefix
            self.seq_5p += self.alleles[0][:common_prefix]
            for i, allele in enumerate(self.alleles):
                self.alleles[i] = allele[common_prefix:]

    def _trim_common_suffix(self):
        """
        Trim the common suffix amongst all alleles.
        """
        minlength = min(map(len, self.alleles))
        common_suffix = 0
        for i in range(1, minlength+1):
            if len(set(allele[-i] for allele in self.alleles)) > 1:
                # Not all alleles match at this site, so common suffix ends.
                break
            common_suffix = i

        # Remove common prefix from all alleles.
        if common_suffix:
            self.log.append('trim common suffix')
            self.seq_3p = self.alleles[0][-common_suffix:] + self.seq_3p
            for i, allele in enumerate(self.alleles):
                self.alleles[i] = allele[:-common_suffix]

    def _left_align(self):
        """
        Align variant as far to the left as possible.
        """
        # Left-aligning only makes sense for INDELs.
        if self.molecular_class != "INDEL":
            return

        # Identify the inserted or deleted sequence.
        alleles_with_seq = [i for i, allele in enumerate(self.alleles)
                            if allele]

        # Can only left-align biallelic, non ins-plus-del indels.
        if len(alleles_with_seq) == 1:
            i = alleles_with_seq[0]
            allele = self.alleles[i]
            offset = len(self.seq_5p)
            offset2, _, allele = justify_indel(
                offset, offset, allele, self.seq_5p, 'left')
            delta = offset - offset2
            if delta > 0:
                self.position.chrom_start -= delta
                self.position.chrom_stop -= delta
                self.seq_5p = self.seq_5p[:-delta]
                seq = self.ref_allele + self.seq_3p
                self.seq_3p = seq[:delta] + self.seq_3p
                self.alleles[i] = allele

    def _1bp_pad(self):
        """
        Ensure no alleles are the empty string by padding to the left 1bp.
        """
        # Padding is only required for INDELs.
        if self.molecular_class != "INDEL":
            return

        # Pad sequences with one 5-prime base before the mutation event.
        empty_seq = any(1 for allele in self.alleles if not allele)
        uniq_starts = set(allele[0] for allele in self.alleles if allele)
        if empty_seq or len(uniq_starts) > 1:
            self.log.append('1bp pad')
            for i, allele in enumerate(self.alleles):
                self.alleles[i] = self.seq_5p[-1] + self.alleles[i]

            self.seq_5p = self.seq_5p[:-1]
            self.position.chrom_start -= 1

        if len(set(a[0] for a in self.alleles)) != 1:
            raise AssertionError(
                "All INDEL alleles should start with same base.")

    def _set_1based_position(self):
        """
        Convert to 1-based end-inclusive coordinates.
        """
        self.position.chrom_start += 1

    @property
    def molecular_class(self):
        for allele in self.alleles:
            if len(allele) != 1:
                return 'INDEL'
        return 'SNP'

    @property
    def ref_allele(self):
        return self.alleles[0]

    @property
    def alt_alleles(self):
        return sorted(self.alleles[1:])

    @property
    def variant(self):
        return (self.position.chrom, self.position.chrom_start,
                self.ref_allele, self.alt_alleles)
