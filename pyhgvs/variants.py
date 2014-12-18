"""
Methods for manipulating genetic variants.
"""


from .models import Position


_COMP = dict(A='T', C='G', G='C', T='A', N='N',
             a='t', c='g', g='c', t='a', n='n')


def revcomp(seq):
    """Reverse complement."""
    return ''.join(_COMP[base] for base in reversed(seq))


def get_sequence(genome, chrom, start, end, is_forward_strand=True):
    """Return a sequence for the genomic region.

    Coordinates are 0-based, end-exclusive.
    """
    # Prevent fetching negative coordinates.
    start = max(start, 0)

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

    # No justification needed for empty indel.
    if len(indel) == 0:
        return start, end, indel

    if justify == 'left':
        while start > 0 and seq[start - 1] == indel[-1]:
            seq_added = seq[start - 1]
            indel = seq_added + indel[:-1]
            start -= 1
            end -= 1

    elif justify == 'right':
        while end < len(seq) and seq[end] == indel[0]:
            seq_added = seq[end]
            indel = indel[1:] + seq_added
            start += 1
            end += 1
    else:
        raise ValueError('unknown justify "%s"' % justify)
    return start, end, indel


def justify_genomic_indel(genome, chrom, start, end, indel, justify,
                          flank_length=20):
    """
    start, end: 0-based, end-exclusive coordinates of 'indel'.
    """
    ref_len = end - start

    while True:
        seq_start = max(start - flank_length, 0)
        indel_len = len(indel)
        fetch_len = indel_len + 2 * flank_length
        seq = get_sequence(
            genome, chrom, seq_start, seq_start + fetch_len)
        seq_end = seq_start + len(seq)
        if seq_end <= end and justify == 'right':
            # Indel is at end of chromosome, cannot justify right any further.
            return start, end, indel
        chrom_end = seq_end if seq_end < seq_start + fetch_len else 1e100

        # Get coordinates of indel within seq.
        indel_start = flank_length
        indel_end = flank_length + indel_len
        indel_start, indel_end, indel = justify_indel(
            indel_start, indel_end, indel, seq, justify)

        # Get indel coordinates with chrom.
        start = seq_start + indel_start
        end = start + ref_len
        if ((indel_start > 0 or seq_start == 0) and
                (indel_end < len(seq) or seq_end == chrom_end)):
            return start, end, indel
        # Since indel was justified to edge of seq, see if more justification
        # can be done.


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
    end = start + len(ref_sequence)
    position = Position(
        chrom=chrom,
        chrom_start=start,
        chrom_stop=end,
        is_forward_strand=True)
    return NormalizedVariant(position, ref_sequence, alt_sequences,
                             genome=genome)


class NormalizedVariant(object):
    """
    Normalizes variant representation to match GATK/VCF.
    """

    def __init__(self, position, ref_allele, alt_alleles,
                 seq_5p='', seq_3p='', genome=None):
        """
        position: a 0-index genomic Position.
        ref_allele: the reference allele sequence.
        alt_alleles: a list of alternate allele sequences.
        seq_5p: 5 prime flanking sequence of variant.
        seq_3p: 3 prime flanking sequence of variant.
        genome: a pygr compatible genome object (optional).
        """
        self.position = position
        self.alleles = [ref_allele] + list(alt_alleles)
        self.seq_5p = seq_5p
        self.seq_3p = seq_3p
        self.genome = genome
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
            self.position.chrom_stop -= common_suffix
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

            if self.genome:
                start, end, allele = justify_genomic_indel(
                    self.genome, self.position.chrom,
                    self.position.chrom_start, self.position.chrom_stop,
                    allele, 'left')
                self.position.chrom_start = start
                self.position.chrom_stop = end
                flank_length = 30
                self.seq_5p = get_sequence(self.genome, self.position.chrom,
                                           start - flank_length, start)
                self.seq_3p = get_sequence(self.genome, self.position.chrom,
                                           end, end + flank_length)
                self.alleles[i] = allele
            else:
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
        empty_seq = any(not allele for allele in self.alleles)
        uniq_starts = set(allele[0] for allele in self.alleles if allele)
        if empty_seq or len(uniq_starts) > 1:
            # Fetch more 5p flanking sequence if needed.
            if self.genome and self.seq_5p == '':
                start = self.position.chrom_start
                self.seq_5p = get_sequence(
                    self.genome, self.position.chrom, start - 5, start)

            self.log.append('1bp pad')
            if self.seq_5p:
                for i, allele in enumerate(self.alleles):
                    self.alleles[i] = self.seq_5p[-1] + self.alleles[i]

                self.seq_5p = self.seq_5p[:-1]
                self.position.chrom_start -= 1
            else:
                # According to VCF standard, if there is no 5prime sequence,
                # use 3prime sequence instead.
                assert self.seq_3p
                for i, allele in enumerate(self.alleles):
                    self.alleles[i] = self.alleles[i] + self.seq_3p[0]

                self.seq_3p = self.seq_3p[1:]
                self.position.chrom_stop += 1

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
