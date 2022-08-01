"""
Methods for manipulating HGVS names

Recommendations for the HGVS naming standard:
http://www.hgvs.org/mutnomen/standards.html

Definition of which transcript to use coding variants:
ftp://ftp.ncbi.nih.gov/refseq/H_sapiens/RefSeqGene/LRG_RefSeqGene

"""
from __future__ import absolute_import
from __future__ import unicode_literals

import re

from .models.cdna import CDNACoord, CDNA_START_CODON, CDNA_STOP_CODON
from .models.genome import GenomeSubset
from .models.hgvs_name import HGVSName, InvalidHGVSName
from .models.variants import justify_indel, normalize_variant, revcomp


def get_genomic_sequence(genome, chrom, start, end):
    """
    Return a sequence for the genomic region.

    start, end: 1-based, end-inclusive coordinates of the sequence.
    """
    if start > end:
        return ''
    else:
        return str(genome[str(chrom)][start - 1:end]).upper()


def get_allele(hgvs, genome, transcript=None):
    """Get an allele from a HGVSName, a genome, and a transcript."""
    chrom, start, end = hgvs.get_ref_coords(transcript)
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

    # Sometimes we need to retrieve alt from reference
    # Eg NC_000001.11:g.169549811=
    if hgvs.mutation_type == "=":
        alt = ref

    if hgvs.mutation_type in _indel_mutation_types:
        if hgvs.mutation_type == 'dup':
            # No alt supplied: NM_000492.3:c.1155_1156dup
            # Number used:     NM_004119.2(FLT3):c.1794_1811dup18
            # We *know* what the sequence is for "dup18", but not for "ins18"
            if not hgvs.alt_allele or re.match("^N+$", hgvs.alt_allele):
                alt = get_alt_from_sequence(hgvs, genome, transcript)

        # Left-pad alternate allele.
        alt = ref[0] + alt
    return chrom, start, end, ref, alt


def get_alt_from_sequence(hgvs, genome, transcript):
    """ returns ready for VCF """

    chrom, start, end = hgvs.get_raw_coords(transcript)
    return get_genomic_sequence(genome, chrom, start, end)


def matches_ref_allele(hgvs, genome, transcript=None):
    """Return True if reference allele matches genomic sequence."""
    is_forward_strand = transcript.tx_position.is_forward_strand if transcript else True
    ref, _ = hgvs.get_ref_alt(is_forward_strand, raw_dup_alleles=True)  # get raw values so dup isn't always True
    chrom, start, end = hgvs.get_raw_coords(transcript)
    genome_ref = get_genomic_sequence(genome, chrom, start, end)
    return genome_ref == ref


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
                    flank_length=30, normalize=True, lazy=False,
                    indels_start_with_same_base=True):
    """
    Parse an HGVS name into (chrom, start, end, ref, alt)

    hgvs_name: HGVS name to parse.
    genome: pygr compatible genome object.
    transcript: Transcript corresponding to HGVS name.
    normalize: If True, normalize allele according to VCF standard.
    lazy: If True, discard version information from incoming transcript/gene.
    indels_start_with_same_base: When normalizing, don't strip common prefix from indels
    """
    hgvs = HGVSName(hgvs_name)

    # Determine transcript.
    if hgvs.kind in ('c', 'n') and not transcript:
        if '.' in hgvs.transcript and lazy:
            hgvs.transcript, _ = hgvs.transcript.split('.')
        elif '.' in hgvs.gene and lazy:
            hgvs.gene, _ = hgvs.gene.split('.')
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

    chrom, start, _, ref, alt = get_vcf_allele(hgvs, genome, transcript)
    if normalize:
        nv = normalize_variant(chrom, start, ref, [alt], genome,
                               flank_length=flank_length,
                               indels_start_with_same_base=indels_start_with_same_base)
        chrom, start, ref, [alt] = nv.variant
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
    if mutation_type == 'ins':
        # Insert uses coordinates around the insert point.
        offset_start = offset - 1
        offset_end = offset
    else:
        offset_start = offset
        offset_end = offset + len(ref) - 1

    if not transcript:
        # Use genomic coordinate when no transcript is available.
        hgvs.kind = 'g'
        hgvs.start = offset_start
        hgvs.end = offset_end
    else:
        # Use cDNA coordinates.
        if transcript.is_coding:
            hgvs.kind = 'c'
        else:
            hgvs.kind = 'n'

        is_single_base_indel = (
            (mutation_type == 'ins' and len(alt) == 1) or
            (mutation_type in ('del', 'delins', 'dup') and len(ref) == 1))

        if mutation_type == '>' or (use_counsyl and is_single_base_indel):
            # Use a single coordinate.
            hgvs.cdna_start = transcript.genomic_to_cdna_coord(offset)
            hgvs.cdna_end = hgvs.cdna_start
        else:
            if transcript.strand == '-':
                offset_start, offset_end = offset_end, offset_start
            hgvs.cdna_start = transcript.genomic_to_cdna_coord(offset_start)
            hgvs.cdna_end = transcript.genomic_to_cdna_coord(offset_end)

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
