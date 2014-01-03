#!/usr/bin/env python
"""
Example usage of HGVS library.

To run the script, first begin by opening a terminal in the root directory
of this software package. Note, the root directory should contain
`setup.py`.

Second, obtain genome sequence in FASTA format, which is required in
example. Genome sequence can be fetched using the following commands:

  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
  tar zxvf chromFa.tar.gz
  cat chr*.fa > hg19.fa
  rm chr*.fa chromFa.tar.gz

This example script can be run using:

  python examples/example1.py

The following output should be displayed:

  chr11 17496508 T C
  NM_000352.3(ABCC8):c.215A>G
  ('NM_000352.3', 'c', '>', CDNACoord(215, -10), CDNACoord(215, -10), 'A', 'G')

"""

import hgvs
import hgvs.utils
from pygr.seqdb import SequenceFileDB

# Read genome sequence using pygr.
genome = SequenceFileDB('hg19.fa')

# Read RefSeq transcripts into a python dict.
with open('hgvs/data/genes.refGene') as infile:
    transcripts = hgvs.utils.read_transcripts(infile)


# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    return transcripts.get(name)


# Parse the HGVS name into genomic coordinates and alleles.
chrom, offset, ref, alt = hgvs.parse_hgvs_name(
    'NM_000352.3:c.215A>G', genome, get_transcript=get_transcript)
print chrom, offset, ref, alt
# Returns variant in VCF style: ('chr11', 17496508, 'T', 'C')
# Notice that since the transcript is on the negative strand, the alleles
# are reverse complemented during conversion.


# Format an HGVS name.
chrom, offset, ref, alt = ('chr11', 17496508, 'T', 'C')
transcript = get_transcript('NM_000352.3')
hgvs_name = hgvs.format_hgvs_name(
    chrom, offset, ref, alt, genome, transcript)
print hgvs_name
# Returns 'NM_000352.3(ABCC8):c.215A>G'


hgvs_name = hgvs.HGVSName('NM_000352.3:c.215-10A>G')
# fields of the HGVS name are available as attributes:
#
# hgvs_name.transcript = 'NM_000352.3'
# hgvs_name.kind = 'c'
# hgvs_name.mutation_type = '>'
# hgvs_name.cdna_start = hgvs.CDNACoord(215, -10)
# hgvs_name.cdna_end = hgvs.CDNACoord(215, -10)
# hgvs_name.ref_allele = 'A'
# hgvs_name.alt_allele = 'G'

print (hgvs_name.transcript,
       hgvs_name.kind,
       hgvs_name.mutation_type,
       hgvs_name.cdna_start,
       hgvs_name.cdna_end,
       hgvs_name.ref_allele,
       hgvs_name.alt_allele)
