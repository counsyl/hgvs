
from ..variants import revcomp

try:
    from pygr.seqdb import SequenceFileDB
except:
    pass


# Default genome file used for deriving mock sequences.
_genome_file = '/seq-data/seq/misc/genomics/genomes/hg19/hg19.fa'


class MockSequence(object):
    def __init__(self, sequence):
        self.sequence = sequence

    def __neg__(self):
        return MockSequence(revcomp(self.sequence))

    def __str__(self):
        return self.sequence


class MockChromosome(object):
    def __init__(self, name, genome=None):
        self.name = name
        self.genome = genome

    def __getslice__(self, start, end, step=1):
        return self.genome.get_seq(self.name, start, end)


class MockGenome(object):
    def __init__(self, lookup=None):
        self._chroms = {}
        self._lookup = lookup

    def __getitem__(self, chrom):
        if chrom in self._chroms:
            return self._chroms[chrom]
        else:
            chromosome = MockChromosome(chrom, self)
            self._chroms[chrom] = chromosome
            return chromosome

    def get_seq(self, chrom, start, end):
        return self._lookup[(chrom, start, end)]


class MockGenomeFile(MockGenome):
    def __init__(self, filename=_genome_file):
        super(MockGenomeFile, self).__init__()
        self._genome = SequenceFileDB(filename)

    def get_seq(self, chrom, start, end):
        seq = self._genome[chrom][start:end]
        print ((chrom, start, end), str(seq))
        return seq
