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