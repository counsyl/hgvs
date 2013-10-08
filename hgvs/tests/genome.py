
from ..variants import revcomp

try:
    from pygr.seqdb import SequenceFileDB
except:
    SequenceFileDB
    SequenceFileDB = None


class MockSequence(object):
    def __init__(self, sequence):
        self.sequence = sequence

    def __neg__(self):
        """Return reverse complement sequence."""
        return MockSequence(revcomp(self.sequence))

    def __str__(self):
        return self.sequence

    def __repr__(self):
        return 'MockSequence("%s")' % self.sequence


class MockChromosome(object):
    def __init__(self, name, genome=None):
        self.name = name
        self.genome = genome

    def __getslice__(self, start, end, step=1):
        """Return sequence from region [start, end)

        Coordinates are 0-based, end-exclusive."""
        return self.genome.get_seq(self.name, start, end)

    def __repr__(self):
        return 'MockChromosome("%s")' % (self.name)


class MockGenome(object):
    def __init__(self, lookup=None, filename=None, db_filename=None):
        """
        A mock genome object that provides a pygr compatiable interface.

        lookup: a list of ((chrom, start, end), seq) values that define
            a lookup table for genome sequence requests.
        filename: a stream or filename containing a lookup table.
        db_filename: a fasta file to use for genome sequence requests.  All
            requests are recorded and can be writen to a lookup table file
            using the `write` method.
        """
        self._chroms = {}
        self._lookup = lookup if lookup is not None else {}
        self._genome = None

        if db_filename:
            # Use a real genome database.
            if SequenceFileDB is None:
                raise AssertionError('pygr is not available.')
            self._genome = SequenceFileDB(db_filename)
        elif filename:
            # Read genome sequence from lookup table.
            self.read(filename)

    def __getitem__(self, chrom):
        """Return a chromosome by its name."""
        if chrom in self._chroms:
            return self._chroms[chrom]
        else:
            chromosome = MockChromosome(chrom, self)
            self._chroms[chrom] = chromosome
            return chromosome

    def get_seq(self, chrom, start, end):
        """Return a sequence by chromosome name and region [start, end).

        Coordinates are 0-based, end-exclusive."""
        if self._genome:
            # Get sequence from real genome object and save result.
            seq = self._genome[chrom][start:end]
            self._lookup[(chrom, start, end)] = str(seq)
            return seq
        else:
            # Use lookup table to fetch genome sequence.
            return MockSequence(self._lookup[(chrom, start, end)])

    def read(self, filename):
        """Read a sequence lookup table from a file."""
        if isinstance(filename, basestring):
            with open(filename) as infile:
                return self.read(infile)
        else:
            infile = filename

        for line in infile:
            tokens = line.rstrip().split('\t')
            chrom, start, end, seq = tokens
            self._lookup[(chrom, int(start), int(end))] = seq

    def write(self, filename):
        """Write a sequence lookup table to file."""
        if isinstance(filename, basestring):
            with open(filename, 'w') as out:
                return self.write(out)
        else:
            out = filename

        for (chrom, start, end), seq in self._lookup.iteritems():
            out.write('\t'.join(map(str, [chrom, start, end, seq])) + '\n')
