from __future__ import absolute_import
from __future__ import unicode_literals

import itertools
import os

from ..variants import revcomp


try:
    from pyfaidx import Genome as SequenceFileDB
    # Allow pyflakes to ignore redefinition in except clause.
    SequenceFileDB
except ImportError:
    SequenceFileDB = None


class MockGenomeError(Exception):
    pass


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

    def __getitem__(self, n):
        """Return sequence from region [start, end)

        Coordinates are 0-based, end-exclusive."""
        if isinstance(n, slice):
            return self.genome.get_seq(self.name, n.start, n.stop)
        else:
            return self.genome.get_seq(self.name, n, n+1)

    def __repr__(self):
        return 'MockChromosome("%s")' % (self.name)


class MockGenome(object):
    def __init__(self, lookup=None, filename=None, db_filename=None,
                 default_seq=None):
        """
        A mock genome object that provides a pygr compatible interface.

        lookup: a list of ((chrom, start, end), seq) values that define
            a lookup table for genome sequence requests.
        filename: a stream or filename containing a lookup table.
        db_filename: a fasta file to use for genome sequence requests.  All
            requests are recorded and can be writen to a lookup table file
            using the `write` method.
        default_seq: if given, this base will always be returned if
            region is unavailable.
        """
        self._chroms = {}
        self._lookup = lookup if lookup is not None else {}
        self._genome = None
        self._default_seq = default_seq

        if db_filename:
            # Use a real genome database.
            if SequenceFileDB is None:
                raise ValueError('pygr is not available.')
            self._genome = SequenceFileDB(db_filename)
        elif filename:
            # Read genome sequence from lookup table.
            self.read(filename)

    def __contains__(self, chrom):
        """Return True if genome contains chromosome."""
        return chrom in (self._genome or self._chroms)

    def __getitem__(self, chrom):
        """Return a chromosome by its name."""
        if chrom not in self._chroms:
            self._chroms[chrom] = MockChromosome(chrom, self)
        return self._chroms[chrom]

    def get_seq(self, chrom, start, end):
        """Return a sequence by chromosome name and region [start, end).

        Coordinates are 0-based, end-exclusive.
        """
        if self._genome:
            # Get sequence from real genome object and save result.
            seq = self._genome[chrom][start:end]
            self._lookup[(chrom, start, end)] = str(seq)
            return seq
        else:
            # Use lookup table to fetch genome sequence.
            try:
                return MockSequence(self._lookup[(chrom, start, end)])
            except KeyError:
                if self._default_seq:
                    # Generate default sequence.
                    return ''.join(itertools.islice(
                        itertools.cycle(self._default_seq),
                        None, end - start))
                else:
                    raise MockGenomeError(
                        'Sequence not in test data: %s:%d-%d' %
                        (chrom, start, end))

    def read(self, filename):
        """Read a sequence lookup table from a file.

        filename: a filename string or file stream.
        """
        if hasattr(filename, 'read'):
            infile = filename
        else:
            with open(filename) as infile:
                return self.read(infile)

        for line in infile:
            tokens = line.rstrip().split('\t')
            chrom, start, end, seq = tokens
            self._lookup[(chrom, int(start), int(end))] = seq
            if chrom not in self._lookup:
                self._chroms[chrom] = MockChromosome(chrom, self)

    def write(self, filename):
        """Write a sequence lookup table to file."""
        if hasattr(filename, 'write'):
            out = filename
        else:
            with open(filename, 'w') as out:
                return self.write(out)

        for (chrom, start, end), seq in self._lookup.items():
            out.write('\t'.join(map(str, [chrom, start, end, seq])) + '\n')


class MockGenomeTestFile(MockGenome):
    def __init__(self, lookup=None, filename=None, db_filename=None,
                 default_seq=None, create_data=False):
        if not create_data:
            db_filename = None
        super(MockGenomeTestFile, self).__init__(
            lookup=lookup, db_filename=db_filename,
            filename=filename,
            default_seq=default_seq)

        self._filename = filename
        self._create_data = (db_filename is not None)

        if self._create_data and os.path.exists(filename):
            # Clear output file when creating data.
            os.remove(filename)

    def get_seq(self, chrom, start, end):
        seq = super(MockGenomeTestFile, self).get_seq(chrom, start, end)

        # Save each query in append mode.
        if self._create_data:
            with open(self._filename, 'a') as out:
                out.write('\t'.join(map(str, [chrom, start, end, seq])) + '\n')
        return seq
