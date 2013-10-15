#!/usr/bin/env python

from distutils.core import setup

setup(
    name='hgvs',
    version='0.8',
    description='HGVS name parsing and formatting',
    long_description = """\
This library provides a simple to use Python API for parsing, formatting, and
normalizing variant names specified in the standard recommended by the
Human Genome Variation Society (HGVS).""",
    author='Matt Rasmussen',
    author_email='rasmus@counsyl.com',

    packages=['hgvs', 'hgvs.tests'],
    scripts=[],
)
