#!/usr/bin/env python

from distutils.core import setup

setup(
    name='hgvs',
    version='0.8',
    description='HGVS name parsing and formatting',
    long_description = """
            """,
    author='Matt Rasmussen',
    author_email='rasmus@counsyl.com',

    packages=['hgvs', 'hgvs.tests'],
    scripts=[],
)
