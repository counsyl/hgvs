#!/usr/bin/env python

from setuptools import setup
import sys

description = ("This library provides a simple to use Python API for parsing, "
               "formatting, and normalizing variant names specified in the "
               "standard recommended by the Human Genome Variation Society "
               "(HGVS).")


def main():
    python_version = sys.version_info
    if python_version < (2, 6):
        print ("This library requires Python version >=2.6, "
               "You have version %d.%d" % python_version[:2])
        sys.exit(1)

    setup(
        name='pyhgvs',
        version='0.9.5',
        description='HGVS name parsing and formatting',
        long_description=description,
        author='Matt Rasmussen',
        author_email='rasmus@counsyl.com',
        packages=['pyhgvs', 'pyhgvs.tests'],
        include_package_data=True,
        package_data={
            '': ['requirements-dev.txt'],
        },
        scripts=[],
        tests_require=[l.strip() for l in open('requirements-dev.txt') if l],
    )

if __name__ == '__main__':
    main()
