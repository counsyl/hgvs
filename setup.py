#!/usr/bin/env python
from __future__ import print_function
from __future__ import unicode_literals

from setuptools import setup
import sys

description = ("This library provides a simple to use Python API for parsing, "
               "formatting, and normalizing variant names specified in the "
               "standard recommended by the Human Genome Variation Society "
               "(HGVS).")


def main():
    python_version = sys.version_info
    if python_version < (2, 6):
        print(("This library requires Python version >=2.6, "
               "You have version %d.%d" % python_version[:2]))
        sys.exit(1)

    setup(
        name='pyhgvs',
        version='0.11.1',
        description='HGVS name parsing and formatting',
        long_description=description,
        author='Matt Rasmussen',
        author_email='rasmus@counsyl.com',
        packages=[str('pyhgvs'), str('pyhgvs.tests')],
        include_package_data=True,
        scripts=[],
        install_requires=['pip>=1.2'],
        tests_require=[
          'flake8==2.2.5',
          'nose==1.3.7',
          'pyfaidx==0.5.0.1',
        ],
    )


if __name__ == '__main__':
    main()
