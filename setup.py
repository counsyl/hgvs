from setuptools import setup, Extension
import distutils.core
import sys
import os

description = ("This library provides a simple to use Python API for parsing, "
               "formatting, and normalizing variant names specified in the "
               "standard recommended by the Human Genome Variation Society "
               "(HGVS). This is a fork of the Counsyl package, updated to use "
               "pyfaidx for indexed fasta access.")

setup(
        name='hgvs',
        version='0.8',
        description='HGVS name parsing and formatting',
        long_description=description,
        author='Matthew Shirley',
        license = 'MIT',
        author_email='mdshw5@gmail.com',
        packages=['hgvs', 'hgvs.tests'],
        install_requires=['pyfaidx'],
        classifiers = [
                "Development Status :: 3 - Alpha",
                "License :: OSI Approved :: MIT License",
                "Environment :: Console",
                "Intended Audience :: Science/Research",
                "Natural Language :: English",
                "Operating System :: Unix",
                "Programming Language :: Python :: 3.3",
                "Topic :: Scientific/Engineering :: Bio-Informatics"
        ]
)