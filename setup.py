"""Setup module for FAMLI"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Set the long description
long_description = """
# FAMLI
Functional Analysis of Metagenomes by Likelihood Inference

Authors: 

  * Samuel Minot, Ph.D.
  * Jonathan Golob, M.D., Ph.D.

### Introduction

The goal of this work is to improve the accuracy of identifying protein-coding sequences
from short-read shotgun metagenomic sequencing data. The core challenge we consider here
is that of 'multi-mapping' reads – short nucleotide sequences that are equally similar to
multiple different reference protein sequences. In other domains such multi-mapping reads can
be dealt with in a variety of ways. For example, in the realm of taxonomic identification
it can be appropriate to assign them to the lowest common ancestor (LCA) of both references. 

However in the case of mapping short reads to a database of protein sequences (or peptides) we can not
assume that there is an underlying directed acyclic graph structure (e.g. a taxonomy). Peptides
can evolve by duplication events, homologous recombination, and other means of sharing highly conserved
domains (leading to shared short reads). If one simply includes all peptides for which there is a read,
we find the false positives outnumber the true positive by as much as 1000:1. 

We developed a method to iteratively assign shared reads to the most likely true peptides, bringing the 
precision (TP / (TP+FP)) to close to 90%. To do so, we used the following principles:

"""

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='famli',
    version='1.2',
    description='Functional Analysis of Metagenomes by Likelihood Inferrence',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/FredHutch/FAMLI',
    author='Samuel Minot',
    author_email='sminot@fredhutch.org',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
    ],
    keywords='metagenomics microbiome science',
    packages=find_packages(exclude=['tests']),
    install_requires=[
        "pandas>=0.20.3",
        "biopython>=1.70",
        "numpy>=1.13.1",
        "scipy>=0.19.1",
        "awscli>=1.11.146",
        "boto3>=1.4.7",
        "python-dateutil==2.6.0"
    ],
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/fredhutch/famli/issues',
        'Source': 'https://github.com/fredhutch/famli/',
    },
    entry_points={
        'console_scripts': [
            'famli = famli.run_famli:main',
        ],
    },
)
