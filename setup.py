"""Setup module for FAMLI"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='famli',
    version='1.0',
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