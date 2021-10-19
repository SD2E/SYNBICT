#!/usr/bin/env python

from setuptools import setup, find_packages

install_requires = [
    'biopython>=1.78',
    'dnaplotlib>=1.0',
    'flashtext>=2.7',
    'sbol2>=1.3',
    'matplotlib>=1.5.0'
]

setup(
    name='SYNBICT',
    version='1.5.1',
    description='Synthetic Biology Curation Tools (SYNBICT)',
    long_description='Synthetic Biology Curation Tools (SYNBICT) is a Python tool suite for automation-assisted annotation, curation, and functional inference for genetic designs.',
    author='Nicholas Roehner',
    author_email='nicholasroehner@gmail.com',
    url='https://github.com/SD2E/SYNBICT',
    download_url='https://github.com/SD2E/SYNBICT/archive/refs/tags/v1.4.tar.gz',
    packages=find_packages(),
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 3 :: Only"
    ]
)
