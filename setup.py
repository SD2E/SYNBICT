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
    version='1.4',
    packages=find_packages(),
    install_requires=install_requires,
    classifiers=[
        "Programming Language :: Python :: 3 :: Only"
    ]
)
