#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name='skin-mb-chen',
    version='0.1.0',
    description='Tools used in skin microbiome data analysis, Chen lab @UCSB',
    url='https://github.com/ynshen/skin_wound_microbiome',
    author='Yuning Shen',
    author_email='ynshen23@gmail.com',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
        'numpy',
        'pandas',
        'matplotlib',
        'biom-format',
        'biopython',
        'dendropy'
    ],
    zip_safe=False
)
