#!/usr/bin/env python

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

def get_version():
    with open(path.join(here, "fusionfusion/version.py"), encoding = 'utf-8') as hin:
        for line in hin:
            if line.startswith("__version__"):
                version = line.partition('=')[2]
                return version.strip().strip('\'"')
    raise ValueError('Could not find version.')

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name = 'fusionfusion',
    version = get_version(),
    description='Python tools for extracting highly confident fusion transcripts from the results of several RNA-seq alignment tools.',
    long_description=long_description, 
    long_description_content_type='text/markdown',  
    url = 'https://github.com/friend1ws/fusionfusion',
    author = 'Yuichi Shiraishi',
    author_email = 'friend1ws@gamil.com',
    license = 'GPLv3',

    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],

    packages = find_packages(exclude = ['tests']),

    install_requires = [],
    entry_points = {'console_scripts': ['fusionfusion = fusionfusion:main']}

)

