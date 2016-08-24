#!/usr/bin/env python

from distutils.core import setup

setup(name='fusionfusion',
      version='0.2.0beta',
      description='Python tools for extracting highly confident fusion transcripts from the results of several RNA-seq alignment tools.',
      author='Yuichi Shiraishi',
      author_email='friend1ws@gamil.com',
      url='https://github.com/friend1ws/fusionfusion',
      package_dir = {'': 'lib'},
      packages=['fusionfusion'],
      scripts=['fusionfusion'],
      license='GPL-3'
     )

