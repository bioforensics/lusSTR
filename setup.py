#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from setuptools import setup

setup(
    name='lusSTR',
    description='Convert ForenSeq sequence strings to a compact representation',
    packages=['lusSTR'],
    package_data={
        'lusSTR': ['lusSTR/str_markers.json']
    },
    include_package_data=True,
    install_requires=['pandas>=1.0'],
    entry_points={
        'console_scripts': ['lusstr = lusSTR.cli:main']
    },
    classifiers=[
        'Environment :: Console',
        'Framework :: IPython',
        'Framework :: Jupyter',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Legal Industry',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
        'Topic : Scientific/Engineering :: Bio-Informatics',
    ],
    zip_safe=True,
)
