#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from setuptools import setup
import versioneer

desc = 'Tool for converting NGS sequence data of forensic STR loci to various annotation styles'
setup(
    name='lusSTR',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=desc,
    packages=['lusSTR', 'lusSTR.tests'],
    package_data={
        'lusSTR': [
             'lusSTR/str_markers.json', 'lusSTR/tests/data/*', 'lusSTR/tests/data/STRait_Razor_test_output/*',
             'lusSTR/tests/data/UAS_bulk_input/*'
         ]
    },
    include_package_data=True,
    install_requires=['pandas>=1.0', 'xlrd=1.2.0'],
    entry_points={
        'console_scripts': ['lusstr = lusSTR.__main__:main']
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
