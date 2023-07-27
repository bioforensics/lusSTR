# -------------------------------------------------------------------------------------------------
# Copyright (c) 2020, DHS.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import glob
from setuptools import setup
import versioneer

desc = "Tool for processing NGS sequence data of forensic STR loci and SNPs for use in probabilistic genotyping software"
setup(
    name="lusSTR",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=desc,
    packages=["lusSTR", "lusSTR.cli", "lusSTR.tests"],
    package_data={
        "lusSTR": [
            "lusSTR/data/*",
            "lusSTR/tests/data/*",
            "lusSTR/tests/data/STRait_Razor_test_output/*",
            "lusSTR/tests/data/UAS_bulk_input/*",
            "lusSTR/tests/data/snps/*",
            "lusSTR/tests/data/RU_stutter_test/*",
            "lusSTR/tests/data/NGS_stutter_test/*",
            "lusSTR/tests/data/kinsnps/*",
            "lusSTR/tests/data/lusstr_output/*",
            "lusSTR/workflows/*",
            "lusSTR/wrappers/*",
        ]
    },
    include_package_data=True,
    install_requires=["pandas>=1.0,<2.0", "openpyxl>=3.0.6", "snakemake>=7.22.0", "pyyaml>=6.0"],
    entry_points={"console_scripts": ["lusstr = lusSTR.cli:main"]},
    scripts=glob.glob("lusSTR/scripts/*"),
    classifiers=[
        "Environment :: Console",
        "Framework :: IPython",
        "Framework :: Jupyter",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Legal Industry",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
        "Topic : Scientific/Engineering :: Bio-Informatics",
    ],
    zip_safe=True,
)
