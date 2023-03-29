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

from setuptools import setup
import versioneer

desc = "Tool for converting NGS sequence data of forensic STR loci to various annotation styles"
setup(
    name="lusSTR",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description=desc,
    packages=["lusSTR", "lusSTR.tests"],
    package_data={
        "lusSTR": [
            "lusSTR/str_markers.json",
            "lusSTR/snp_data.json",
            "lusSTR/tests/data/*",
            "lusSTR/tests/data/STRait_Razor_test_output/*",
            "lusSTR/tests/data/UAS_bulk_input/*",
            "lusSTR/tests/data/snps/*",
            "lusSTR/tests/data/RU_stutter_test/*",
            "lusSTR/filters.json",
            "lusSTR/tests/data/NGS_stutter_test/*",
        ]
    },
    include_package_data=True,
    install_requires=["pandas>=1.0", "openpyxl>=3.0.6"],
    entry_points={"console_scripts": ["lusstr = lusSTR.cli:main"]},
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
