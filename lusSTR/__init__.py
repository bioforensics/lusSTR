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

from pkg_resources import resource_filename
from lusSTR import cli
from lusSTR._version import get_versions

__version__ = get_versions()["version"]
del get_versions


def snakefile(workflow="strs"):
    return resource_filename("lusSTR", f"workflows/{workflow}.smk")


def wrapper(label):
    return resource_filename("lusSTR", f"wrappers/{label}.py")

from . import _version
__version__ = _version.get_versions()['version']
