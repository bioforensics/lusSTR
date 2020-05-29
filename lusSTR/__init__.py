#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from lusSTR import annot
from lusSTR import marker
from lusSTR import repeat
from lusSTR import format
from lusSTR import cli
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
