#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import os
from pkg_resources import resource_filename


def data_file(path):
    pathparts = path.split('/')
    relpath = os.path.join('tests', 'data', *pathparts)
    return resource_filename('lusSTR', relpath)