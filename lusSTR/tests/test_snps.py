#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2021, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import filecmp
import os
import pandas as pd
import pytest
import lusSTR
from lusSTR.tests import data_file
import re
from tempfile import NamedTemporaryFile


def test_uas_all(tmp_path):
    inputdb = data_file('snps')
    exp_out = data_file('snps_uas_all.txt')
    obs_out = str(tmp_path / 'uas.txt')
    arglist = ['snps', inputdb, '-o', obs_out, '--type', 'all', '--uas']
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.snps.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True
