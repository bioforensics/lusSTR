# -------------------------------------------------------------------------------------------------
# Copyright (c) 2021, DHS.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import filecmp
import os
import pandas as pd
import pytest
import lusSTR
from lusSTR.tests import data_file
import re
from tempfile import NamedTemporaryFile


def test_uas_all(tmp_path):
    inputdb = data_file("snps")
    exp_out = data_file("snps_uas_all.txt")
    obs_out = str(tmp_path / "uas.txt")
    arglist = ["snps", inputdb, "-o", obs_out, "--type", "all", "--uas"]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.snps.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True


@pytest.mark.parametrize("type, lines", [("i", 189), ("p", 157)])
def test_uas_type(type, lines, tmp_path):
    inputdb = data_file("snps")
    obs_out = str(tmp_path / "uas.txt")
    arglist = ["snps", inputdb, "-o", obs_out, "--type", type, "--uas"]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.snps.main(args)
    with open(obs_out, "r") as fh:
        assert len(fh.readlines()) == lines


def test_sr_all(tmp_path):
    inputdb = data_file("snps")
    exp_out = data_file("snps_sr_all.txt")
    exp_out_full = data_file("snps_sr_all_full_output.txt")
    obs_out = str(tmp_path / "sr.txt")
    obs_out_full = str(tmp_path / "sr_full_output.txt")
    arglist = ["snps", inputdb, "-o", obs_out, "--type", "all"]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.snps.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True
    assert filecmp.cmp(exp_out_full, obs_out_full) is True


@pytest.mark.parametrize("type, lines, full_lines", [("i", 181, 2152), ("p", 158, 2982)])
def test_sr_type(type, lines, full_lines, tmp_path):
    inputdb = data_file("snps")
    obs_out = str(tmp_path / "sr.txt")
    obs_out_full = str(tmp_path / "sr_full_output.txt")
    arglist = ["snps", inputdb, "-o", obs_out, "--type", type]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.snps.main(args)
    with open(obs_out, "r") as fh:
        assert len(fh.readlines()) == lines
    with open(obs_out_full, "r") as fh:
        assert len(fh.readlines()) == full_lines
