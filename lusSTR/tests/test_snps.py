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


@pytest.mark.parametrize(
    "input, filtering", [("snps_uas_all.txt", False), ("snps_uas_filtered.txt", True)]
)
def test_uas_all(input, filtering, tmp_path):
    inputdb = data_file("snps")
    exp_out = data_file(input)
    obs_out = str(tmp_path / "uas.txt")
    if filtering:
        arglist = ["config", "-w", str(tmp_path), "-o", "uas", "--input", inputdb, "--snps"]
    else:
        arglist = [
            "config",
            "-w",
            str(tmp_path),
            "-o",
            "uas",
            "--input",
            inputdb,
            "--snps",
            "--nofiltering",
        ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    convert_arglist = ["snps", "convert", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True


@pytest.mark.parametrize("type, lines", [("i", 131), ("p", 30), ("a", 69), ("p, a", 95)])
def test_uas_type(type, lines, tmp_path):
    inputdb = data_file("snps")
    obs_out = str(tmp_path / "uas.txt")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "uas",
        "--input",
        inputdb,
        "--snps",
        "--snp-type",
        type,
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    convert_arglist = ["snps", "convert", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    with open(obs_out, "r") as fh:
        assert len(fh.readlines()) == lines


def test_sr_all(tmp_path):
    inputdb = data_file("snps")
    exp_out = data_file("snps_sr_all.txt")
    exp_out_full = data_file("snps_sr_all_full_output.txt")
    obs_out = str(tmp_path / "sr.txt")
    obs_out_full = str(tmp_path / "sr_full_output.txt")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "sr",
        "--input",
        inputdb,
        "--snps",
        "--straitrazor",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    convert_arglist = ["snps", "convert", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True
    assert filecmp.cmp(exp_out_full, obs_out_full) is True


@pytest.mark.parametrize(
    "type, lines, full_lines",
    [("i", 181, 2152), ("p", 56, 1628), ("a", 107, 1381), ("p, a", 158, 3008)],
)
def test_sr_type(type, lines, full_lines, tmp_path):
    inputdb = data_file("snps")
    exp_out = data_file("snps_sr_all.txt")
    exp_out_full = data_file("snps_sr_all_full_output.txt")
    obs_out = str(tmp_path / "sr.txt")
    obs_out_full = str(tmp_path / "sr_full_output.txt")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "sr",
        "--input",
        inputdb,
        "--snps",
        "--straitrazor",
        "--snp-type",
        type,
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    convert_arglist = ["snps", "convert", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    with open(obs_out, "r") as fh:
        assert len(fh.readlines()) == lines
    with open(obs_out_full, "r") as fh:
        assert len(fh.readlines()) == full_lines
