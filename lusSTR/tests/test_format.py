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

import filecmp
import lusSTR
from lusSTR.tests import data_file
import pytest
import os
from shutil import copytree
from tempfile import NamedTemporaryFile


def test_format(tmp_path):
    UAStestfile = data_file("snps/Positive Control Sample Details Report 2315.xlsx")
    exp_output = data_file("testformat.csv")
    obs_output = str(tmp_path / "lusstr_output.csv")
    str_path = str(tmp_path)
    config_arglist = ["config", "--input", str(UAStestfile), "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(config_arglist))
    format_arglist = ["strs", "format", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    assert filecmp.cmp(exp_output, obs_output) is True


@pytest.mark.parametrize(
    "input, testoutput",
    [
        ("STRait_Razor_test_output/", "STRait_Razor_test_output.csv"),
        ("STRait_Razor_test_output/A001.txt", "STRaitRazor_output_test_A001.csv"),
    ],
)
def test_format_straitrazor(input, testoutput, tmp_path):
    input_file = data_file(input)
    exp_output = data_file(testoutput)
    obs_output = str(tmp_path / "lusstr_output.csv")
    str_path = str(tmp_path)
    config_arglist = ["config", "--input", str(input_file), "-w", str_path, "--straitrazor"]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(config_arglist))
    format_arglist = ["strs", "format", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    assert filecmp.cmp(exp_output, obs_output) is True


def test_format_sexloci_uas(tmp_path):
    UAStestfile = data_file("snps/Positive Control Sample Details Report 2315.xlsx")
    exp_output = data_file("testformat_uas_sexloci.csv")
    obs_output = str(tmp_path / "lusstr_output_sexloci.csv")
    str_path = str(tmp_path)
    config_arglist = ["config", "--input", str(UAStestfile), "-w", str_path, "--sex"]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(config_arglist))
    format_arglist = ["strs", "format", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    assert filecmp.cmp(exp_output, obs_output) is True


def test_format_sex_loci_straitrazor(tmp_path):
    inputdb = data_file("STRait_Razor_test_output")
    exp_output = data_file("testformat_sr_sexloci.csv")
    obs_output = str(tmp_path / "lusstr_output_sexloci.csv")
    str_path = str(tmp_path)
    config_arglist = ["config", "--input", str(inputdb), "-w", str_path, "--straitrazor", "--sex"]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(config_arglist))
    format_arglist = ["strs", "format", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    assert filecmp.cmp(exp_output, obs_output) is True


def test_uas_directory_with_xy(tmp_path):
    inputdb = data_file("UAS_bulk_input")
    exp_output_sex = data_file("UAS_bulk_test_sexloci.csv")
    obs_sex_output = str(tmp_path / "lusstr_output_sexloci.csv")
    exp_output = data_file("UAS_bulk_test.csv")
    obs_output = str(tmp_path / "lusstr_output.csv")
    str_path = str(tmp_path)
    config_arglist = ["config", "--input", str(inputdb), "-w", str_path, "--sex"]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(config_arglist))
    format_arglist = ["strs", "format", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    assert filecmp.cmp(exp_output_sex, obs_sex_output) is True
    assert filecmp.cmp(exp_output, obs_output) is True
