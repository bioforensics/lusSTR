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
    "input, filtering", [("snps_uas_all.tsv", False), ("snps_uas_filtered.tsv", True)]
)
def test_uas_all(input, filtering, tmp_path):
    inputdb = data_file("snps")
    exp_out = data_file(input)
    obs_out = str(tmp_path / "uas.tsv")
    if filtering:
        arglist = ["config", "-w", str(tmp_path), "-o", "uas", "--input", str(inputdb), "--snps"]
    else:
        arglist = [
            "config",
            "-w",
            str(tmp_path),
            "-o",
            "uas",
            "--input",
            str(inputdb),
            "--snps",
            "--nofiltering",
        ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    format_arglist = ["snps", "format", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True


@pytest.mark.parametrize("type, lines", [("i", 131), ("p", 30), ("a", 69), ("p, a", 95)])
def test_uas_type(type, lines, tmp_path):
    inputdb = data_file("snps")
    obs_out = str(tmp_path / "uas.tsv")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "uas",
        "--input",
        str(inputdb),
        "--snps",
        "--snp-type",
        type,
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    format_arglist = ["snps", "format", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    with open(obs_out, "r") as fh:
        assert len(fh.readlines()) == lines


def test_sr_all(tmp_path):
    inputdb = data_file("snps")
    exp_out = data_file("snps_sr_all.tsv")
    exp_out_full = data_file("snps_sr_all_full_output.tsv")
    obs_out = str(tmp_path / "sr.tsv")
    obs_out_full = str(tmp_path / "sr_full_output.tsv")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "sr",
        "--input",
        str(inputdb),
        "--snps",
        "--analysis-software",
        "straitrazor",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    format_arglist = ["snps", "format", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True
    assert filecmp.cmp(exp_out_full, obs_out_full) is True


@pytest.mark.parametrize(
    "type, lines, full_lines",
    [("i", 181, 2152), ("p", 56, 1628), ("a", 107, 1381), ("p, a", 158, 3008)],
)
def test_sr_type(type, lines, full_lines, tmp_path):
    inputdb = data_file("snps")
    exp_out = data_file("snps_sr_all.tsv")
    exp_out_full = data_file("snps_sr_all_full_output.tsv")
    obs_out = str(tmp_path / "sr.tsv")
    obs_out_full = str(tmp_path / "sr_full_output.tsv")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "sr",
        "--input",
        str(inputdb),
        "--snps",
        "--analysis-software",
        "straitrazor",
        "--snp-type",
        type,
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    format_arglist = ["snps", "format", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    with open(obs_out, "r") as fh:
        assert len(fh.readlines()) == lines
    with open(obs_out_full, "r") as fh:
        assert len(fh.readlines()) == full_lines


@pytest.mark.parametrize(
    "output, filtering",
    [("kinsnps/snps_kin_all.tsv", False), ("kinsnps/snps_kin_filtered.tsv", True)],
)
def test_kintelligence(output, filtering, tmp_path):
    inputdb = data_file("kinsnps/Kin_pos_1ng_reference Sample Report 2023_07_11_13_16_31.xlsx")
    exp_out = data_file(output)
    obs_out = str(tmp_path / "kin.tsv")
    if filtering:
        arglist = [
            "config",
            "-w",
            str(tmp_path),
            "-o",
            "kin",
            "--input",
            str(inputdb),
            "--snps",
            "--kintelligence",
        ]
    else:
        arglist = [
            "config",
            "-w",
            str(tmp_path),
            "-o",
            "kin",
            "--input",
            str(inputdb),
            "--snps",
            "--kintelligence",
            "--nofiltering",
        ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    format_arglist = ["snps", "format", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True


def test_kintelligence_all(tmp_path):
    inputdb = data_file("kinsnps/")
    evid_exp_output = data_file("kinsnps/evidence.tsv")
    ref_exp_output = data_file("kinsnps/reference.tsv")
    evid_obs_output = f"{str(tmp_path)}/kin_snp_evidence.tsv"
    ref_obs_output = f"{str(tmp_path)}/kin_snp_reference.tsv"
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "kin",
        "--input",
        str(inputdb),
        "--snps",
        "--kintelligence",
        "--snp-reference",
        "Kin_pos_reference",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    all_arglist = ["snps", "all", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(all_arglist))
    assert filecmp.cmp(evid_exp_output, evid_obs_output) is True
    assert filecmp.cmp(ref_exp_output, ref_obs_output) is True


def test_multiple_reference_profiles(tmp_path):
    inputdb = data_file("kinsnps/")
    exp_out = data_file("kinsnps/multiplerefs.tsv")
    obs_out = str(tmp_path / "kin_snp_reference.tsv")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "kin",
        "--input",
        str(inputdb),
        "--snps",
        "--kintelligence",
        "--snp-reference",
        "Kin_pos_reference, Kin_pos_1ng",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    all_arglist = ["snps", "all", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(all_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True


def test_snp_bins(tmp_path):
    input_sample = data_file("kinsnps/Kin_pos_1ng Sample Report 2023_07_11_13_16_31.xlsx")
    exp_out = data_file("kinsnps/Kin_pos_1ng_snpsetscombined_evidence.tsv")
    obs_out = str(tmp_path / "evidence_samples/Kin_pos_1ng_snpsetscombined_evidence.tsv")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "kin",
        "--input",
        str(input_sample),
        "--snps",
        "--kintelligence",
        "--separate",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    all_arglist = ["snps", "all", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(all_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True
    for snp_set in range(1, 10):
        path = tmp_path / f"evidence_samples/Kin_pos_1ng_set{snp_set}.tsv"
        assert path.is_file()


def test_uas_version5(tmp_path):
    input_sample = data_file("NA24385 Sample Report 2023_09_07_15_11_11.xlsx")
    exp_out = data_file("kinsnps/uasversion_snp_evidence.tsv")
    obs_out = str(tmp_path / "kin_v5_snp_evidence.tsv")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "kin_v5",
        "--input",
        str(input_sample),
        "--snps",
        "--kintelligence",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    all_arglist = ["snps", "all", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(all_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True


def test_snp_logs(tmp_path):
    inputdb = data_file("snps")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "-o",
        "sr",
        "--input",
        str(inputdb),
        "--snps",
        "--analysis-software",
        "straitrazor",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    format_arglist = ["snps", "format", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(format_arglist))
    for dir in os.listdir(str(tmp_path / "logs/")):
        assert os.path.exists(str(tmp_path / f"logs/{dir}/snp_config.yaml"))
        assert os.path.exists(
            str(tmp_path / f"logs/{dir}/input/Positive Control Sample Details Report 2315.xlsx")
        )
