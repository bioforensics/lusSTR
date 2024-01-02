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
import importlib.resources
import os
import pandas as pd
import pytest
import lusSTR
from lusSTR.scripts.marker import STRMarkerObject
from lusSTR.scripts.repeat import reverse_complement
from lusSTR.tests import data_file
from pathlib import Path
import re
from tempfile import NamedTemporaryFile
import shutil
import yaml


def test_split_sequence_into_two_strings():
    sequence = "TAGATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGGTGTGTGTGTGTG"
    reverse_comp_sequence = reverse_complement(sequence)
    repeat_for_split = "CACA"
    seq1, seq2 = lusSTR.scripts.repeat.split_sequence_into_two_strings(
        reverse_comp_sequence, repeat_for_split
    )
    assert seq1 == "CACACACACACA"
    assert seq2 == "CCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATCTA"


@pytest.mark.parametrize(
    "infile, len_sum, len_uncom, xy_len_sum, xy_len_uncom, pwrseq",
    [
        ("testformat_sr", 897, 913, 9701, 16108, False),
        ("powerseq", 359, 441, 256, 303, True),
    ],
)
def test_convert_full_nocombine(
    infile, len_sum, len_uncom, xy_len_sum, xy_len_uncom, pwrseq, tmp_path
):
    str_path = str(tmp_path / "WD")
    inputfile = data_file(f"{infile}.csv")
    inputfile_sex = data_file(f"{infile}_sexloci.csv")
    obs_out = f"{infile}.csv"
    obs_sex_out = f"{infile}_sexloci.csv"
    if pwrseq is False:
        arglist = [
            "config",
            "-w",
            str_path,
            "--analysis-software",
            "straitrazor",
            "--nocombine",
            "--sex",
            "-o",
            infile,
            "--input",
            "WD",
        ]
    else:
        arglist = [
            "config",
            "-w",
            str_path,
            "--analysis-software",
            "straitrazor",
            "--nocombine",
            "--sex",
            "-o",
            infile,
            "--input",
            "WD",
            "--powerseq",
        ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    shutil.copyfile(inputfile_sex, os.path.join(str_path, obs_sex_out))
    shutil.copyfile(inputfile, os.path.join(str_path, obs_out))
    convert_arglist = ["strs", "convert", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    with open(f"{str_path}/{infile}_no_combined_reads.txt", "r") as fh:
        assert len(fh.readlines()) == len_uncom
    with open(f"{str_path}/{infile}_sexloci_no_combined_reads.txt", "r") as fh:
        assert len(fh.readlines()) == xy_len_uncom
    with open(f"{str_path}/{infile}.txt", "r") as fh:
        assert len(fh.readlines()) == len_sum
    with open(f"{str_path}/{infile}_sexloci.txt", "r") as fh:
        assert len(fh.readlines()) == xy_len_sum


def test_flanks(tmp_path):
    inputfile = data_file("test_flank.csv")
    exp_out = data_file("testflanks_flanks.txt")
    str_path = str(tmp_path / "WD")
    obs_out = str(tmp_path / "WD/testflanks_flanks.txt")
    arglist = [
        "config",
        "-w",
        str_path,
        "-o",
        "testflanks",
        "--analysis-software",
        "straitrazor",
        "--input",
        "WD",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    shutil.copyfile(inputfile, os.path.join(str_path, "testflanks.csv"))
    convert_arglist = ["strs", "convert", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True


@pytest.mark.parametrize(
    "input, exp_length", [("Flanks_testing_file.csv", 952), ("test_FGA_short_seq.csv", 2)]
)
def test_convert_combine(input, exp_length, tmp_path):
    inputfile = data_file(input)
    str_path = str(tmp_path / "WD")
    obs_out = str(tmp_path / "WD/testflanks.txt")
    arglist = [
        "config",
        "-w",
        str_path,
        "-o",
        "testflanks",
        "--analysis-software",
        "straitrazor",
        "--input",
        "WD",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    shutil.copyfile(inputfile, os.path.join(str_path, "testflanks.csv"))
    convert_arglist = ["strs", "convert", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    with open(obs_out, "r") as fh:
        assert len(fh.readlines()) == exp_length


@pytest.mark.parametrize(
    "locus, sequence, software, kit, output",
    [
        (
            "CSF1PO",
            "CTTCCTATCTATCTATCTATCTAATCTATCTATCTT",
            "straitrazor",
            "forenseq",
            "Possible indel or partial sequence",
        ),
        (
            "DYS393",
            "AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATGTATGTCTTTTCTATGAGACATACC",
            "straitrazor",
            "powerseq",
            "UAS region indicates entire sequence; Possible indel or partial sequence",
        ),
    ],
)
def test_indel_flag(locus, sequence, software, kit, output):
    marker = STRMarkerObject(locus, sequence, software, kit=kit)
    assert marker.indel_flag == output


def test_powerseq_flanks(tmp_path):
    inputfile = data_file("powerseq_flanking_convert_test.csv")
    exp_out = data_file("powerseq_flanking_convert_test_flanks.txt")
    str_path = str(tmp_path / "WD")
    obs_out = str(tmp_path / "WD/powerseq_flanks.txt")
    arglist = [
        "config",
        "-w",
        str_path,
        "-o",
        "powerseq",
        "--analysis-software",
        "straitrazor",
        "--input",
        "WD",
        "--powerseq",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    shutil.copyfile(inputfile, os.path.join(str_path, "powerseq.csv"))
    convert_arglist = ["strs", "convert", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    assert filecmp.cmp(exp_out, obs_out) is True


def test_convert_uas_sexloci(tmp_path):
    str_path = str(tmp_path / "WD")
    inputfile = data_file("testformat_uas.csv")
    inputfile_sex = data_file("testformat_uas_sexloci.csv")
    exp_sex_out = data_file("testformat_uas_sexloci.txt")
    obs_sex_out = str(tmp_path / "WD/testformatuas_sexloci.txt")
    arglist = ["config", "-w", str_path, "-o", "testformatuas", "--sex", "--input", "WD"]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    shutil.copyfile(inputfile_sex, os.path.join(str_path, "testformatuas_sexloci.csv"))
    shutil.copyfile(inputfile, os.path.join(str_path, "testformatuas.csv"))
    convert_arglist = ["strs", "convert", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    assert filecmp.cmp(exp_sex_out, obs_sex_out) is True


@pytest.mark.parametrize(
    "input, testoutput, flank_output, kit",
    [
        (
            "testformat_sr",
            "testformat_sr_sexloci.txt",
            "testformat_sr_sexloci_flanks.txt",
            "forenseq",
        ),
        (
            "powerseq_flanking_convert_test",
            "powerseq_flanking_convert_test_sexloci.txt",
            "powerseq_flanking_convert_test_sexloci_flanks.txt",
            "powerseq",
        ),
    ],
)
def test_convert_sr_sexloci(input, testoutput, flank_output, kit, tmp_path):
    str_path = str(tmp_path / "WD")
    inputfile = data_file(f"{input}.csv")
    inputfile_sex = data_file(f"{input}_sexloci.csv")
    exp_sex_out = data_file(testoutput)
    exp_sex_flank_out = data_file(flank_output)
    obs_sex_out = str(tmp_path / "WD/testformatsr_sexloci.txt")
    obs_sex_flank_out = str(tmp_path / "WD/testformatsr_sexloci_flanks.txt")
    if kit == "forenseq":
        arglist = [
            "config",
            "-w",
            str_path,
            "-o",
            "testformatsr",
            "--sex",
            "--input",
            "WD",
            "--analysis-software",
            "straitrazor",
        ]
    else:
        arglist = [
            "config",
            "-w",
            str_path,
            "-o",
            "testformatsr",
            "--sex",
            "--input",
            "WD",
            "--analysis-software",
            "straitrazor",
            "--powerseq",
        ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    shutil.copyfile(inputfile_sex, os.path.join(str_path, "testformatsr_sexloci.csv"))
    shutil.copyfile(inputfile, os.path.join(str_path, "testformatsr.csv"))
    convert_arglist = ["strs", "convert", "-w", str_path]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(convert_arglist))
    assert filecmp.cmp(exp_sex_out, obs_sex_out) is True
    assert filecmp.cmp(exp_sex_flank_out, obs_sex_flank_out) is True


def test_config(tmp_path):
    obs_config = str(tmp_path / "config.yaml")
    exp_config = importlib.resources.files("lusSTR") / "data/config.yaml"
    arglist = ["config", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    assert os.path.exists(obs_config)


def test_config_settings(tmp_path):
    obs_config = str(tmp_path / "config.yaml")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "--analysis-software",
        "straitrazor",
        "--reference",
        "--str-type",
        "ce",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    with open(obs_config, "r") as file:
        data = yaml.safe_load(file)
    assert data["analysis_software"] == "straitrazor"
    assert data["data_type"] == "ce"
    assert data["profile_type"] == "reference"


@pytest.mark.parametrize(
    "command, output, format_out, convert_out, all_out",
    [
        ("format", "lusstr_output.csv", True, False, False),
        ("convert", "lusstr_output.txt", True, True, False),
        ("all", "lusstr_output/Positive_Control_evidence_ngs.csv", True, True, True),
    ],
)
def test_snakemake(command, output, format_out, convert_out, all_out, tmp_path):
    config = str(tmp_path / "config.yaml")
    inputfile = data_file("UAS_bulk_input/Positive Control Sample Details Report 2315.xlsx")
    exp_output = data_file(output)
    obs_output = str(tmp_path / output)
    obs_format_output = str(tmp_path / "lusstr_output.csv")
    obs_convert_output = str(tmp_path / "lusstr_output.txt")
    obs_all_output = str(tmp_path / "lusstr_output/")
    arglist = ["config", "-w", str(tmp_path), "--input", str(inputfile)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    snakemake_arglist = ["strs", command, "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(snakemake_arglist))
    assert os.path.exists(obs_format_output) is format_out
    assert os.path.exists(obs_convert_output) is convert_out
    assert os.path.exists(obs_all_output) is all_out
    assert filecmp.cmp(exp_output, obs_output) is True


@pytest.mark.parametrize(
    "sex",
    [True, False],
)
def test_marker_plots(sex, tmp_path):
    inputfile = data_file("UAS_bulk_input/Positive Control Sample Details Report 2315.xlsx")
    exp_output = str(tmp_path / "MarkerPlots/lusstr_output_Positive_Control_marker_plots.pdf")
    sex_exp = str(tmp_path / "MarkerPlots/lusstr_output_Positive_Control_sexchr_marker_plots.pdf")
    if sex:
        arglist = ["config", "-w", str(tmp_path), "--input", str(inputfile), "--sex"]
    else:
        arglist = ["config", "-w", str(tmp_path), "--input", str(inputfile)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    snakemake_arglist = ["strs", "convert", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(snakemake_arglist))
    assert os.path.exists(exp_output) is True
    if sex:
        assert os.path.exists(sex_exp) is True


def test_genemarker(tmp_path):
    input = data_file("genemarker/2800M_strresults_filtered.csv")
    arglist = [
        "config",
        "-w",
        str(tmp_path),
        "--out",
        "genemarker_test",
        "--input",
        str(input),
        "--sex",
        "--analysis-software",
        "genemarker",
        "--powerseq",
    ]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(arglist))
    snakemake_arglist = ["strs", "convert", "-w", str(tmp_path)]
    lusSTR.cli.main(lusSTR.cli.get_parser().parse_args(snakemake_arglist))
    for ext in [".csv", ".txt", "_flanks.txt", "_sexloci.csv", "_sexloci_flanks.txt"]:
        exp_output = data_file(f"genemarker/genemarker_test{ext}")
        print(exp_output)
        obs_output = str(tmp_path / f"genemarker_test{ext}")
        assert filecmp.cmp(exp_output, obs_output) is True
