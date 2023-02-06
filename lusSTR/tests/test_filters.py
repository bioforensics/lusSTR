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
import json
import lusSTR
from lusSTR.filter_settings import get_filter_metadata_file
from lusSTR.tests import data_file
import pytest
from tempfile import NamedTemporaryFile


with open(get_filter_metadata_file(), "r") as fh:
    filter_marker_data = json.load(fh)


@pytest.mark.parametrize(
    "filter, locus, total_reads, allele_reads, final_reads, pass_filt",
    [
        ("Analytical", "VWA", 318, 19, 318, False),
        ("Analytical", "VWA", 318, 240, 318, True),
        ("Detection", "VWA", 600, 9, 591, False),
        ("Detection", "VWA", 600, 123, 600, True),
    ],
)
def test_thresholds(filter, locus, total_reads, allele_reads, final_reads, pass_filt):
    metadata = filter_marker_data[locus]
    test_total_reads, test_passfilt = lusSTR.filter_settings.thresholds(
        filter, metadata, total_reads, allele_reads
    )
    assert test_total_reads == final_reads
    assert test_passfilt == pass_filt


@pytest.mark.parametrize(
    "perc, perc_stut, reads, forward_threshold", [(0, 0.18, 100, 4), (0.15, 0.21, 100, 15)]
)
def test_forward_stutter_threshold(perc, perc_stut, reads, forward_threshold):
    test_forward_thresh = lusSTR.filter_settings.forward_stut_thresh(perc, perc_stut, reads)
    assert test_forward_thresh == forward_threshold


@pytest.mark.parametrize(
    "all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, ref_reads, al1_ref_reads,"
    "al_reads, called_allele_type, stut_perc",
    [
        (None, 0.18, 0, 18, 100, None, 15, "-1_stutter", 0.15),
        ("real_allele", 0.18, 0, 18, 100, None, 20, "real_allele", None),
        (None, 0.18, 0, 18, 100, None, 20, None, None),
        ("+1_stutter", 0.18, 0, 18, 100, 200, 20, "-1_stutter/+1_stutter", None),
        ("+1_stutter", 0.18, 0, 18, 100, 200, 30, "real_allele", None),
        ("-2_stutter", 0.18, 0, 18, 100, 100, 30, "-1_stutter/-2_stutter", None),
        ("-2_stutter", 0.18, 0, 18, 100, 100, 40, "real_allele", None),
    ],
)
def test_minus1stutter(
    all_type,
    stutter_thresh,
    forward_thresh,
    stutter_thresh_reads,
    ref_reads,
    al1_ref_reads,
    al_reads,
    called_allele_type,
    stut_perc,
):
    test_stutter_type, test_stut_perc = lusSTR.filter_settings.minus1_stutter(
        all_type,
        stutter_thresh,
        forward_thresh,
        stutter_thresh_reads,
        ref_reads,
        al1_ref_reads,
        al_reads,
    )
    assert test_stutter_type == called_allele_type
    assert test_stut_perc == stut_perc


@pytest.mark.parametrize(
    "all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, al1_ref_reads,"
    "ref_reads, al_reads, called_allele_type, stut_perc",
    [
        (None, 0.18, 0, 18, None, 100, 15, "-2_stutter", 0.15),
        ("real_allele", 0.18, 0, 18, None, 100, 20, "real_allele", None),
        (None, 0.18, 0, 18, None, 100, 20, None, None),
        ("+1_stutter", 0.18, 0, 18, 100, 200, 20, "+1_stutter/-2_stutter", None),
        ("+1_stutter", 0.18, 0, 18, 100, 200, 30, "real_allele", None),
        ("-1_stutter", 0.18, 0, 18, 100, 100, 30, "-1_stutter/-2_stutter", None),
        ("-1_stutter", 0.18, 0, 18, 100, 100, 40, "real_allele", None),
    ],
)
def test_minus2stutter(
    all_type,
    stutter_thresh,
    forward_thresh,
    stutter_thresh_reads,
    ref_reads,
    al1_ref_reads,
    al_reads,
    called_allele_type,
    stut_perc,
):
    test_stutter_type, test_stut_perc = lusSTR.filter_settings.minus2_stutter(
        all_type,
        stutter_thresh,
        forward_thresh,
        stutter_thresh_reads,
        ref_reads,
        al1_ref_reads,
        al_reads,
    )
    assert test_stutter_type == called_allele_type
    assert test_stut_perc == stut_perc


@pytest.mark.parametrize(
    "all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads,"
    "al_reads, called_allele_type, stut_perc",
    [
        (None, 0.18, 0, 100, None, 3, "+1_stutter", 0.03),
        ("real_allele", 0.18, 0, 100, None, 20, "real_allele", None),
        (None, 0.18, 0, 100, None, 20, None, None),
        ("-1_stutter", 0.18, 0, 100, 200, 3, "-1_stutter/+1_stutter", None),
        ("-1_stutter", 0.18, 0, 100, 200, 50, "real_allele", None),
        ("-2_stutter", 0.18, 0, 100, 100, 3, "+1_stutter/-2_stutter", None),
        ("-2_stutter", 0.18, 0, 100, 100, 40, "real_allele", None),
    ],
)
def test_plus1stutter(
    all_type,
    stutter_thresh,
    forward_thresh,
    ref_reads,
    al1_ref_reads,
    al_reads,
    called_allele_type,
    stut_perc,
):
    test_stutter_type, test_stut_perc = lusSTR.filter_settings.plus1_stutter(
        all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads, al_reads
    )
    assert test_stutter_type == called_allele_type
    assert test_stut_perc == stut_perc


def test_EFMoutput_format(tmp_path):
    input_file = data_file("test_stutter.txt")
    exp_out = data_file("RU_stutter_test/test_filtering_EFMoutput.csv")
    exp_info_out = data_file("RU_stutter_test/test_filtering_EFMoutput_sequence_info.csv")
    obs_out = str(tmp_path / "test_output.csv")
    obs_info_out = str(tmp_path / "test_output_sequence_info.csv")
    arglist = ["filter", "-o", obs_out, "--output-type", "efm", "--info", input_file]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True
    assert filecmp.cmp(exp_info_out, obs_info_out) is True


@pytest.mark.parametrize(
    "outputdir, datatype", [("RU_stutter_test/", "ce"), ("NGS_stutter_test/", "ngs")]
)
def test_STRmixoutput_format(outputdir, datatype, tmp_path):
    input_file = data_file("test_stutter.txt")
    exp_out = data_file(f"{outputdir}Sample1_{datatype}.csv")
    exp_info_out = data_file(f"{outputdir}STRmix_Files_sequence_info.csv")
    obs_outdir = str(tmp_path / "STRmix_Files")
    obs_out = str(tmp_path / f"STRmix_Files/Sample1_{datatype}.csv")
    obs_info_out = str(tmp_path / f"STRmix_Files/sequence_info.csv")
    arglist = [
        "filter",
        "-o",
        obs_outdir,
        "--output-type",
        "strmix",
        "--info",
        "--data-type",
        datatype,
        input_file,
    ]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True
    assert filecmp.cmp(exp_info_out, obs_info_out) is True


def test_stdout(capsys):
    input_file = data_file("test_stutter.txt")
    output = data_file("RU_stutter_test/test_filtering_EFMoutput.csv")
    arglist = ["filter", "--output-type", "efm", input_file]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    with open(output, "r") as fh:
        exp_out = fh.read().strip()
    terminal = capsys.readouterr()
    obs_out = terminal.out.strip()
    assert obs_out == exp_out


def test_nofilters(tmp_path):
    input_file = data_file("test_stutter.txt")
    exp_out = data_file("NGS_stutter_test/Sample1_nofilter.csv")
    obs_out = str(tmp_path / "Sample1_ngs.csv")
    arglist = [
        "filter",
        "-o",
        str(tmp_path),
        "--output-type",
        "strmix",
        "--no-filters",
        "--data-type",
        "ngs",
        input_file,
    ]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True


def test_flags(tmp_path):
    input_file = data_file("test_stutter.txt")
    exp_out = data_file("RU_stutter_test/Flagged_Loci.csv")
    obs_outdir = str(tmp_path / "RU_stutter_test")
    obs_out = str(tmp_path / "RU_stutter_test/Flagged_Loci.csv")
    arglist = ["filter", "-o", obs_outdir, "--output-type", "strmix", "--info", input_file]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True


def test_efm_reference(tmp_path):
    input_file = data_file("test_references.txt")
    exp_out = data_file("RU_stutter_test/EFM_test_reference.csv")
    obs_out = str(tmp_path / "test_output.csv")
    obs_efm_out = str(tmp_path / "test_output_reference.csv")
    arglist = [
        "filter",
        "-o",
        obs_out,
        "--output-type",
        "efm",
        "--profile-type",
        "reference",
        input_file,
    ]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_efm_out) is True


def test_strmix_reference(tmp_path):
    input_file = data_file("test_references.txt")
    exp_out = data_file("RU_stutter_test/Positive_Control_reference.csv")
    obs_out = str(tmp_path / "Positive_Control_reference_ce.csv")
    arglist = [
        "filter",
        "-o",
        str(tmp_path),
        "--output-type",
        "strmix",
        "--profile-type",
        "reference",
        input_file,
    ]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True


def test_D7(tmp_path):
    input_file = data_file("test_D7.txt")
    exp_out = data_file("D7_microvariant_flagged.csv")
    obs_out = str(tmp_path / "Flagged_Loci.csv")
    arglist = ["filter", "-o", str(tmp_path), "--output-type", "strmix", "--info", input_file]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out)


def test_ngs_reference_error(capsys):
    input_file = data_file("test_stutter.txt")
    arglist = [
        "filter",
        "--output-type",
        "strmix",
        "--data-type",
        "ngs",
        "--profile-type",
        "reference",
        input_file,
    ]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    exp_error = "Cannot create reference file from ngs data. Abort!"
    with pytest.raises(ValueError, match=exp_error) as fnfe:
        fullpath = lusSTR.filter.main(args)


@pytest.mark.parametrize(
    "ref_bracket, quest_bracket, stutter, actual_call",
    [
        ("TCTA [TCTG]3 [TCTA]10 TCCA TCTA", "TCTA [TCTG]3 [TCTA]9 TCCA TCTA", -1, -1),
        ("TCTA [TCTG]3 [TCTA]10 TCCA TCTA", "TCTA [TCTG]4 [TCTA]8 TCCA TCTA", -1, None),
        ("[TAGA]4 TGA [TAGA]10 TAGG [TGTG]2 TG", "[TAGA]4 TGA [TAGA]11 TAGG [TGTG]2 TG", 1, 1),
        ("[TAGA]4 TGA [TAGA]10 TAGG [TGTG]2 TG", "[TAGA]6 TGA [TAGA]9 TAGG [TGTG]2 TG", 1, None),
        ("[TTTC]3 TTTT [CTTT]17 CTCC [TTCC]2", "[TTTC]3 TTTT [CTTT]15 CTCC [TTCC]2", -2, -2),
        ("[TTTC]3 TTTT [CTTT]17 CTCC [TTCC]2", "[TTTC]4 TTTT [CTTT]14 CTCC [TTCC]2", -2, None),
    ],
)
def test_ngs_stutter(ref_bracket, quest_bracket, stutter, actual_call):
    test_stutter = lusSTR.filter_settings.bracketed_stutter_id(ref_bracket, quest_bracket, stutter)
    assert test_stutter == actual_call
