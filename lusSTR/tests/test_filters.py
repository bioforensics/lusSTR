#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import filecmp
import json
import lusSTR
from lusSTR.filter_settings import get_filter_metadata_file
from lusSTR.tests import data_file
import pytest
from tempfile import NamedTemporaryFile


with open(get_filter_metadata_file(), 'r') as fh:
    filter_marker_data = json.load(fh)


@pytest.mark.parametrize('filter, locus, total_reads, allele_reads, final_reads, pass_filt', [
    ('Analytical', 'VWA', 318, 19, 318, False),
    ('Analytical', 'VWA', 318, 240, 318, True),
    ('Detection', 'VWA', 600, 9, 591, False),
    ('Detection', 'VWA', 600, 123, 600, True)
])
def test_thresholds(filter, locus, total_reads, allele_reads, final_reads, pass_filt):
    metadata = filter_marker_data[locus]
    test_total_reads, test_passfilt = lusSTR.filter_settings.thresholds(
        filter, metadata, total_reads, allele_reads
    )
    assert test_total_reads == final_reads
    assert test_passfilt == pass_filt


@pytest.mark.parametrize('perc, perc_stut, reads, forward_threshold', [
    (0, 0.18, 100, 4),
    (0.15, 0.21, 100, 15)
])
def test_forward_stutter_threshold(perc, perc_stut, reads, forward_threshold):
    test_forward_thresh = lusSTR.filter_settings.forward_stut_thresh(perc, perc_stut, reads)
    assert test_forward_thresh == forward_threshold


@pytest.mark.parametrize(
    'all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, ref_reads, al1_ref_reads,'
    'al_reads, called_allele_type, stut_perc', [
        (None, 0.18, 0, 18, 100, None, 15, '-1_stutter', 0.15),
        (None, 0.18, 0, 18, 100, None, 20, 'real_allele', None),
        ('+1_stutter', 0.18, 0, 18, 100, 200, 20, '-1_stutter/+1_stutter', None),
        ('+1_stutter', 0.18, 0, 18, 100, 200, 30, 'real_allele', None),
        ('-2_stutter', 0.18, 0, 18, 100, 100, 30, '-1_stutter/-2_stutter', None),
        ('-2_stutter', 0.18, 0, 18, 100, 100, 40, 'real_allele', None)
    ]
)
def test_minus1stutter(
    all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, ref_reads,
    al1_ref_reads, al_reads, called_allele_type, stut_perc
):
    test_stutter_type, test_stut_perc = lusSTR.filter_settings.minus1_stutter(
        all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, ref_reads, al1_ref_reads,
        al_reads)
    assert test_stutter_type == called_allele_type
    assert test_stut_perc == stut_perc


@pytest.mark.parametrize(
    'all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, al1_ref_reads,'
    'ref_reads, al_reads, called_allele_type, stut_perc', [
        (None, 0.18, 0, 18, None, 100, 15, '-2_stutter', 0.15),
        (None, 0.18, 0, 18, None, 100, 20, 'real_allele', None),
        ('+1_stutter', 0.18, 0, 18, 100, 200, 20, '+1_stutter/-2_stutter', None),
        ('+1_stutter', 0.18, 0, 18, 100, 200, 30, 'real_allele', None),
        ('-1_stutter', 0.18, 0, 18, 100, 100, 30, '-1_stutter/-2_stutter', None),
        ('-1_stutter', 0.18, 0, 18, 100, 100, 40, 'real_allele', None),
    ]
)
def test_minus2stutter(
    all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, ref_reads,
    al1_ref_reads, al_reads, called_allele_type, stut_perc
):
    test_stutter_type, test_stut_perc = lusSTR.filter_settings.minus2_stutter(
        all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, ref_reads, al1_ref_reads,
        al_reads)
    assert test_stutter_type == called_allele_type
    assert test_stut_perc == stut_perc


@pytest.mark.parametrize(
    'all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads,'
    'al_reads, called_allele_type, stut_perc', [
        (None, 0.18, 0, 100, None, 3, '+1_stutter', 0.03),
        (None, 0.18, 0, 100, None, 20, 'real_allele', None),
        ('-1_stutter', 0.18, 0, 100, 200, 3, '-1_stutter/+1_stutter', None),
        ('-1_stutter', 0.18, 0, 100, 200, 50, 'real_allele', None),
        ('-2_stutter', 0.18, 0, 100, 100, 3, '+1_stutter/-2_stutter', None),
        ('-2_stutter', 0.18, 0, 100, 100, 40, 'real_allele', None)
    ]
)
def test_plus1stutter(
    all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads,
    al_reads, called_allele_type, stut_perc
):
    test_stutter_type, test_stut_perc = lusSTR.filter_settings.plus1_stutter(
        all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads,
        al_reads)
    assert test_stutter_type == called_allele_type
    assert test_stut_perc == stut_perc


def test_EFMoutput_format(tmp_path):
    input_file = data_file('test_stutter.txt')
    exp_out = data_file('RU_stutter_test/test_filtering_EFMoutput.csv')
    exp_info_out = data_file('RU_stutter_test/test_filtering_EFMoutput_sequence_info.csv')
    obs_out = str(tmp_path / 'test_output.csv')
    obs_info_out = str(tmp_path / 'test_output_sequence_info.csv')
    arglist = ['filter', '-o', obs_out, '--output-type', 'efm', '--info', input_file]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True
    assert filecmp.cmp(exp_info_out, obs_info_out) is True


def test_STRmixoutput_format(tmp_path):
    input_file = data_file('test_stutter.txt')
    exp_out = data_file('RU_stutter_test/Sample1.csv')
    exp_info_out = data_file('RU_stutter_test/STRmix_Files_sequence_info.csv')
    obs_outdir = str(tmp_path / 'STRmix_Files')
    obs_out = str(tmp_path / 'STRmix_Files/Sample1.csv')
    obs_info_out = f'{obs_outdir}/STRmix_Files_sequence_info.csv'
    arglist = ['filter', '-o', obs_outdir, '--output-type', 'strmix', '--info', input_file]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True
    assert filecmp.cmp(exp_info_out, obs_info_out) is True


def test_stdout(capsys):
    input_file = data_file('test_stutter.txt')
    output = data_file('RU_stutter_test/test_filtering_EFMoutput.csv')
    arglist = ['filter', '--output-type', 'efm', input_file]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    with open(output, 'r') as fh:
        exp_out = fh.read().strip()
    terminal = capsys.readouterr()
    obs_out = terminal.out.strip()
    assert obs_out == exp_out


def test_flags(tmp_path):
    input_file = data_file('test_stutter.txt')
    exp_out = data_file('RU_stutter_test/Flagged_Loci.csv')
    obs_outdir = str(tmp_path / 'RU_stutter_test')
    obs_out = str(tmp_path / 'RU_stutter_test/Flagged_Loci.csv')
    arglist = ['filter', '-o', obs_outdir, '--output-type', 'strmix', '--info', input_file]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True


def test_efm_reference(tmp_path):
    input_file = data_file('test_references.txt')
    exp_out = data_file('RU_stutter_test/EFM_test_reference.csv')
    obs_out = str(tmp_path / 'test_output.csv')
    obs_efm_out = str(tmp_path / 'test_output_reference.csv')
    arglist = [
        'filter', '-o', obs_out, '--output-type', 'efm',
        '--profile-type', 'reference', input_file
    ]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_efm_out) is True


def test_strmix_reference(tmp_path):
    input_file = data_file('test_references.txt')
    exp_out = data_file('RU_stutter_test/Positive_Control_reference.csv')
    obs_out = str(tmp_path / 'Positive_Control_reference.csv')
    arglist = [
        'filter', '-o', str(tmp_path), '--output-type', 'strmix',
        '--profile-type', 'reference', input_file
    ]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.filter.main(args)
    assert filecmp.cmp(exp_out, obs_out) is True
