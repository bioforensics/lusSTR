#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import filecmp
import lusSTR
from lusSTR.tests import data_file
import pytest
import os
from shutil import copytree
from tempfile import NamedTemporaryFile


def test_format():
    UAStestfile = data_file('snps/Positive Control Sample Details Report 2315.xlsx')
    formatoutput = data_file('testformat.csv')
    with NamedTemporaryFile(suffix='.csv') as outfile:
        arglist = ['format', UAStestfile, '-o', outfile.name, '--uas']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        assert filecmp.cmp(formatoutput, outfile.name) is True


def test_format_stdout(capsys):
    UAStestfile = data_file('snps/Positive Control Sample Details Report 2315.xlsx')
    formatoutput = data_file('testformat.csv')
    arglist = ['format', UAStestfile, '--uas']
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.format.main(args)
    with open(formatoutput, 'r') as fh:
        exp_out = fh.read().strip()
    terminal = capsys.readouterr()
    obs_out = terminal.out.strip()
    assert obs_out == exp_out


@pytest.mark.parametrize('input, testoutput', [
    ('STRait_Razor_test_output', 'STRait_Razor_test_output.csv'),
    ('STRait_Razor_test_output/A001.txt', 'STRaitRazor_output_test_A001.csv')
])
def test_format_straitrazor(input, testoutput):
    with NamedTemporaryFile() as outfile:
        inputdb = data_file(input)
        testformat = data_file(testoutput)
        arglist = ['format', inputdb, '-o', outfile.name]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        assert filecmp.cmp(testformat, outfile.name) is True


def test_format_sexloci_uas():
    UAStestfile = data_file('snps/Positive Control Sample Details Report 2315.xlsx')
    formatoutput = data_file('testformat_uas_sexloci.csv')
    with NamedTemporaryFile(suffix='.csv') as outfile:
        arglist = ['format', UAStestfile, '-o', outfile.name, '--uas', '--include-sex']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_sexloci.csv'
        assert filecmp.cmp(formatoutput, outfile_name_output) is True


def test_format_sex_loci_straitrazor(tmp_path):
    inputdb = data_file('STRait_Razor_test_output')
    exp_out = data_file('testformat_sr_sexloci.csv')
    obs_out = str(tmp_path / 'sr.csv')
    obs_out_sex = str(tmp_path / 'sr_sexloci.csv')
    arglist = ['format', inputdb, '-o', obs_out, '--include-sex']
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.format.main(args)
    assert filecmp.cmp(exp_out, obs_out_sex) is True


def test_uas_directory_autosomal_only(tmp_path):
    inputdb = data_file('UAS_bulk_input')
    copydb = str(tmp_path / 'UAS_bulk_input')
    copytree(inputdb, copydb)
    bogusfile = os.path.join(copydb, 'bogusfile.txt')
    with open(bogusfile, 'w') as fh:
        pass
    exp_out_auto = data_file('UAS_bulk_test.csv')
    obs_out_auto = str(tmp_path / 'format_output.csv')
    arglist = ['format', '-o', obs_out_auto, '--uas', copydb]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.format.main(args)
    assert filecmp.cmp(exp_out_auto, obs_out_auto) is True


def test_uas_directory_with_xy(tmp_path):
    inputdb = data_file('UAS_bulk_input')
    exp_out_auto = data_file('UAS_bulk_test.csv')
    exp_out_sex = data_file('UAS_bulk_test_sexloci.csv')
    obs_out_auto = str(tmp_path / 'format_output.csv')
    obs_out_sex = str(tmp_path / 'format_output_sexloci.csv')
    arglist = ['format', '-o', obs_out_auto, '--uas', '--include-sex', inputdb]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.format.main(args)
    assert filecmp.cmp(exp_out_auto, obs_out_auto) is True
    assert filecmp.cmp(exp_out_sex, obs_out_sex) is True
