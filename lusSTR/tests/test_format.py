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
import os
from lusSTR.tests import data_file
from tempfile import NamedTemporaryFile


def test_format():
    UAStestfile = data_file('UAS_Sample_Details_Report_test.xlsx')
    formatoutput = data_file('testformat.csv')
    with NamedTemporaryFile(suffix='.csv') as outfile:
        arglist = ['format', UAStestfile, '-o', outfile.name, '--uas']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        assert filecmp.cmp(formatoutput, outfile.name) is True


def test_format_stdout(capsys):
    UAStestfile = data_file('UAS_Sample_Details_Report_test.xlsx')
    formatoutput = data_file('testformat.csv')
    arglist = ['format', UAStestfile, '--uas']
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.format.main(args)
    with open(formatoutput, 'r') as fh:
        exp_out = fh.read().strip()
    terminal = capsys.readouterr()
    obs_out = terminal.out.strip()
    assert obs_out == exp_out


def test_format_straitrazor():
    with NamedTemporaryFile() as outfile:
        inputdb = data_file('STRait_Razor_test_output/')
        testformat = data_file('STRait_Razor_test_output.csv')
        arglist = ['format', inputdb, '-o', outfile.name]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        assert filecmp.cmp(testformat, outfile.name) is True


def test_format_sexloci_uas():
    UAStestfile = data_file('UAS_Sample_Details_Report_test.xlsx')
    formatoutput = data_file('testformat_uas_sexloci.csv')
    with NamedTemporaryFile(suffix='.csv') as outfile:
        arglist = ['format', UAStestfile, '-o', outfile.name, '--uas', '--include-sex']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_sexloci.csv'
        assert filecmp.cmp(formatoutput, outfile_name_output) is True


def test_format_sex_loci_straitrazor():
    with NamedTemporaryFile() as outfile:
        inputdb = data_file('STRait_Razor_test_output/')
        testformat = data_file('testformat_sr_sexloci.csv')
        arglist = ['format', inputdb, '-o', outfile.name, '--include-sex']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_sexloci.csv'
        assert filecmp.cmp(testformat, outfile_name_output) is True


def test_uas_directory():
    with NamedTemporaryFile() as outfile:
        inputdb = data_file('UAS_bulk_input/')
        testformat_auto = data_file('UAS_bulk_test.csv')
        testformat_sex = data_file('UAS_bulk_test_sexloci.csv')
        arglist = ['format', inputdb, '-o', outfile.name, '--uas', '--include-sex']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        assert filecmp.cmp(testformat_auto, outfile.name) is True
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_sexloci.csv'
        assert filecmp.cmp(testformat_sex, outfile_name_output) is True
