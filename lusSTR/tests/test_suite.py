#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import filecmp
import os
import pandas as pd
import pytest
import lusSTR
from lusSTR.repeat import reverse_complement
from lusSTR.tests import data_file
import re
from tempfile import NamedTemporaryFile


def test_format():
    UAStestfile = data_file('UAS_Sample_Details_Report_test.xlsx')
    formatoutput = data_file('testformat.csv')
    with NamedTemporaryFile(suffix='.csv') as outfile:
        arglist = ['format', UAStestfile, '-o', outfile.name, '--uas']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        assert filecmp.cmp(formatoutput, outfile.name) is True


def test_split_sequence_into_two_strings():
    sequence = 'TAGATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGGTGTGTGTGTGTG'
    reverse_comp_sequence = reverse_complement(sequence)
    repeat_for_split = 'CACA'
    seq1, seq2 = lusSTR.annot.split_sequence_into_two_strings(reverse_comp_sequence,
                                                              repeat_for_split)
    assert seq1 == 'CACACACACACA'
    assert seq2 == 'CCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATCTA'


def test_annotate_uas():
    with NamedTemporaryFile() as outfile:
        os.unlink(outfile.name)
        inputfile = data_file('2800M_formatted_uas.csv')
        testanno = data_file('2800M_uas_anno.txt')
        arglist = ['annotate', inputfile, '-o', outfile.name, '--kit', 'forenseq', '--uas']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        assert filecmp.cmp(testanno, outfile.name) is True


def test_annotate_full_nocombine():
    with NamedTemporaryFile() as outfile:
        inputfile = data_file('2800M_formatted_full.csv')
        testfullanno = data_file('2800M_full_anno.txt')
        arglist = [
            'annotate', inputfile, '-o', outfile.name, '--kit', 'forenseq', '--nocombine'
        ]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_no_combined_reads.txt'
        import subprocess
        subprocess.check_call(['cp', outfile_name_output, 'FLARFY'])
        assert filecmp.cmp(testfullanno, outfile_name_output) is True


def test_format_straitrazor():
    with NamedTemporaryFile() as outfile:
        inputdb = data_file('STRait_Razor_test_output/')
        testformat = data_file('STRait_Razor_test_output.csv')
        arglist = ['format', inputdb, '-o', outfile.name]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.format.main(args)
        assert filecmp.cmp(testformat, outfile.name) is True


def test_flank_anno():
    with NamedTemporaryFile(suffix='.txt') as outfile:
        inputfile = data_file('Flanks_testing_file.csv')
        testflanks = data_file('testflanks_flanks_anno.txt')
        arglist = ['annotate', inputfile, '-o', outfile.name, '--kit', 'forenseq']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_flanks_anno.txt'
        assert filecmp.cmp(testflanks, outfile_name_output) is True


def test_annotate_combine():
    with NamedTemporaryFile() as outfile:
        inputfile = data_file('Flanks_testing_file.csv')
        arglist = ['annotate', inputfile, '-o', outfile.name, '--kit', 'forenseq']
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        with open(outfile.name, 'r') as fh:
            assert len(fh.readlines()) == 952
