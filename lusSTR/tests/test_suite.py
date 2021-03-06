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
from lusSTR.marker import STRMarkerObject
from lusSTR.repeat import reverse_complement
from lusSTR.tests import data_file
import re
from tempfile import NamedTemporaryFile


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
        testfullanno = data_file('2800M_full_anno_no_combined_reads.txt')
        arglist = [
            'annotate', inputfile, '-o', outfile.name, '--kit', 'forenseq', '--nocombine'
        ]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_no_combined_reads.txt'
        assert filecmp.cmp(testfullanno, outfile_name_output) is True


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


def test_FGA_short_seq():
    with NamedTemporaryFile(suffix='.txt') as outfile:
        input = data_file('test_FGA_short_seq.csv')
        arglist = [
            'annotate', input, '-o', outfile.name, '--kit', 'forenseq'
        ]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        with open(outfile.name, 'r') as fh:
            assert len(fh.readlines()) == 1


@pytest.mark.parametrize('locus, sequence, uas, kit, output', [
    (
        'CSF1PO', 'CTTCCTATCTATCTATCTATCTAATCTATCTATCTT', False, 'forenseq',
        'Possible indel or partial sequence'
    ),
    (
        'DYS393', 'AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATGTATGTCTTTTCTATGAGACATACC',
        False, 'powerseq',
        'UAS region indicates entire sequence; Possible indel or partial sequence'
    )
])
def test_indel_flag(locus, sequence, uas, kit, output):
    marker = STRMarkerObject(
        locus, sequence, uas=uas, kit=kit
    )
    assert marker.indel_flag == output


def test_powerseq_flanking_anno():
    with NamedTemporaryFile(suffix='.txt') as outfile:
        input = data_file('powerseq_flanking_anno_test.csv')
        test_powerseq = data_file('powerseq_flanking_anno_test_flanks_anno.txt')
        arglist = [
            'annotate', input, '-o', outfile.name, '--kit', 'powerseq'
        ]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_flanks_anno.txt'
        assert filecmp.cmp(test_powerseq, outfile_name_output) is True


def test_annotate_uas_sexloci():
    with NamedTemporaryFile() as outfile:
        os.unlink(outfile.name)
        inputfile = data_file('testformat_uas.csv')
        testanno = data_file('testformat_uas_sexloci.txt')
        arglist = [
            'annotate', inputfile, '-o', outfile.name, '--kit', 'forenseq', '--uas',
            '--include-sex'
            ]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_sexloci.txt'
        assert filecmp.cmp(testanno, outfile_name_output) is True


@pytest.mark.parametrize('inputfile, testoutput, flank_output, kit', [
    (
        'testformat_sr.csv', 'testformat_sr_sexloci.txt', 'testformat_sr_sexloci_flanks_anno.txt',
        'forenseq'
    ),
    (
        'powerseq_flanking_anno_test.csv', 'powerseq_flanking_anno_test_sexloci.txt',
        'powerseq_flanking_anno_test_sexloci_flanks_anno.txt', 'powerseq'
    )
])
def test_annotate_sr_sexloci(inputfile, testoutput, flank_output, kit):
    with NamedTemporaryFile() as outfile:
        os.unlink(outfile.name)
        inputfile = data_file(inputfile)
        testanno = data_file(testoutput)
        flankanno = data_file(flank_output)
        arglist = [
            'annotate', inputfile, '-o', outfile.name, '--kit', kit, '--include-sex'
            ]
        args = lusSTR.cli.get_parser().parse_args(arglist)
        lusSTR.annot.main(args)
        outfile_name = os.path.splitext(outfile.name)[0]
        outfile_name_output = f'{outfile_name}_sexloci.txt'
        assert filecmp.cmp(testanno, outfile_name_output) is True
        flank_outfile = f'{outfile_name}_sexloci_flanks_anno.txt'
        assert filecmp.cmp(flankanno, flank_outfile) is True
