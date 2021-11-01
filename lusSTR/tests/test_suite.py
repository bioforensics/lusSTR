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


@pytest.mark.parametrize('infile, len_sum, len_uncom, xy_len_sum, xy_len_uncom, kit', [
    ('testformat_sr.csv', 897, 913, 9701, 16108, 'forenseq'),
    ('powerseq.csv', 353, 441, 256, 303, 'powerseq')
])
def test_annotate_full_nocombine(infile, len_sum, len_uncom, xy_len_sum, xy_len_uncom, kit):
    inputfile = data_file(infile)
    with NamedTemporaryFile() as outfile:
        arglist_nocomb = [
            'annotate', inputfile, '-o', outfile.name, '--kit', kit, '--nocombine', '--include-sex'
        ]
        args_nocomb = lusSTR.cli.get_parser().parse_args(arglist_nocomb)
        lusSTR.annot.main(args_nocomb)
        outfile_name = os.path.splitext(outfile.name)[0]
        with open(f'{outfile_name}_no_combined_reads.txt', 'r') as fh:
            assert len(fh.readlines()) == len_uncom
        with open(f'{outfile_name}_sexloci_no_combined_reads.txt', 'r') as fh:
            assert len(fh.readlines()) == xy_len_uncom
        arglist_comb = ['annotate', inputfile, '-o', outfile.name, '--kit', kit, '--include-sex']
        args_comb = lusSTR.cli.get_parser().parse_args(arglist_comb)
        lusSTR.annot.main(args_comb)
        with open(outfile.name, 'r') as fh:
            assert len(fh.readlines()) == len_sum
        with open(f'{outfile_name}_sexloci.txt', 'r') as fh:
            assert len(fh.readlines()) == xy_len_sum


def test_flank_anno():
    with NamedTemporaryFile(suffix='.txt') as outfile:
        inputfile = data_file('test_flank.csv')
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
            assert len(fh.readlines()) == 0


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


@pytest.mark.parametrize('flag', 'sex', [
    ('', ''),
    ('--include-sex', '_sexloci')
])
def separate_output(tmp_path, flag, sex):
    inputfile = data_file('UAS_bulk_test.csv')
    outputfile = str(tmp_path / 'UAS_bulk_test.txt')
    arglist = ['annotate', inputfile, '-o', outputfile, '--separate', flag]
    args = lusSTR.cli.get_parser().parse_args(arglist)
    lusSTR.annot.main(args)
    assert os.file.exists(
        f'{tmp_path}/Separated_lusstr_Files/UAS_bulk_test/Positive_Control{sex}.txt'
    )
    assert os.file.exists(
        f'{tmp_path}/Separated_lusstr_Files/UAS_bulk_test/Positive_Control2{sex}.txt'
    )
