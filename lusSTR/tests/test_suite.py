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


@pytest.mark.parametrize('sequence, repeat_list, output', [
    (
        'AGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACAGACAGAT',
        ['AGAT', 'AGAC'], 'AGAC [AGAT]11 [AGAC]6 AGAT'
    ),
    (
        'TAGATAGATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGG',
        ['TCTA', 'CATA', 'TCTG', 'CACA', 'CCTA'],
        'TAGATAGATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGG'
    )
])
def test_collapse_all_repeats(sequence, repeat_list, output):
    final_output = lusSTR.annot.collapse_all_repeats(sequence, repeat_list)
    assert final_output == output


def test_split_by_n():
    sequence = 'AGGTAGGTAGGTCGAACGAATTGG'
    blocks = list(lusSTR.annot.split_by_n(sequence, n=4))
    assert blocks == [
        'AGGT', 'AGGT', 'AGGT', 'CGAA', 'CGAA', 'TTGG'
    ]


def test_sequence_to_bracketed_form():
    sequence = (
        'TCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCTATCTA'
    )
    repeats = ['TCTA', 'TCTG']
    final_output = lusSTR.annot.sequence_to_bracketed_form(sequence, 6, repeats)
    assert final_output == '[TCTA]3 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11'


def test_extract():
    s = '[ATCT]3 ATGT [ATCT]12'
    repeat = 'ATCT'
    final_output = lusSTR.annot.extract(s, repeat)
    assert str(final_output) == '12'


def test_split_sequence_into_two_strings():
    sequence = 'TAGATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGGTGTGTGTGTGTG'
    reverse_comp_sequence = lusSTR.annot.rev_complement_anno(sequence)
    repeat_for_split = 'CACA'
    seq1, seq2 = lusSTR.annot.split_sequence_into_two_strings(reverse_comp_sequence,
                                                              repeat_for_split)
    assert seq1 == 'CACACACACACA'
    assert seq2 == 'CCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATCTA'


def test_rev_complement_anno():
    sequence = 'TAGATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGGTGTGTGTGTGTG'
    final_output = lusSTR.annot.rev_complement_anno(sequence)
    assert final_output == (
        'CACACACACACACCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATCTA'
    )


def test_rev_comp_uas_output_bracket():
    foward_strand = '[AGGT]3 [CGAA]2 TTGG'
    rev_comp_bracket = lusSTR.annot.rev_comp_uas_output_bracket(foward_strand, 4)
    assert rev_comp_bracket == 'CCAA [TTCG]2 [ACCT]3'


def test_collapse_repeats_by_length():
    sequence = 'TCTATCTATCTATCTATCTATCTATCTATATATCTATCTATCTATCTA'
    assert lusSTR.annot.collapse_repeats_by_length(sequence, 4) == '[TCTA]7 TATA [TCTA]4'


@pytest.mark.parametrize('sequence, bracket_form', [
    (
        'TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCTATCTATATCTA',
        '[TCTA]5 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]10 TA TCTA'
    ),
    (
        'TCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCTAACTATCTA',
        '[TCTA]4 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]9 ACTA TCTA'
    ),
    (
        'TCTATCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATC'
        'TATCTATCTATCTATCTATCTATCTATCTA',
        '[TCTA]6 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]8'
    ),
    (
       'TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCT'
       'ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAACTA',
       '[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11 TA ACTA'
    ),
    (
       'TCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCT'
       'ATCTATCTATCTATCTATCTATCTATCTAACTATCTA',
       '[TCTA]4 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]9 ACTA TCTA'
    )
])
def test_D21_bracket(sequence, bracket_form):
    repeats = ['TCTA', 'TCTG']
    final_output = lusSTR.annot.D21_bracket(sequence, 6, repeats)
    assert final_output == bracket_form


def test_D19_annotation():
    sequence = (
        'AAGGAAAAGGTAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAGAGAAGAAGAAAGAGAG'
    )
    repeats = ['TCTA', 'TCTG']
    repeat_for_split = 'CCTT'
    reverse_comp_sequence = lusSTR.annot.rev_complement_anno(sequence)
    final_output = lusSTR.annot.D19_annotation(reverse_comp_sequence, repeats, repeat_for_split)
    assert final_output == 'CT CTCT TTCT TCTT CTCT [CCTT]14 CCTA CCTT TT CCTT'


def test_D1_annotation():
    sequence = 'TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGTGTATGTG'
    reverse_comp_sequence = lusSTR.annot.rev_complement_anno(sequence)
    repeats = ['TCTA', 'CATA', 'TCTG', 'CACA', 'CCTA']
    repeat_for_split = 'CACA'
    final_output = lusSTR.annot.D1_annotation(reverse_comp_sequence, repeats, repeat_for_split)
    print(final_output)
    assert final_output == 'CA CATA CACA [TCTA]11'


@pytest.mark.parametrize('sequence, bracket_form', [
    (
        'AAAAGAAAGAAAAGAAAAGAAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGA',
        'AAAAG [AAAGA]3 A [AAAGA]11'
    ),
    ('GAAAAGAAAAGAAAAGA', 'GA [AAAGA]3'),
    ('AAAAGAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAAAAAAGAAAAGA', 'AAAAG [AAAGA]7 AAAAA [AAAGA]2'),
    ('AAAAGAAAAAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGA', 'AAAAG AAAAA [AAAGA]8')
])
def test_PentaD_annotation(sequence, bracket_form):
    repeats = ['AAAGA']
    no_of_repeat_bases = 5
    final_output = lusSTR.annot.PentaD_annotation(sequence, no_of_repeat_bases, repeats)
    assert final_output == bracket_form


def test_FGA_anno():
    sequence = (
        'TTTCTTTCTTTCTTTCTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTGTCTGTCTGTCTTTCTTTCTTTCTTTCTT'
        'TCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTCCTTCCTTCCTTTCTTTCTTTCTCCTTCCTTCCTTCCTTCC'
    )
    repeats = ['AAAG', 'GAAA', 'GAAG', 'ACAG', 'AAAA']
    reverse_comp_sequence = lusSTR.annot.rev_complement_anno(sequence)
    final_output = lusSTR.annot.FGA_anno(reverse_comp_sequence, repeats)
    assert final_output == '[GGAA]4 GGAG [AAAG]3 [GAAG]3 [AAAG]15 [ACAG]3 [AAAG]9 AA AAAA [GAAA]4'


@pytest.mark.parametrize('sequence, n, n_sub_out, allele', [
    (
        'TCTATCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA', 4, 0, 17
    ),
    (
        'AGATTAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACAGACAGAT', 4,
        0, '21.1'
    ),
    (
        'TAGATAGATAGATAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGGTGTGTGTGTG', 4, 10, '14.3'
    ),
    (
        'AAAAGAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAA'
        'AGAAAAGA', 5, 5, 17
    )
])
def test_traditional_str_allele(sequence, n, n_sub_out, allele):
    assert lusSTR.annot.traditional_str_allele(sequence, n, n_sub_out) == allele


@pytest.mark.parametrize(
    'forward_bracket, lus, sec, tert, locus, lus_allele, sec_allele, tert_allele, str_allele', [
        (
            '[TCTA]5 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]13 TA TCAA', 'TCTA',
            'TCTA', 'TCTG', 'D21S11', '13', '5', '5', '31.2'
        ),
        ('[TATC]6 TGTC [TATC]4', 'TATC', 'TGTC', '', 'D7S820', '6', '1', '0', 11),
        ('TAGA TCGA [TAGA]10', 'TAGA', '', '', 'D9S1122', '10', '', '', 12),
        ('[GATA]5 GACA [GATA]4', 'GATA', 'GACA', '', 'D16S539', '5', '1', '', 10),
        (
            '[TATC]8 [AATC]2 [ATCT]3 TTCT GTCT GTC', 'TATC', 'ATCT', 'GTCT', 'D13S317', '8', '3',
            '1', '15.3'
        )
    ])
def test_lus_anno(forward_bracket, lus, sec, tert, locus, lus_allele, sec_allele, tert_allele,
                  str_allele):
    lus_out, sec_out, tert_out = lusSTR.annot.lus_anno(
        forward_bracket, lus, sec, tert, locus, str_allele
    )
    assert str(lus_out) == lus_allele
    assert str(sec_out) == sec_allele
    assert str(tert_out) == tert_allele


def test_THO1():
    sequence = 'AATGAATGAATGAATGAATGATGATGAATGAATGAATG'
    repeats = ['AATG']
    final_output = lusSTR.annot.TH01_annotation(sequence, repeats)
    assert final_output == '[AATG]5 ATG ATG [AATG]3'


@pytest.mark.parametrize('forward_bracket, lus_allele, sec_allele, tert_allele', [
    ('[TCTA]5 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]13 TA TCAA', '13', '5', '5'),
    ('[TCTA]6 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]8', '8', '6', '5'),
    ('[TCTA]4 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]12', '12', '4', '5')
])
def test_D21_anno(forward_bracket, lus_allele, sec_allele, tert_allele):
    lus = 'TCTA'
    sec = 'TCTA'
    tert = 'TCTG'
    final_lus, final_sec, final_tert = lusSTR.annot.lus_D21_anno(forward_bracket, lus, sec, tert)
    assert str(final_lus) == lus_allele
    assert str(final_sec) == sec_allele
    assert str(final_tert) == tert_allele


def test_D21_lus_sec():
    sequence = '[TCTA]4 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]10'
    repeat = 'TCTA'
    tert = 'TCTG'
    lus_out, sec_out = lusSTR.annot.D21_lus_sec(sequence, repeat, tert)
    assert str(lus_out) == '10'
    assert str(sec_out) == '4'


@pytest.mark.parametrize('sequence, uas_seq, front, back', [
    (
        'CTATGCATCTATCTATCTATCTATCTATCTATCTATCTATCTAATGGTTA',
        'ATCTATCTATCTATCTATCTATCTATCTATCTATCT',
        6,
        8,
    ),
    (
        'TCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTCCC',
        'TCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA',
        0,
        5,
    ),
])
def test_full_seq_to_uas(sequence, uas_seq, front, back):
    uas_sequence = lusSTR.annot.full_seq_to_uas(sequence, front, back)
    assert uas_sequence == uas_seq


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
