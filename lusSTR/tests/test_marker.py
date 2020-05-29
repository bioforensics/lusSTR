#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import pytest
import lusSTR
from lusSTR.marker import STRMarkerObject
from lusSTR.repeat import reverse_complement as rc


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
    marker = STRMarkerObject('D21S11', sequence, uas=True)
    assert marker.annotation == bracket_form


def test_D19_annotation():
    sequence = (
        'AAGGAAAAGGTAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAGAGAAGAAGAAAGAGAG'
    )
    marker = STRMarkerObject('D19S433', sequence, uas=True)
    assert marker.annotation == 'CT CTCT TTCT TCTT CTCT [CCTT]14 CCTA CCTT TT CCTT'


def test_D1_annotation():
    marker = STRMarkerObject('D1S1656', 'TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGTGTATGTG', uas=True)
    assert marker.annotation == 'CA CATA CACA [TCTA]11'


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
    marker = STRMarkerObject('PentaD', sequence, uas=True)
    assert marker.annotation == bracket_form


def test_FGA_anno():
    sequence = (
        'TTTCTTTCTTTCTTTCTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTGTCTGTCTGTCTTTCTTTCTTTCTTTCTT'
        'TCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTCCTTCCTTCCTTTCTTTCTTTCTCCTTCCTTCCTTCCTTCC'
    )
    annotation = '[GGAA]4 GGAG [AAAG]3 [GAAG]3 [AAAG]15 [ACAG]3 [AAAG]9 AA AAAA [GAAA]4'
    marker = STRMarkerObject('FGA', sequence, uas=True)
    print(annotation)
    print(marker.annotation)
    assert marker.annotation == annotation


def test_THO1():
    marker = STRMarkerObject('TH01', 'AATGAATGAATGAATGAATGATGATGAATGAATGAATG', uas=True)
    assert marker.annotation == '[AATG]5 ATG ATG [AATG]3'


@pytest.mark.parametrize('sequence, lus_allele, sec_allele, tert_allele', [
    (
        'TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCAA', '13', '5', '5'
    ),
    (
        'TCTATCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATC'
        'TATCTATCTATCTATCTATCTATCTATCTA', '8', '6', '5'
    ),
    (
        'TCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCTATCTATCTA', '12', '4', '5'
    ),
])
def test_D21_anno(sequence, lus_allele, sec_allele, tert_allele):
    marker = STRMarkerObject('D21S11', sequence, uas=True)
    lus, sec, tert = marker.designation
    assert str(lus) == lus_allele
    assert str(sec) == sec_allele
    assert str(tert) == tert_allele


def test_D21_lus_sec():
    sequence = (
        'TCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCTATCTA'
    )
    marker = STRMarkerObject('D21S11', sequence, uas=True)
    lus, sec, tert = marker.designation
    assert str(lus) == '10'
    assert str(sec) == '4'
    assert str(tert) == '6'


@pytest.mark.parametrize('locus, sequence, forward_bracket, lus, sec, tert', [
        (
            'D21S11', 'TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATC'
            'TATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCAA',
            '[TCTA]5 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]13 TA TCAA',
            '13', '5', '5'
        ),
        (
            'D7S820', rc('TATCTATCTATCTATCTATCTATCTGTCTATCTATCTATCTATC'), '[TATC]6 TGTC [TATC]4',
            '6', '1', '0'
        ),
        (
            'D9S1122', rc('TAGATCGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA'),
            'TAGA TCGA [TAGA]10', '10', 'None', 'None'
        ),
        (
            'D16S539', rc('GATAGATAGATAGATAGATAGACAGATAGATAGATAGATA'), '[GATA]5 GACA [GATA]4',
            '5', '1', 'None'
        ),
        (
            'D13S317', 'TATCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTC',
            '[TATC]8 [AATC]2 [ATCT]3 TTCT GTCT GTC', '8', '3', '1'
        )
    ])
def test_annotation_and_lus(locus, sequence, forward_bracket, lus, sec, tert):
    marker = STRMarkerObject(locus, sequence, uas=True)
    assert marker.annotation == forward_bracket
    lus_out, sec_out, tert_out = marker.designation
    assert str(lus_out) == lus
    assert str(sec_out) == sec
    assert str(tert_out) == tert


def test_strobj_CSF1PO():
    marker = STRMarkerObject(
        'CSF1PO', 'CTTCCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAATCTATCTATCTT',
        uas=False, kit='forenseq'
    )
    assert marker.uas_sequence == 'AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT'
    assert marker.forward_sequence == 'ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT'
    assert marker.annotation == '[ATCT]12'
    assert marker.annotation_reverse == '[AGAT]12'
    assert marker.canonical == 12
    assert marker.designation == ('12', '0', None)


def test_strobj_D10S1248():
    marker = STRMarkerObject(
        'D10S1248', 'TTGAACAAATGAGTGAGTGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAATGAAGA'
        'CAATACAACCAGAGTT', uas=False, kit='forenseq'
    )
    assert marker.uas_sequence == 'GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA'
    assert marker.forward_sequence == 'GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA'
    assert marker.annotation == '[GGAA]13'
    assert marker.annotation_reverse == '[GGAA]13'
    assert marker.canonical == 13
    assert marker.designation == ('13', None, None)
