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
    uas_sequence = (
        'AAGGAAAAGGTAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAGAGAAGAAGAAAGAGAG'
    )
    marker = STRMarkerObject('D19S433', uas_sequence, uas=True)
    assert marker.annotation == 'CT CTCT TTCT TCTT CTCT [CCTT]14 CCTA CCTT TT CCTT'


def test_D1_annotation():
    uas_sequence = 'TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGTGTATGTG'
    marker = STRMarkerObject('D1S1656', uas_sequence, uas=True)
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
    marker = STRMarkerObject('PENTA D', sequence, uas=True)
    assert marker.annotation == bracket_form


def test_FGA_anno():
    uas_sequence = (
        'TTTCTTTCTTTCTTTCTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTGTCTGTCTGTCTTTCTTTCTTTCTTTCTT'
        'TCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTCCTTCCTTCCTTTCTTTCTTTCTCCTTCCTTCCTTCCTTCC'
    )
    annotation = '[GGAA]4 GGAG [AAAG]3 [GAAG]3 [AAAG]15 [ACAG]3 [AAAG]9 AA AAAA [GAAA]4'
    marker = STRMarkerObject('FGA', uas_sequence, uas=True)
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
            'D7S820', 'GATAGATAGATAGATAGACAGATAGATAGATAGATAGATAGATA',
            '[TATC]6 TGTC [TATC]4', '6', '1', '0'
        ),
        (
            'D9S1122', 'TAGATCGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA',
            'TAGA TCGA [TAGA]10', '10', 'None', 'None'
        ),
        (
            'D16S539', 'GATAGATAGATAGATAGATAGACAGATAGATAGATAGATA', '[GATA]5 GACA [GATA]4',
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
    assert marker.annotation_uas == '[AGAT]12'
    assert marker.canonical == 12, ' '
    assert marker.designation == ('12', '0', None)
    assert marker.flank_5p == 'CT TCCT'


def test_strobj_D10S1248():
    marker = STRMarkerObject(
        'D10S1248', 'TTGAACAAATGAGTGAGTGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAATGAAGA'
        'CAATACAACCAGAGTT', uas=False, kit='forenseq'
    )
    assert marker.uas_sequence == 'GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA'
    assert marker.forward_sequence == 'GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA'
    assert marker.annotation == '[GGAA]13'
    assert marker.annotation_uas == '[GGAA]13'
    assert marker.canonical == 13, ' '
    assert marker.designation == ('13', None, None)


def test_strobj_D1S1656():
    marker = STRMarkerObject(
        'D1S1656', 'TTCAGAGAAATAGAATCACTAGGGAACCAAATATATATACATACAATTAAACACACACACACCTATCTATCTATCTAT'
        'CTATCTATCTATCTATCTATCTATCTATCTA', uas=False, kit='forenseq'
    )
    assert marker.uas_sequence == 'TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGGTGTGTGTGTG'
    assert marker.forward_sequence == 'CACACACACACCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA'
    assert marker.annotation == 'CA [CACA]2 CCTA [TCTA]11'
    assert marker.annotation_uas == '[TAGA]11 TAGG [TGTG]2 TG'
    assert marker.canonical == 12, ' '
    assert marker.designation == ('11', '1', '0')


def test_strobj_D5S818():
    marker = STRMarkerObject(
        'D5S818', 'TATTTATACCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTTCAAAAT',
        uas=False, kit='forenseq'
    )
    assert marker.uas_sequence == 'AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAG'
    assert marker.forward_sequence == 'CTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT'
    assert marker.annotation == 'CTCT [ATCT]12'
    assert marker.annotation_uas == '[AGAT]12 AGAG'
    assert marker.canonical == 12, ' '
    assert marker.designation == ('12',  None, None)


def test_strobj_D16S539():
    marker = STRMarkerObject(
        'D16S539', 'TCCTCTTCCCTAGATCAATACAGACAGACAGACAGGTGGATAGATAGATAGATTGATTGATAGATAGATAGATAGATA'
        'TCATTGAAAGACAAAACAGAGATGGATGATAGATAC',
        uas=False, kit='forenseq'
    )
    assert marker.uas_sequence == 'GATAGATAGATAGATTGATTGATAGATAGATAGATAGATA'
    assert marker.forward_sequence == 'GATAGATAGATAGATTGATTGATAGATAGATAGATAGATA'
    assert marker.flank_5p == 'TC CTCT T CCCT AGAT CAAT [ACAG]4 GTG'
    assert marker.flank_3p == 'TCAT TGAA AGAC AAA A CAGA [GATG]2 ATA GA T AC'
    assert marker.annotation == '[GATA]3 [GATT]2 [GATA]5'


def test_strobj_D7S820():
    marker = STRMarkerObject(
        'D7S820', 'TATTTAGTGAGATAAAAAAAAAACTATCAATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCGTT'
        'AGTTCGTTCTAAACTAT',
        uas=False, kit='forenseq'
    )
    assert marker.uas_sequence == 'GATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGATTGATAGTTTT'
    assert marker.forward_sequence == 'AAAACTATCAATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC'
    assert marker.flank_5p == 'T ATTT AGTG AGAT AAAAAA'
    assert marker.flank_3p == 'GTTA [GTTC]2 TAAA CTAT'
    assert marker.annotation == 'A AAAC TATC AATC TGTC [TATC]10'


def test_strobj_D3S1358():
    sequence = 'TCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTAACTAACTATCTATCTA'
    marker = STRMarkerObject('D3S1358', sequence, uas=True, kit='forenseq')
    assert marker.forward_sequence == sequence
    assert marker.uas_sequence == sequence
    assert marker.annotation == 'TCTA [TCTG]2 [TCTA]9 [ACTA]2 [TCTA]2'


def test_strobj_D19S433_newformat():
    marker = STRMarkerObject(
        'D19S433',
        'AATAAAAATCTTCTCTCTTTCTTCCTCTCTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTACCT'
        'TCCTTCCTTCAACAGAATCTTATTCTGTTGCCCAGGCTGGAGTCCAGTGTTACAATTATAGCT', uas=False,
        kit='forenseq'
    )
    assert marker.annotation == 'CT CTCT TTCT TCCT CTCT [CCTT]12 CCTA [CCTT]3'


def test_strobj_D21S11_newformat():
    marker = STRMarkerObject(
        'D21S11',
        'AAATATGTGAGTCAATTCCCCAAGTGAATTGCCTTCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCT'
        'ATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCT'
        'ACTATCTAT', uas=False, kit='forenseq'
    )
    assert marker.annotation == (
        '[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11 TA'
    )


def test_strobj_FGA_newformat():
    marker = STRMarkerObject(
        'FGA', 'CCAGCAAAAAAGAAAGAAAGAGAAAAAAGAAAGAAAGAAA', uas=False, kit='forenseq'
    )
    assert marker.annotation == 'AAAA [AGAA]3 A'


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_3p, kit', [
    (
        'TGATTTTTGCAGGTGTTCACTGCAAGCCATGCCTGGTTAAACTACTGTGCCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTT'
        'TTCTTTTCTTTTCTTTCTTTTTAAAACTT', '[CTTTT]10 CTTTC TTTT', '10', '10', None, None,
        'TAAAA CTT', 'forenseq'
    ),
    (
        'ATGCTCTGTGATTTTTGCAGGTGTTCACTGCAAGCCATGCCTGGTTAAACTACTGTGCCTTTTCTTTTCTTTTCTTTTCTTTTCTTTT'
        'CTTTTCTTTTCTTTTCTTTTCTTTCTTTTTAAAACTTTTTACTTCAGTAGAATTTTGGGG', '[CTTTT]10 CTTTC TTTT',
        '10', '10', None, None, 'TAAAA CTT TTTAC TTCAG TAGAA TTTTG GGG', 'powerseq'
    )
])
def test_strobj_DYS643(sequence, bracketed, conc, lus, sec, tert, flank_3p, kit):
    marker = STRMarkerObject('DYS643', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_3p == flank_3p


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_5p, kit', [
    (
        'ATCAATCAATGAATGGATAAAGAAAATGTGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATACATACAT'
        'AGATAGATACATACATAGATAGATAGATAGAGATTCTATGCAAAGTGAGAAGCCA',
        '[TAGA]12 [TACA]2 [TAGA]2 [TACA]2 [TAGA]4', '22', '12', None, None,
        'A [TCAA]2 TGAA TGGA TAAA GAAA ATGT GA', 'forenseq'
    ),
    (
        'CCAAATATCCATCAATCAATGAATGGATAAAGAAAATGTGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAC'
        'ATACATAGATAGATACATACATAGATAGATAGATAGAG', '[TAGA]11 [TACA]2 [TAGA]2 [TACA]2 [TAGA]4',
        '21', '11', None, None, 'CC AAAT ATCC A [TCAA]2 TGAA TGGA TAAA GAAA ATGT GA', 'powerseq'
    )
])
def test_strobj_DYS635(sequence, bracketed, conc, lus, sec, tert, flank_5p, kit):
    marker = STRMarkerObject('DYS635', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5p


def test_strobj_DYS612():
    marker = STRMarkerObject(
        'DYS612',
        'TTTCACACAGGTTCAGAGGTTTGCCTCCTCCTCCTCCTCTTTCTTCTTCTTCTCCTTCTTCTTCTTCTTCTTCTTCTT'
        'CTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTGTCACTTTTCCAAATTATTTTCTTTT', uas=False,
        kit='forenseq')
    assert marker.annotation == '[CCT]5 CTT [TCT]4 CCT [TCT]24'
    assert marker.canonical == 29
    assert marker.designation == (24, 4, None)
    assert marker.flank_3p == 'G TCA CTT TTC CAA [ATT]2 TTC TTT T'


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_3p, kit', [
    (
        'AAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGCCAAGACAAATACG'
        'CTTATTACTCCCATCTCCT', '[AAAG]17', '17', '17', None, None,
        'AAAA AGCC AAGA CAAA TACG CTTA TTAC TCCC ATCT CCT', 'forenseq'
    ),
    (
        'TAAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGCCAAGACA'
        'AATACGCTTATTACTCCCATCTCCTCCTTCATCTCCAGGAAATGAGACTG', '[AAAG]18', '18', '18', None, None,
        'AAAA AGCC AAGA CAAA TACG CTTA TTAC TCCC ATCT CCT CCTT CATC TCCA GGAA ATGA GACT G',
        'powerseq'
    )
])
def test_strobj_DYS576(sequence, bracketed, conc, lus, sec, tert, flank_3p, kit):
    marker = STRMarkerObject('DYS576', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_3p == flank_3p


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_5p, kit', [
    (
        'TAATAAGGTAGACATAGCAATTAGGTAGGTAAAGAGGAAGATGATAGATGATTAGAAAGATGATAGATAGATAGATAGATAGATAGAT'
        'AGATAGATAGATAGATAGATAGATAGAAAAAATCTACATAAACAAAATCACAAATGGAAAAGGGGACATTACCA', '[GATA]13',
        '13', '13', None, None,
        'TA ATAA GGTA GACA TAGC AATT [AGGT]2 AAAG AGGA AGAT GATA GATG ATTA GAAA GAT', 'forenseq'
    ),
    (
        'GATGTAAAGAACTATAAAAAGATTAATACAACAAAAATTTGGTAATCTGAAATAATAAGGTAGACATAGCAATTAGGTAGGTAAAGAG'
        'GAAGATGATAGATGATTAGAAAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAAAATCTA'
        'C', '[GATA]13', '13', '13', None, None, 'GATG TAAA GAAC TATA AAAA GATT AATA CAAC AAAA '
        'ATTT GGTA ATCT GAAA TA ATAA GGTA GACA TAGC AATT [AGGT]2 AAAG AGGA AGAT GATA GATG ATTA '
        'GAAA GAT', 'powerseq'
    )
])
def test_strobj_DYS549(sequence, bracketed, conc, lus, sec, tert, flank_5p, kit):
    marker = STRMarkerObject('DYS549', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5p


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_3p, kit', [
    (
        'TAACTATATAACTATGTATTATCTATCAATCTTCTACCTATCATCTTTCTAGCTAGCTATCATCTATCTATCTATCTATCTATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCATCTATCATCTTCTATTGTTT', '[TATC]13', '13', '13', None, None,
        'ATCT ATCA TCTT CTAT TGTT T', 'forenseq'
    ),
    (
        'CTACCTAATATTTATCTATATCATTCTAATTATGTCTCTTCTAACTATATAACTATGTATTATCTATCAATCTTCTACCTATCATCTT'
        'TCTAGCTAGCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCATCTTCTATTGTTTGG'
        'TTGAGTTAAGAACTGATCATGAATAAATACATTTCATTG', '[TATC]12', '12', '12', None, None, 'ATCT ATCA'
        ' TCTT CTAT TGTT T GGTT GAGT TAAG AACT GATC ATGA ATAA ATAC ATTT CATT G', 'powerseq'
    )
])
def test_strobj_DYS533(sequence, bracketed, conc, lus, sec, tert, flank_3p, kit):
    marker = STRMarkerObject('DYS533', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_3p == flank_3p


def test_strobj_DYS522():
    marker = STRMarkerObject(
        'DYS522', 'AGTTAGAGAGAGAGATGATGGATAGATAAATAGATAGATGATAGATGAATAGATAGGCGGGTAATAGATTTTATATAG'
        'ATAGATGATAGCTAGATAATGGATAGACATAGGTGACAGATGATAAATACATAGATAAATAGATGATAGATAGATAGATAGATAGATA'
        'GATAGATAGATAGATAGATAGATAGACAGATGTCCACCATGAGGTTC', uas=False, kit='forenseq'
    )
    assert marker.annotation == 'ATA GATG [ATAG]12'
    assert marker.canonical == 12
    assert marker.designation == ('12', None, None)


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, kit', [
    (
        'AAGGTGATAGATATACAGATAGATAGATACATAGGTGGAGACAGATAGATGATAAATAGAAGATAGATAGATAGATAGATAGATAGAT'
        'AGATAGATAGATAGATAGAAAGTATAAGTAAAGAGATGAT', '[GATA]2 TACA [GATA]3 CATA GGTG GAGA CAGA '
        'TAGA TGAT AAAT AGAA [GATA]11', '11', '11', None, None, 'forenseq'
    ),
    (
        'ATAAATAGAAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAGTATAAGTAAAGAGATGA',
        'ATAA ATAG AA [GATA]12', '12', '12', None, None, 'powerseq'
    )
])
def test_strobj_DYS439(sequence, bracketed, conc, lus, sec, tert, kit):
    marker = STRMarkerObject('DYS439', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
