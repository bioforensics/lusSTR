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


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert', [
    (
        'TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGATACATAGATAATACAGATGAGAGTTGGATACA'
        'GAAGTAGGTATAATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACACA'
        'CACATAGATAATACAGAT', '[TAGA]11 [CAGA]3 TACA TAGA TAAT ACAG ATGA GAGT TGGA TACA GAAG TAGG'
        ' TATA ATGA [TAGA]13 [CAGA]4', 31, 13, 11, None
    )
])
def test_strobj_DYS389II(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject('DYS389II', sequence, uas=False, kit='forenseq')
    assert marker.annotation == bracketed
    assert marker.canonical == conc
    assert marker.designation == (lus, sec, tert)


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
        'AAATAGAAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAGTATAAGTAAAGAGATGATGG',
        'AAAT AGAA [GATA]12', '12', '12', None, None, 'powerseq'
    )
])
def test_strobj_DYS439(sequence, bracketed, conc, lus, sec, tert, kit):
    marker = STRMarkerObject('DYS439', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_3p, kit', [
    (
        'ATGCCCATCCGGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGTCTGTCTGTCTATCTATCTATCTATCATCTATCATC'
        'TGTGAATGATGTCTATCTACTTATCTATGAATGATATTTATCTGTGGTTATCTATCTATCTATATCATCTGTGAATGACAGGGTCTTC'
        'CTCTG', '[TCTA]9 [TCTG]3 [TCTA]4', '16', '9', None, None, 'TCA TCTA TCAT CTGT GAAT GATG '
        '[TCTA]2 CTTA TCTA TGAA TGAT ATTT ATCT GTGG TTAT [CTAT]3 ATCA TCTG TGAA TGAC AGGG TCTT '
        'CCTC TG', 'forenseq'
    ),
    (
        'CCATCCGGTCTATCTATCTATCTATCTATCTATCTATCTATCTGTCTGTCTATCTATCTATCTATCATCTATCATCTGTGAATGATGT'
        'CTATCTACTTATCTATGAATGATATTTATCTGTGGTTATCTATCTATCTATA', '[TCTA]8 [TCTG]2 [TCTA]4', '14',
        '8', None, None, 'TCA TCTA TCAT CTGT GAAT GATG [TCTA]2 CTTA TCTA TGAA TGAT ATTT ATCT '
        'GTGG TTAT [CTAT]3 A', 'powerseq'
    )
])
def test_strobj_DYS437(sequence, bracketed, conc, lus, sec, tert, flank_3p, kit):
    marker = STRMarkerObject('DYS437', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_3p == flank_3p


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_3p, kit', [
    (
        'TTAAACCTACCAATCCCATTCCTTAGTAAATAATAATAATAATAATAATAATAATAATAATAATAATAAATAAATGGTGATACAAGAA'
        'AAAAATTTGTTTTCCTTCTTGGCTTTTAAATAACAAACACTTGAAATCAAATTAGTTGTTTTTAAAAGCTAGATTAATGAAGAA',
        '[ATA]13', '13', '13', None, None, 'AAT AAA TGG TGA TAC AAG [AAA]2 ATT TGT TTT CCT TCT '
        'TGG CTT TTA AAT AAC AAA CAC TTG AAA TCA AAT TAG TTG TTT TTA AAA GCT AGA TTA ATG AAG AA',
        'forenseq'
    ),
    (
        'TAAATAATAATAATAATAATAATAATAATAATAATAATAATAAATAAATGGTGATACAAGAAAAAAATTTGTTTTCCTTCTTGGCTTT'
        'TAAATAACAAACACTTGAAATCAAATTAGTT', '[ATA]13', '13', '13', None, None, 'AAT AAA TGG TGA '
        'TAC AAG [AAA]2 ATT TGT TTT CCT TCT TGG CTT TTA AAT AAC AAA CAC TTG AAA TCA AAT TAG TT',
        'powerseq'
    )
])
def test_strobj_DYS392(sequence, bracketed, conc, lus, sec, tert, flank_3p, kit):
    marker = STRMarkerObject('DYS392', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_3p == flank_3p


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_3p, kit', [
    (
        'ATATCTGTCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGCCTATCTGCCTGCCTAC'
        'CTATCCCTCTAT', '[TCTG]3 [TCTA]13 TCTG', '13', '13', None, None, 'CCTA TCT [GCCT]2 ACCT '
        'ATCC CTCT AT', 'forenseq'
    ),
    (
        'ATATCTGTCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGCCTATCTGCCTGCCTACCTATCCCTCTAT'
        'GGCAATTGCTTGCAACCAGGGAGATTTTAT', '[TCTG]3 [TCTA]10 TCTG', '10', '10', None, None,
        'CCTA TCT [GCCT]2 ACCT ATCC CTCT AT GGCA ATTG CTTG CAAC CAGG GAGA TTTT AT', 'powerseq'
    )
])
def test_strobj_DYS391(sequence, bracketed, conc, lus, sec, tert, flank_3p, kit):
    marker = STRMarkerObject('DYS391', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_3p == flank_3p


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_3p, kit', [
    (
        'TGGTCAATCTCTGCACCTGGAAATAGTGGCTGGGGCACCAGGAGTAATACTTCGGGCCATGGCCATGTAGTGAGGACAAGGAGTCCAT'
        'CTGGGTTAAGGAGAGTGTCACTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTATCTATCTATCTA',
        '[TCTA]11 CCTA [TCTA]3', '14', '11', None, None, '', 'forenseq'
    ),
    (
        'CTGGGTTAAGGAGAGTGTCACTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTATCTATCTATCTAAAA'
        'CACTATATATATATAACACTATATATATAATACTATATATATATTAAAAAACACTAT', '[TCTA]11 CCTA [TCTA]3',
        '14', '11', None, None, 'AA ACAC [TATA]3 ACAC [TATA]2 TA ATAC [TATA]2 TATT AAAA AACA '
        'CTAT', 'powerseq'
    )
])
def test_strobj_DYS19(sequence, bracketed, conc, lus, sec, tert, flank_3p, kit):
    marker = STRMarkerObject('DYS19', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_3p == flank_3p


def test_strobj_DYS458():
    marker = STRMarkerObject(
        'DYS458', 'GAAAGAAAGAAAAGGAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA'
        'GAAAGAAAGGAGGGTGGGCGTGGTGGCTCATGCTTGTAATGCCAGAACTTTGGGAGGCCGAGGT', uas=False,
        kit='powerseq'
    )
    assert marker.annotation == (
        '[GAAA]3 AG GAAG [GAAA]17 GGAG GGTG GGCG TGGT GGCT CATG CTTG TAAT GCCA GAAC TTTG GGAG '
        'GCCG AGGT'
    )
    assert marker.canonical == 17
    assert marker.designation == ('17', None, None)


def test_strobj_HPRTB():
    marker = STRMarkerObject(
        'HPRTB',
        'CTAGAACTTATCTTCTGTAAATCTGTCTCTATTTCCATCTCTGTCTCCATCTTTGTCTCTATCTCTATCTGTCTATCTCTATCTATCT'
        'ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAAAGCAAATTCATGCCCTTCTCCTATTT', uas=False,
        kit='forenseq'
    )
    assert marker.annotation == '[ATCT]12'
    assert marker.canonical == 12
    assert marker.designation == ('12', None, None)
    assert marker.flank_5p == (
        'CTAG AACT TATC TTCT GTAA ATCT GTCT CTAT TTCC ATCT CTGT CTCC ATCT TTGT CTCT ATCT CTAT '
        'CTGT CTAT C TCT'
    )


def test_strobj_DXS8378():
    marker = STRMarkerObject(
        'DXS8378',
        'AGTGAGCTGAGATGGTGCCACTGAACTCCAGCCTGGGCGACAAGAGCGAAACTCCAACTCAAAAAATAAATAAATAAAATATAGATAG'
        'ATAGATAGATAGATAGATAGATAGATAGATAGATAGTGACCTGCCAGGAGCAGGGGACCACCGGGTTGCCTAAGGAGGGGTGAACTGT'
        'CCCAGGATGGAAATGAAACA', uas=False, kit='forenseq'
    )
    assert marker.annotation == '[ATAG]11'
    assert marker.canonical == 11
    assert marker.designation == ('11', None, None)
    assert marker.flank_5p == (
        'AGTG AGCT GAGA TGGT GCCA C TGAA CTCC AGCC TGGG CGAC AAGA GC G AAA CTCC AACT CAAA [AAAT]3'
        ' AAAA T'
    )
    assert marker.flank_3p == (
        'TGAC CTGC CAGG AGCA GGGG ACCA CC G GGTT GCCT AAGG AGGG GTGA ACTG TCCC AGGA TGGA AATG '
        'AAAC A'
    )


def test_strobj_DXS7132():
    marker = STRMarkerObject(
        'DXS7132',
        'TCCAGAGAAACAGAACCAATAGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGTCAGATAGA'
        'TGAGAGGGGATTTATTATGAAAATGGGCTCACACTATTAAGGAGGCTAAGAAGTTCCACAGTAT', uas=False,
        kit='forenseq'
    )
    assert marker.annotation == '[TAGA]13'
    assert marker.canonical == 13
    assert marker.designation == ('13', None, None)
    assert marker.flank_3p == (
        'CAGT C [AGAT]2 GAGA GGGG ATTT ATTA TGAA AAT G GGCT CACA CTAT TAAG GAGG CTAA GAAG TTCC '
        'ACAG TAT'
    )


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert', [
    (
        'AAGAAAGAAAGAGAAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAA'
        'AAGAGAATAGAAAAGAAGAGAAGAGAAAAGAGAAAAGAAAAAAGAAAAGAAA', '[AAGA]3 gaaagga [AAGA]17 AAAG '
        'AGAA TAGA AAAG AAGA GAAG AGAA AAGA GAAA AGAA AAAA GAAA AGAA A', '21', '17', '0', None
    ),
    (
        'AAGAAAGAAAGAGAAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAA'
        'AGAAAGGAAGAAAGAAAGAAAGGAAGAAAAGAGAATAGAAAAGAAGAGAAGAGAAAAGAGAAAAGAAAAAAGAAAAGAAA',
        '[AAGA]3 gaaagga [AAGA]18 AAGG [AAGA]3 AAGG AAGA AAAG AGAA TAGA AAAG AAGA GAAG AGAA AAGA '
        'GAAA AGAA AAAA GAAA AGAA A', '28', '18', '1', None
    )
])
def test_strobj_DXS10135(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject('DXS10135', sequence, uas=False, kit='forenseq')
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_5p, flank_3p', [
    (
        'TGTGTGTGCATGCATACACACACAGAGAGAGAGAGAGAGAAAAAGAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAG'
        'AAAGAAAGAAAGAAAGAAAGGAAGAAAGAAAGGAAGAAAATAGAACAAATCAGCTTATATTCAGTATTTTTTTAGTATTTTCTGTGTC'
        'AGCTC', '[AAGA]15 AAGG [AAGA]2', '18', '15', '0', None, '[TGTG]2 CATG CATA CACA CA C '
        '[AGAG]4 AAAA AG', 'AAGG A AAGG AAGA AAAT AGAA CAAA TCAG CTTA TATT CAGT ATTT TTTT AGTA '
        'TTTT CTGT GTCA G TC'
    ),
    (
        'TGTGTGTGCATGCATACACACACAGAGAGAGAGAGAGAAAAAGAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAA'
        'AGAAAGAAAGAAAGAAAGGAAGAAAGAAAGGAAGAAAATAGAACAAATCAGCTTATATTCAGTATTTTTTTAGTATTTTCTGTGTCAG'
        'CTC', 'GA [AAGA]14 AAGG [AAGA]2', '17.2', '14', '0', None, '[TGTG]2 CATG CATA CACA CA C '
        '[AGAG]3 AGAA AAAG AA', 'AAGG A AAGG AAGA AAAT AGAA CAAA TCAG CTTA TATT CAGT ATTT TTTT '
        'AGTA TTTT CTGT GTCA G TC'
    )
])
def test_strobj_DXS10074(sequence, bracketed, conc, lus, sec, tert, flank_5p, flank_3p):
    marker = STRMarkerObject('DXS10074', sequence, uas=False, kit='forenseq')
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5p
    assert marker.flank_3p == flank_3p


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, kit', [
    (
        'CTATCTATCTATCTATTCATCCATCTAATCTATCCATTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTACCT'
        'ACCTATCTATCTATAGATCTATCTATCTATCT', 'C [TATC]3 TATT [CATC]2 TAAT CTAT CCAT [TCTA]11 '
        '[CCTA]3 [TCTA]2 TAGA [TCTA]3 TCT', '11', '11', None, None, 'forenseq'
    ),
    (
        'CATTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTACCTACCTATCTATCTATAGATCTATCTATCTATCTTA'
        'AATTTGGAAATTCTCCTCAGCATAACATTTTAATGATGATTCCTAGGATACAAGTGATGTGCTGAAAGTATCAATGTGTATCAGAAAA'
        'CCAACATCTCTGCTTAGGTCTC', 'CAT [TCTA]11 CCTA CCTA CCTA [TCTA]2 TAGA [TCTA]3 TCT', '11',
        '11', None, None, 'powerseq'
    ),
    (
        'CATTCTATCTATCTATCTATCCATCTATCTATCTATCTATCTATCTACCTACCTACCTATCTATCTATAGATCTATCTATCTATCTTA'
        'AATTTGGAAATTCTCCTCAGCATAACATTTTAATGATGATTCCTAGGATACAAGTGATGTGCTGAAAGTATCAATGTGTATCAGAAAA'
        'CCAACATCTCTGCTTAGGTCTC', 'CAT [TCTA]4 TCCA [TCTA]6 CCTA CCTA CCTA [TCTA]2 TAGA [TCTA]3 '
        'TCT', '11', '6', None, None, 'powerseq'
    )
])
def test_strobj_Y_GATA_H4(sequence, bracketed, conc, lus, sec, tert, kit):
    marker = STRMarkerObject('Y-GATA-H4', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, flank_5p, kit', [
    (
        'AGAAATGGATGACAGTAAAATGAAAACATTGCAATGTGTATACTCAGAAACAAGGAAAGATAGATAGATGATAGATAGATAGATAGAC'
        'AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACAGACAGACAGACAGATAGAT',
        '[TAGA]4 CAGA [TAGA]11 [CAGA]8 TAGA T', '24', '11', '8', '0', 'AG AAAT GGAT GACA GTAA '
        'AATG AAAA CATT GCAA TGTG TATA CTCA GAAA CAAG GAAA [GATA]2 GATG A', 'forenseq'
    ),
    (
        'AACAAGGAAAGATAGATAGATGATAGATAGATAGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAC'
        'AGACAGACAGACAGACAGACAGACAGACAGATAGATAGAATATATTATGGGGTACCAAAATGCAGGGCCCAAAAATGTGTAAAATATA'
        'TGTGT', '[TAGA]4 CAGA [TAGA]11 [CAGA]8 [TAGA]2', '24', '11', '8', '0',
        'AA CAAG GAAA [GATA]2 GATG A', 'powerseq'
    )
])
def test_strobj_DYS390(sequence, bracketed, conc, lus, sec, tert, flank_5p, kit):
    marker = STRMarkerObject('DYS390', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5p


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert, kit', [
    (
        'TTTCTTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTATTTCTTTCCCTTCCTTCCTT'
        'CCTTCCTTCCTTTCTTTCTCTTTCCTCTTTCTCTTTCTTCTCTTTCTTTCTTTTTCTCTTTTTCTCTTTCTTTCTTTTTTACTTTCTT'
        'TCTCCTTCCTTCCTTCCTTTCTGAATTTCATTTCTTTTCTTT', 'TTTC TTTT TCTC [TTTC]13 TTTA [TTTC]2 '
        '[CCTT]6 TCTT TCTC TTTC CTCT TTCT CTTT CTTC [TCTT]3 TTTC TCTT TTTC [TCTT]3 TTTT ACTT TCTT'
        ' TCTC [CTTC]3 CTT', '16', '13', None, None, 'forenseq'
    ),
    (
        'TTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCCCTTCCTTCCTTCCTTCCTTCCTTTCTT'
        'TCTCTTTCCTCTTTCTCTTTCTTCTCTTTCTTTCTTTTTCTCTTTTTCTCTTTCTTTCTTTTTTACTTTCTTTCTCCTTCCTTCCTTC'
        'CTTTCTGAATTTCATTTCTTTTCTTT', 'TTTT TCTC [TTTC]13 [CCTT]6 TCTT TCTC TTTC CTCT TTCT CTTT '
        'CTTC [TCTT]3 TTTC TCTT TTTC [TCTT]3 TTTT ACTT TCTT TCTC [CTTC]3 CTT', '13', '13',
        None, None, 'powerseq'
    )
])
def test_strobj_DYS385(sequence, bracketed, conc, lus, sec, tert, kit):
    marker = STRMarkerObject('DYS385A-B', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert', [
    (
        'GAGATAGAGACATGGATAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATATA'
        'GAGATAGAGAGATAGAGATAGAGATAGATAGATAGAGAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATA'
        'GAGA', '[AGAGAT]11 [ATAGAG]2 [AGATAG]3 ATAGAT AGAGAA [AGAGAT]8', 19, 11, 8, None
    ),
    (
        'GAGATAGAGACATGGATAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGA'
        'GATATAGAGATAGAGAGATAGAGATAGAGATAGATAGATAGAGAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGAGAGAGATA'
        'GAGATAGAGA', '[AGAGAT]12 [ATAGAG]2 [AGATAG]3 ATAGAT AGAGAA [AGAGAT]5 AGAGAG [AGAGAT]2',
        20, 12, 2, None
    )
])
def test_strobj_DYS448(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject('DYS448', sequence, uas=False, kit='forenseq')
    assert marker.annotation == bracketed
    assert marker.canonical == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert', [
    (
        'TGCATGCACATACACATAACTAGATAGACTGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGA'
        'CAGATAGATA', '[TAGA]2 CTGA CAGA [TAGA]10 [CAGA]4 TAGA', 18, 10, 2, None
    ),
    (
        'TGCATGCACATACACATATCTAGATCGACTGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGA'
        'CAGATAGATA', 'TAGA TCGA CTGA CAGA [TAGA]10 [CAGA]4 TAGA', 18, 10, 1, None
    )
])
def test_strobj_DXS10103(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject('DXS10103', sequence, uas=False, kit='forenseq')
    assert marker.annotation == bracketed
    assert marker.canonical == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize('sequence, bracketed, conc, lus, sec, tert', [
    (
        'TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGATACATAGATAATACAGATGAGAGTTGGATACA'
        'GAAGTAGGTATAATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACACA'
        'CACATAGATAATACAGAT', '[TAGA]11 [CAGA]3 TACA TAGA TAAT ACAG ATGA GAGT TGGA TACA GAAG TAGG'
        ' TATA ATGA [TAGA]13 [CAGA]4', 31, 13, 11, None
    )
])
def test_strobj_DYS389II(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject('DYS389II', sequence, uas=False, kit='forenseq')
    assert marker.annotation == bracketed
    assert marker.canonical == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize('sequence, bracketed, lus, sec, tert, flank_5, flank_3, kit', [
    (
        'GTCTCAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGAGAGAGGAAAGAAAGAGAAAAAG'
        'AAAAGAAATAGTAGCAACTGTTATTGTAAGA', '[AGAA]14 AAAG AGAG AG', '14', '1', None, 'G TCTC',
        'GA [AAGA]2 GAAA AAGA AAAG AAAT AGTA GCA A CTGT TATT GT AAGA', 'forenseq'
    ),
    (
        'AGGCTGCAGTGAGCCATGTTCATGCCACTGCACTTCACTCTGAGTGACAAATTGAGACCTTGTCTCAGAAAGAAAGAAAGAAAGAAAG'
        'AAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGAGAGAGGAAAGAAAGAGAAAAAGAAAAGAAATAG'
        'TAGCAACTGTTATTG', '[AGAA]18 AAAG AGAG AG', '18', '1', None, 'A GGCT GCAG TGAG CCAT GTTC '
        'ATGC CACT GCAC TTCA CTCT GAGT GACA AATT GAGA CCTT G TCTC', 'GA [AAGA]2 GAAA AAGA AAAG '
        'AAAT AGTA GCA A ACTG TTAT TG', 'powerseq'
    )
])
def test_strobj_D18S51(sequence, bracketed, lus, sec, tert, flank_5, flank_3, kit):
    marker = STRMarkerObject('D18S51', sequence, uas=False, kit=kit)
    assert marker.annotation == bracketed
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5
    assert marker.flank_3p == flank_3


@pytest.mark.parametrize('locus, sequence, bracketed, conc, lus, sec, tert, flank_5, flank_3', [
    (
        'CSF1PO', 'CTAAGTACTTCCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAATCTATCTATCTTCTA'
        'TCTATGAAGGCAGTTACTGTTAATATCTTCATTTTACAGGTAGGAAAACTGAGACACAGGGTGGTTAGCAACCTGCTAGTCCTTGGCA'
        'GACTCAG', '[ATCT]12', '12', '12', '0', None, 'CTA AGTA CT TCCT', 'A [ATCT]3 T [CTAT]2 '
        'GAAG GCAG TTAC TGTT AATA TCTT CATT TTAC AGGT AGGA AAAC TGAG ACAC AGGG TGGT TAG CA ACCT '
        'GCTA GTCC TTGG CAGA CTCA G'
    ),
    (
        'D10S1248', 'CCCCAGGACCAATCTGGTCACAAACATATTAATGAATTGAACAAATGAGTGAGTGGAAGGAAGGAAGGAAGGAAGG'
        'AAGGAAGGAAGGAAGGAAGGAAGGAAGGAA', '[GGAA]13', '13', '13', None, None, 'CCCC AGGA CCAA '
        'TCTG GTCA CAAA CATA TTAA TGAA TT GAAC AAAT [GAGT]2', ''
    ),
    (
        'D13S317', 'TTCTTTAGTGGGCATCCGTGACTCTCTGGACTCTGACCCATCTAACGCCTATCTGTATTTACAAATACATTATCTA'
        'TCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTCTTTTTGGGCTGCCTATGGCTCAACCCA'
        'AGTTGAAGGAGGAGATTT', '[TATC]9 [AATC]2 [ATCT]3 TTCT GTCT GTC', '9', '9', '3', '1', 'TT '
        'CTT TAGT GGGC ATCC G TGAC TCTCT GGAC TC TGAC CCAT CTAA C G CCT ATCT GTAT TTAC AAAT ACAT',
        'TTTT TGGG CTGC CTAT GGCT CAAC CCAA GTTG AAGG AGGA GATT T'
    ),
    (
        'D16S539', 'GTGCACAAATCTAAATGCAGAAAAGCACTGAAAGAAGAATCCAGAAAACCACAGTTCCCATTTTTATATGGGAGCAA'
        'ACAAAGGCAGATCCCAAGCTCTTCCTCTTCCCTAGATCAATACAGACAGACAGACAGGTGGATAGATAGATAGATAGATAGATAGATA'
        'GATAGATATCAT', '[GATA]9', '9', '9', '0', None, 'GT GCAC AAAT CTAA ATGC AGAA AAGC ACTG '
        'AAAG AAGA ATCC AG AAAA CCAC AGTT CCCA TTTT TATA TGGG AG [CAAA]2 GGCA GATC CCAA G CTCT TC'
        ' CTCT T CCCT AGAT CAAT [ACAG]4 GTG', 'TCAT'
    ),
    (
        'D19S433', 'AAGTTCTTTAGCAGTGATTTCTGATATTTTGGTGCACCCATTACCCGAATAAAAATCTTCTCTCTTTCTTCCTCTCT'
        'CCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTACCTTCTTTCCTTCAACAGAATCTTATTCTGTTGCCCAGGC'
        'TGGAGTGCA', 'CT CTCT TTCT TCCT CTCT [CCTT]11 CCTA CCTT CTTT CCTT', '13', '11', '1', '0',
        'AAG TTCT TTAG CAGT GATT TCTG ATAT TTTG GTGC ACCC ATTA CCCG AATA AAAA TCTT', 'CAAC AGAA '
        'TCTT ATTC TGTT GCCC AGGC TGGA GTGC A'
    ),
    (
        'D1S1656', 'GAAATAGAATCACTAGGGAACCAAATATATATACATACAATTAAACACACACACATCTATCTATCTATCTATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCTACATCACACAGTTGACCCTTGA', 'CA [CACA]2 [TCTA]13', '13', '13',
        '0', '0', 'G AAAT AGAA TCAC TAGG GAAC CAAA [TATA]2 CATA CAAT TAAA', 'CATC AC ACA GTTG '
        'ACCC TTGA'
    ),
    (
        'D21S11', 'TGAATTGCCTTCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATC'
        'TATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCTATCGTCTATCTATCCAGTCT'
        'ATCTACCTCCTATTAGTCT', '[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11 TA'
        ' TCTA', '31.2', '11', '5', '6', 'TG AATT GCCT', 'TC G [TCTA]2 T CCAG [TCTA]2 CCTC CTAT T'
        ' AGTC T'
    ),
    (
        'D22S1045', 'CCTTCTTATAGCTGCTATGGGGGCTAGATTTTCCCCGATGATAGTAGTCTCATTATTATTATTATTATTATTATTA'
        'TTATTATTATTATTACTATTATTGTTATAAAAATATTG', '[ATT]13 ACT [ATT]2', '16', '13', None, None,
        'CCT TCT TAT AGC TGC TAT GGG GGC TAG ATT TTC CCC [GAT]2 [AGT]2 CTC', 'GTT ATA AAA ATA TTG'
    ),
    (
        'D2S1338', 'CTAGCATGGTACCTGCAGGTGGCCCATAATCATGAGTTATTCAGTAAGTTAAAGGATTGCAGGAGGGAAGGAAGGAC'
        'GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGCAGGCAGGCAGGCAGGCAGGCAGGCAAGGCCAAGCCAT'
        'TTCTGTTTCCAAATCCACTGGCTCCCTCCCACAGCT', '[GGAA]2 GGAC [GGAA]12 [GGCA]7', '22', '12', '1',
        '7', 'CT AGCA TGGT ACCT GCAG GTGG CCCA TAAT C ATGA GTTA TTCA GTAA GTTA AAGG ATTG CAG GAG',
        'AGGC CAAG CCAT TT CTGT TTCC AAAT CCAC TGGC [TCCC]2 ACAG CT'
    ),
    (
        'D2S441', 'TGCACCCAACATTCTAACAAAAGGCTGTAACAAGGGCTACAGGAATCATGAGCCAGGAACTGTGGCTCATCTATGAAA'
        'ACTTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCA', '[TCTA]10', '10', '10', '0', None,
        'TG CACC CAAC ATTC TAAC AAAA GGCT GTAA CAAG GGCT ACAG GAA T CATG AG CCAG G AACT GTGG CTCA'
        ' TCTA TGAA AACT', 'TATC A'
    ),
    (
        'D3S1358', 'TGCCCACTTCTGCCCAGGGATCTATTTTTCTGTGGTGTGTATTCCCTGTGCCTTTGGGGGCATCTCTTATACTCATG'
        'AAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAT'
        'CTATGAGACAGGGTCTTGCTCTGTC', 'TCTA [TCTG]3 [TCTA]13', '17', '13', '3', None, 'TGCC CACT '
        'TCTG CCCA GGGA TCTA TTTT TCTG TGGT GTGT ATTC CCTG TGCC TTTG GGGG CATC TCTT ATAC TCAT '
        'GAAA TCAA CAGA GGCT TGCA TGTA', 'TGAG ACAG GGTC TTGC TC TGTC'
    ),
    (
        'D5S818', 'AACATTTGTATCTTTATCTGTATCCTTATTTATACCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC'
        'TATCTATCTTCAAAATATTACGTAAGGATACCAAAGAGGAAAATCACCCTTGTCACATACTTGCTATTAAAATATACTTTTATTAGTA'
        'CA', 'CTCT [ATCT]12', '12', '12', None, None, 'AA CATT TGTA TCTT TATC TGTA TCCT T ATTT '
        'ATAC', 'TCAA AAT ATTA CG TAAG GATA CCAA AGAG GAAA ATCA CCCT TGTC ACAT ACTT GCTA T TAAA '
        'ATAT ACTT TTAT TAGT ACA'
    ),
    (
        'D7S820', 'AGAATTGCACCAAATATTGGTAATTAAATGTTTACTATAGACTATTTAGTGAGATAAAAAAAAACTATCAATCTGTCT'
        'ATCTATCTATCTATCTATCTATCTATCTATCGTTAGTTCGTTCTAAACTATGACAAGTGTTCTATCATACCCTTTATATATATTAACC'
        'TTAAAATAACTC', 'AAAC TATC AATC TGTC [TATC]8', '8', '8', '1', 0, 'AGAA TTGC ACCA A '
        'ATAT TGGT AATT AAAT GTTT ACTA T AGAC T ATTT AGTG AGAT AAAAAA', 'GTTA [GTTC]2 TAAA CTAT '
        'GACA AGTG TTCT ATCA TACC CTTT [ATAT]2 TAAC CTTA AAAT AACT C'
    ),
    (
        'FGA', 'GTCTGAAATCGAAAATATGGTTATTGAAGTAGCTGCTGAGTGATTTGTCTGTAATTGCCAGCAAAAAAGAAAGGAAGAAAG'
        'GAAGGAAGGAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAGAAAAAAGAAAGAAAGAAA',
        '[GGAA]2 GGAG [AAAG]12 AGAA AAAA [GAAA]3', '20', '12', '3', '0', 'G TCTG AAAT CGAA AATA '
        'TGGT TATT GAAG TAGC TGCT GAGT GATT TGTC TGTA ATTG CCA GCAA AAAA GAAA GGAA GAAA', ''
    )
])
def test_new_power_config(
    locus, sequence, bracketed, conc, lus, sec, tert, flank_5, flank_3
):
    marker = STRMarkerObject(locus, sequence, uas=False, kit='powerseq')
    assert marker.annotation == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5
    assert marker.flank_3p == flank_3
