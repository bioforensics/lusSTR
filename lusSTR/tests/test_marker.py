# -------------------------------------------------------------------------------------------------
# Copyright (c) 2020, DHS.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import pytest
import lusSTR
from lusSTR.scripts.marker import STRMarkerObject


@pytest.mark.parametrize(
    "sequence, bracket_form",
    [
        (
            "TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATC"
            "TATCTATCTATCTATCTATCTATCTATCTATCTATATCTA",
            "[TCTA]5 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]10 TA TCTA",
        ),
        (
            "TCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATC"
            "TATCTATCTATCTATCTATCTATCTATCTAACTATCTA",
            "[TCTA]4 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]9 ACTA TCTA",
        ),
        (
            "TCTATCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATC"
            "TATCTATCTATCTATCTATCTATCTATCTA",
            "[TCTA]6 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]8",
        ),
        (
            "TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCT"
            "ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATAACTA",
            "[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]11 TA ACTA",
        ),
        (
            "TCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCT"
            "ATCTATCTATCTATCTATCTATCTATCTAACTATCTA",
            "[TCTA]4 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]9 ACTA TCTA",
        ),
    ],
)
def test_D21_bracket(sequence, bracket_form):
    marker = STRMarkerObject("D21S11", sequence, "uas")
    assert marker.convert == bracket_form


def test_D19_convert():
    uas_sequence = (
        "AAGGAAAAGGTAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAGAGAAGAAGAAAGAGAG"
    )
    marker = STRMarkerObject("D19S433", uas_sequence, "uas")
    assert marker.convert == "CT CTCT TTCT TCTT CTCT [CCTT]14 CCTA CCTT TT CCTT"


def test_D1_convert():
    uas_sequence = "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGTGTATGTG"
    marker = STRMarkerObject("D1S1656", uas_sequence, "uas")
    assert marker.convert == "CA CATA CACA [TCTA]11"


@pytest.mark.parametrize(
    "sequence, bracket_form",
    [
        (
            "AAAAGAAAGAAAAGAAAAGAAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGA",
            "AAAAG [AAAGA]3 A [AAAGA]11",
        ),
        ("GAAAAGAAAAGAAAAGA", "GA [AAAGA]3"),
        (
            "AAAAGAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAAAAAAGAAAAGA",
            "AAAAG [AAAGA]7 AAAAA [AAAGA]2",
        ),
        ("AAAAGAAAAAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGA", "AAAAG AAAAA [AAAGA]8"),
    ],
)
def test_PentaD_convert(sequence, bracket_form):
    marker = STRMarkerObject("PENTA D", sequence, "uas")
    assert marker.convert == bracket_form


def test_FGA_convert():
    uas_sequence = (
        "TTTCTTTCTTTCTTTCTTTTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTGTCTGTCTGTCTTTCTTTCTTTCTTTCTT"
        "TCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTCCTTCCTTCCTTTCTTTCTTTCTCCTTCCTTCCTTCCTTCC"
    )
    convert = "[GGAA]4 GGAG [AAAG]3 [GAAG]3 [AAAG]15 [ACAG]3 [AAAG]9 AA AAAA [GAAA]4"
    marker = STRMarkerObject("FGA", uas_sequence, "uas")
    assert marker.convert == convert


def test_THO1():
    marker = STRMarkerObject("TH01", "AATGAATGAATGAATGAATGATGATGAATGAATGAATG", "uas")
    assert marker.convert == "[AATG]5 ATG ATG [AATG]3"


@pytest.mark.parametrize(
    "sequence, lus_allele, sec_allele, tert_allele",
    [
        (
            "TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATC"
            "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCAA",
            "13",
            "5",
            "5",
        ),
        (
            "TCTATCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATC"
            "TATCTATCTATCTATCTATCTATCTATCTA",
            "8",
            "6",
            "5",
        ),
        (
            "TCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATC"
            "TATCTATCTATCTATCTATCTATCTATCTATCTATCTA",
            "12",
            "4",
            "5",
        ),
        (
            "TCTATCTATCTATCTATCTATCTGTCTGTCTTTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATAT"
            "CTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA",
            "11",
            "5",
            "3",
        ),
    ],
)
def test_D21_convert(sequence, lus_allele, sec_allele, tert_allele):
    marker = STRMarkerObject("D21S11", sequence, "uas")
    lus, sec, tert = marker.designation
    assert str(lus) == lus_allele
    assert str(sec) == sec_allele
    assert str(tert) == tert_allele


def test_D21_lus_sec():
    sequence = (
        "TCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATC"
        "TATCTATCTATCTATCTATCTATCTATCTATCTA"
    )
    marker = STRMarkerObject("D21S11", sequence, "uas")
    lus, sec, tert = marker.designation
    assert str(lus) == "10"
    assert str(sec) == "4"
    assert str(tert) == "6"


@pytest.mark.parametrize(
    "locus, sequence, forward_bracket, lus, sec, tert",
    [
        (
            "D21S11",
            "TCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATC"
            "TATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCAA",
            "[TCTA]5 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]13 TA TCAA",
            "13",
            "5",
            "5",
        ),
        (
            "D7S820",
            "GATAGATAGATAGATAGACAGATAGATAGATAGATAGATAGATA",
            "[TATC]6 TGTC [TATC]4",
            "6",
            "1",
            "0",
        ),
        (
            "D9S1122",
            "TAGATCGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGA",
            "TAGA TCGA [TAGA]10",
            "10",
            "None",
            "None",
        ),
        (
            "D16S539",
            "GATAGATAGATAGATAGATAGACAGATAGATAGATAGATA",
            "[GATA]5 GACA [GATA]4",
            "5",
            "1",
            "None",
        ),
        (
            "D13S317",
            "TATCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTC",
            "[TATC]8 [AATC]2 [ATCT]3 TTCT GTCT GTC",
            "8",
            "3",
            "1",
        ),
        (
            "D13S317",
            "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGT"
            "CTTTTTGGGCTGCCTATATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTC",
            "[TATC]10 [AATC]2 [ATCT]3 TTCT [GTCT]2 TTTT GGGC TGCC TA [TATC]7 [AATC]2 [ATCT]3 TTCT"
            " GTCT GTC",
            "10",
            "3",
            "2",
        ),
    ],
)
def test_convert_and_lus(locus, sequence, forward_bracket, lus, sec, tert):
    marker = STRMarkerObject(locus, sequence, "uas")
    assert marker.convert == forward_bracket
    lus_out, sec_out, tert_out = marker.designation
    assert str(lus_out) == lus
    assert str(sec_out) == sec
    assert str(tert_out) == tert


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert",
    [
        (
            "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGATACATAGATAATACAGATGAGAGTTGGATACA"
            "GAAGTAGGTATAATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACACA"
            "CACATAGATAATACAGAT",
            "[TAGA]11 [CAGA]3 TACA TAGA TAAT ACAG ATGA GAGT TGGA TACA GAAG TAGG"
            " TATA ATGA [TAGA]13 [CAGA]4",
            31,
            13,
            11,
            None,
        )
    ],
)
def test_strobj_DYS389II(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject("DYS389II", sequence, "straitrazor", kit="forenseq")
    assert marker.convert == bracketed
    assert marker.canonical == conc
    assert marker.designation == (lus, sec, tert)


def test_strobj_CSF1PO():
    marker = STRMarkerObject(
        "CSF1PO",
        "CTTCCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAATCTATCTATCTT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.uas_sequence == "AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT"
    assert marker.forward_sequence == "ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT"
    assert marker.convert == "[ATCT]12"
    assert marker.convert_uas == "[AGAT]12"
    assert marker.canonical == 12, " "
    assert marker.designation == ("12", "0", None)
    assert marker.flank_5p == "CT TCCT"


def test_strobj_D10S1248():
    marker = STRMarkerObject(
        "D10S1248",
        "TTGAACAAATGAGTGAGTGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAATGAAGA"
        "CAATACAACCAGAGTT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.uas_sequence == "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA"
    assert marker.forward_sequence == "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA"
    assert marker.convert == "[GGAA]13"
    assert marker.convert_uas == "[GGAA]13"
    assert marker.canonical == 13, " "
    assert marker.designation == ("13", None, None)


def test_strobj_D1S1656():
    marker = STRMarkerObject(
        "D1S1656",
        "TTCAGAGAAATAGAATCACTAGGGAACCAAATATATATACATACAATTAAACACACACACACCTATCTATCTATCTAT"
        "CTATCTATCTATCTATCTATCTATCTATCTA",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.uas_sequence == "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGGTGTGTGTGTG"
    assert marker.forward_sequence == "CACACACACACCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA"
    assert marker.convert == "CA [CACA]2 CCTA [TCTA]11"
    assert marker.convert_uas == "[TAGA]11 TAGG [TGTG]2 TG"
    assert marker.canonical == 12, " "
    assert marker.designation == ("11", "1", "0")


def test_strobj_D5S818():
    marker = STRMarkerObject(
        "D5S818",
        "TATTTATACCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTTCAAAAT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.uas_sequence == "AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAG"
    assert marker.forward_sequence == "CTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCT"
    assert marker.convert == "CTCT [ATCT]12"
    assert marker.convert_uas == "[AGAT]12 AGAG"
    assert marker.canonical == 12, " "
    assert marker.designation == ("12", None, None)


def test_strobj_D16S539():
    marker = STRMarkerObject(
        "D16S539",
        "TCCTCTTCCCTAGATCAATACAGACAGACAGACAGGTGGATAGATAGATAGATTGATTGATAGATAGATAGATAGATA"
        "TCATTGAAAGACAAAACAGAGATGGATGATAGATAC",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.uas_sequence == "GATAGATAGATAGATTGATTGATAGATAGATAGATAGATA"
    assert marker.forward_sequence == "GATAGATAGATAGATTGATTGATAGATAGATAGATAGATA"
    assert marker.flank_5p == "TC CTCT T CCCT AGAT CAAT [ACAG]4 GTG"
    assert marker.flank_3p == "TCAT TGAA AGAC AAA A CAGA [GATG]2 ATA GA T AC"
    assert marker.convert == "[GATA]3 [GATT]2 [GATA]5"


def test_strobj_D7S820():
    marker = STRMarkerObject(
        "D7S820",
        "TATTTAGTGAGATAAAAAAAAAACTATCAATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCGTT"
        "AGTTCGTTCTAAACTAT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.uas_sequence == "GATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGATTGATAGTTTT"
    assert marker.forward_sequence == "AAAACTATCAATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC"
    assert marker.flank_5p == "T ATTT AGTG AGAT AAAAAA"
    assert marker.flank_3p == "GTTA [GTTC]2 TAAA CTAT"
    assert marker.convert == "A AAAC TATC AATC TGTC [TATC]10"


def test_strobj_D3S1358():
    sequence = "TCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTAACTAACTATCTATCTA"
    marker = STRMarkerObject("D3S1358", sequence, "uas", kit="forenseq")
    assert marker.forward_sequence == sequence
    assert marker.uas_sequence == sequence
    assert marker.convert == "TCTA [TCTG]2 [TCTA]9 [ACTA]2 [TCTA]2"


def test_strobj_D19S433_newformat():
    marker = STRMarkerObject(
        "D19S433",
        "AATAAAAATCTTCTCTCTTTCTTCCTCTCTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTACCT"
        "TCCTTCCTTCAACAGAATCTTATTCTGTTGCCCAGGCTGGAGTCCAGTGTTACAATTATAGCT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "CT CTCT TTCT TCCT CTCT [CCTT]12 CCTA [CCTT]3"


def test_strobj_D21S11_newformat():
    marker = STRMarkerObject(
        "D21S11",
        "AAATATGTGAGTCAATTCCCCAAGTGAATTGCCTTCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCT"
        "ATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCT"
        "ACTATCTAT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == ("[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]11 TA")


def test_strobj_FGA_newformat():
    marker = STRMarkerObject(
        "FGA", "CCAGCAAAAAAGAAAGAAAGAGAAAAAAGAAAGAAAGAAA", "straitrazor", kit="forenseq"
    )
    assert marker.convert == "AAAA [AGAA]3 A"


def test_strobj_DYS643_foren():
    marker = STRMarkerObject(
        "DYS643",
        "TGATTTTTGCAGGTGTTCACTGCAAGCCATGCCTGGTTAAACTACTGTGCCTTTTCTTTTCTTTTCTTTTCTTTTCTT"
        "TTCTTTTCTTTTCTTTTCTTTTCTTTCTTTTTAAAACTT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[CTTTT]10 CTTTC TTTT"
    assert str(marker.canonical) == "10"
    assert marker.designation == ("10", None, None)
    assert marker.flank_3p == "TAAAA CTT"


def test_strobj_DYS635_foren():
    marker = STRMarkerObject(
        "DYS635",
        "ATCAATCAATGAATGGATAAAGAAAATGTGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAG"
        "ATACATACATAGATAGATACATACATAGATAGATAGATAGAGATTCTATGCAAAGTGAGAAGCCA",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[TAGA]12 [TACA]2 [TAGA]2 [TACA]2 [TAGA]4"
    assert str(marker.canonical) == "22"
    assert marker.designation == ("12", None, None)
    assert marker.flank_5p == "A [TCAA]2 TGAA TGGA TAAA GAAA ATGT GA"


def test_strobj_DYS612():
    marker = STRMarkerObject(
        "DYS612",
        "TTTCACACAGGTTCAGAGGTTTGCCTCCTCCTCCTCCTCTTTCTTCTTCTTCTCCTTCTTCTTCTTCTTCTTCTTCTT"
        "CTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTGTCACTTTTCCAAATTATTTTCTTTT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[CCT]5 CTT [TCT]4 CCT [TCT]24"
    assert marker.canonical == 29
    assert marker.designation == (24, 4, None)
    assert marker.flank_3p == "G TCA CTT TTC CAA [ATT]2 TTC TTT T"


def test_strobj_DYS576_foren():
    marker = STRMarkerObject(
        "DYS576",
        "AAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGCCAA"
        "GACAAATACGCTTATTACTCCCATCTCCT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[AAAG]17"
    assert str(marker.canonical) == "17"
    assert marker.designation == ("17", None, None)
    assert marker.flank_3p == "AAAA AGCC AAGA CAAA TACG CTTA TTAC TCCC ATCT CCT"


def test_strobj_DYS549_foren():
    marker = STRMarkerObject(
        "DYS549",
        "TAATAAGGTAGACATAGCAATTAGGTAGGTAAAGAGGAAGATGATAGATGATTAGAAAGATGATAGATAGATAGATAG"
        "ATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAAAATCTACATAAACAAAATCACAAATGGAAAAGGGGACATTACCA",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[GATA]13"
    assert str(marker.canonical) == "13"
    assert marker.designation == ("13", None, None)
    assert marker.flank_5p == (
        "TA ATAA GGTA GACA TAGC AATT [AGGT]2 AAAG AGGA AGAT GATA GATG ATTA GAAA GAT"
    )


def test_strobj_DYS533():
    marker = STRMarkerObject(
        "DYS533",
        "TAACTATATAACTATGTATTATCTATCAATCTTCTACCTATCATCTTTCTAGCTAGCTATCATCTATCTATCTATCTA"
        "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCATCTTCTATTGTTT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[TATC]13"
    assert str(marker.canonical) == "13"
    assert marker.designation == ("13", None, None)
    assert marker.flank_3p == "ATCT ATCA TCTT CTAT TGTT T"


def test_strobj_DYS522():
    marker = STRMarkerObject(
        "DYS522",
        "AGTTAGAGAGAGAGATGATGGATAGATAAATAGATAGATGATAGATGAATAGATAGGCGGGTAATAGATTTTATATAG"
        "ATAGATGATAGCTAGATAATGGATAGACATAGGTGACAGATGATAAATACATAGATAAATAGATGATAGATAGATAGATAGATAGATA"
        "GATAGATAGATAGATAGATAGATAGACAGATGTCCACCATGAGGTTC",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "ATA GATG [ATAG]12"
    assert marker.canonical == 12
    assert marker.designation == ("12", None, None)


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert, kit",
    [
        (
            "AAGGTGATAGATATACAGATAGATAGATACATAGGTGGAGACAGATAGATGATAAATAGAAGATAGATAGATAGATAGATAGATAGAT"
            "AGATAGATAGATAGATAGAAAGTATAAGTAAAGAGATGAT",
            "[GATA]2 TACA [GATA]3 CATA GGTG GAGA CAGA " "TAGA TGAT AAAT AGAA [GATA]11",
            "11",
            "11",
            None,
            None,
            "forenseq",
        ),
        (
            "AAATAGAAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAGTATAAGTAAAGAGATGATGG",
            "AAAT AGAA [GATA]12",
            "12",
            "12",
            None,
            None,
            "powerseq",
        ),
    ],
)
def test_strobj_DYS439(sequence, bracketed, conc, lus, sec, tert, kit):
    marker = STRMarkerObject("DYS439", sequence, "straitrazor", kit=kit)
    assert marker.convert == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)


def test_strobj_DYS437_foren():
    marker = STRMarkerObject(
        "DYS437",
        "ATGCCCATCCGGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGTCTGTCTGTCTATCTATCTATCTATC"
        "ATCTATCATCTGTGAATGATGTCTATCTACTTATCTATGAATGATATTTATCTGTGGTTATCTATCTATCTATATCATCTGTGAATGA"
        "CAGGGTCTTCCTCTG",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[TCTA]9 [TCTG]3 [TCTA]4"
    assert str(marker.canonical) == "16"
    assert marker.designation == ("9", None, None)
    assert marker.flank_3p == (
        "TCA TCTA TCAT CTGT GAAT GATG [TCTA]2 CTTA TCTA TGAA TGAT ATTT ATCT GTGG TTAT [CTAT]3 "
        "ATCA TCTG TGAA TGAC AGGG TCTT CCTC TG"
    )


def test_strobj_DYS392_foren():
    marker = STRMarkerObject(
        "DYS392",
        "TTAAACCTACCAATCCCATTCCTTAGTAAATAATAATAATAATAATAATAATAATAATAATAATAATAAATAAATGGT"
        "GATACAAGAAAAAAATTTGTTTTCCTTCTTGGCTTTTAAATAACAAACACTTGAAATCAAATTAGTTGTTTTTAAAAGCTAGATTAAT"
        "GAAGAA",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[ATA]13"
    assert str(marker.canonical) == "13"
    assert marker.designation == ("13", None, None)
    assert marker.flank_3p == (
        "AAT AAA TGG TGA TAC AAG [AAA]2 ATT TGT TTT CCT TCT TGG CTT TTA AAT AAC AAA CAC TTG AAA "
        "TCA AAT TAG TTG TTT TTA AAA GCT AGA TTA ATG AAG AA"
    )


def test_strobj_DYS391_foren():
    marker = STRMarkerObject(
        "DYS391",
        "ATATCTGTCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGCCTATCT"
        "GCCTGCCTACCTATCCCTCTAT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[TCTG]3 [TCTA]13 TCTG"
    assert str(marker.canonical) == "13"
    assert marker.designation == ("13", None, None)
    assert marker.flank_3p == "CCTA TCT [GCCT]2 ACCT ATCC CTCT AT"


def test_strobj_DYS19_foren():
    marker = STRMarkerObject(
        "DYS19",
        "TGGTCAATCTCTGCACCTGGAAATAGTGGCTGGGGCACCAGGAGTAATACTTCGGGCCATGGCCATGTAGTGAGGACAA"
        "GGAGTCCATCTGGGTTAAGGAGAGTGTCACTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTATCTATC"
        "TATCTA",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[TCTA]11 CCTA [TCTA]3"
    assert str(marker.canonical) == "14"
    assert marker.designation == ("11", None, None)
    assert marker.flank_3p == ""


def test_strobj_HPRTB():
    marker = STRMarkerObject(
        "HPRTB",
        "CTAGAACTTATCTTCTGTAAATCTGTCTCTATTTCCATCTCTGTCTCCATCTTTGTCTCTATCTCTATCTGTCTATCTCTATCTATCT"
        "ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAAAGCAAATTCATGCCCTTCTCCTATTT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[ATCT]12"
    assert marker.canonical == 12
    assert marker.designation == ("12", None, None)
    assert marker.flank_5p == (
        "CTAG AACT TATC TTCT GTAA ATCT GTCT CTAT TTCC ATCT CTGT CTCC ATCT TTGT CTCT ATCT CTAT "
        "CTGT CTAT C TCT"
    )


def test_strobj_DXS8378():
    marker = STRMarkerObject(
        "DXS8378",
        "AGTGAGCTGAGATGGTGCCACTGAACTCCAGCCTGGGCGACAAGAGCGAAACTCCAACTCAAAAAATAAATAAATAAAATATAGATAG"
        "ATAGATAGATAGATAGATAGATAGATAGATAGATAGTGACCTGCCAGGAGCAGGGGACCACCGGGTTGCCTAAGGAGGGGTGAACTGT"
        "CCCAGGATGGAAATGAAACA",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[ATAG]11"
    assert marker.canonical == 11
    assert marker.designation == ("11", None, None)
    assert marker.flank_5p == (
        "AGTG AGCT GAGA TGGT GCCA C TGAA CTCC AGCC TGGG CGAC AAGA GC G AAA CTCC AACT CAAA [AAAT]3"
        " AAAA T"
    )
    assert marker.flank_3p == (
        "TGAC CTGC CAGG AGCA GGGG ACCA CC G GGTT GCCT AAGG AGGG GTGA ACTG TCCC AGGA TGGA AATG "
        "AAAC A"
    )


def test_strobj_DXS7132():
    marker = STRMarkerObject(
        "DXS7132",
        "TCCAGAGAAACAGAACCAATAGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGTCAGATAGA"
        "TGAGAGGGGATTTATTATGAAAATGGGCTCACACTATTAAGGAGGCTAAGAAGTTCCACAGTAT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == "[TAGA]13"
    assert marker.canonical == 13
    assert marker.designation == ("13", None, None)
    assert marker.flank_3p == (
        "CAGT C [AGAT]2 GAGA GGGG ATTT ATTA TGAA AAT G GGCT CACA CTAT TAAG GAGG CTAA GAAG TTCC "
        "ACAG TAT"
    )


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert",
    [
        (
            "AAGAAAGAAAGAGAAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAA"
            "AAGAGAATAGAAAAGAAGAGAAGAGAAAAGAGAAAAGAAAAAAGAAAAGAAA",
            "[AAGA]3 gaaagga [AAGA]17 AAAG "
            "AGAA TAGA AAAG AAGA GAAG AGAA AAGA GAAA AGAA AAAA GAAA AGAA A",
            "21",
            "17",
            "0",
            None,
        ),
        (
            "AAGAAAGAAAGAGAAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAA"
            "AGAAAGGAAGAAAGAAAGAAAGGAAGAAAAGAGAATAGAAAAGAAGAGAAGAGAAAAGAGAAAAGAAAAAAGAAAAGAAA",
            "[AAGA]3 gaaagga [AAGA]18 AAGG [AAGA]3 AAGG AAGA AAAG AGAA TAGA AAAG AAGA GAAG AGAA AAGA "
            "GAAA AGAA AAAA GAAA AGAA A",
            "28",
            "18",
            "1",
            None,
        ),
    ],
)
def test_strobj_DXS10135(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject("DXS10135", sequence, "straitrazor", kit="forenseq")
    assert marker.convert == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert, flank_5p, flank_3p",
    [
        (
            "TGTGTGTGCATGCATACACACACAGAGAGAGAGAGAGAGAAAAAGAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAG"
            "AAAGAAAGAAAGAAAGAAAGGAAGAAAGAAAGGAAGAAAATAGAACAAATCAGCTTATATTCAGTATTTTTTTAGTATTTTCTGTGTC"
            "AGCTC",
            "[AAGA]15 AAGG [AAGA]2",
            "18",
            "15",
            "0",
            None,
            "[TGTG]2 CATG CATA CACA CA C " "[AGAG]4 AAAA AG",
            "AAGG A AAGG AAGA AAAT AGAA CAAA TCAG CTTA TATT CAGT ATTT TTTT AGTA "
            "TTTT CTGT GTCA G TC",
        ),
        (
            "TGTGTGTGCATGCATACACACACAGAGAGAGAGAGAGAAAAAGAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAA"
            "AGAAAGAAAGAAAGAAAGGAAGAAAGAAAGGAAGAAAATAGAACAAATCAGCTTATATTCAGTATTTTTTTAGTATTTTCTGTGTCAG"
            "CTC",
            "GA [AAGA]14 AAGG [AAGA]2",
            "17.2",
            "14",
            "0",
            None,
            "[TGTG]2 CATG CATA CACA CA C " "[AGAG]3 AGAA AAAG AA",
            "AAGG A AAGG AAGA AAAT AGAA CAAA TCAG CTTA TATT CAGT ATTT TTTT "
            "AGTA TTTT CTGT GTCA G TC",
        ),
    ],
)
def test_strobj_DXS10074(sequence, bracketed, conc, lus, sec, tert, flank_5p, flank_3p):
    marker = STRMarkerObject("DXS10074", sequence, "straitrazor", kit="forenseq")
    assert marker.convert == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5p
    assert marker.flank_3p == flank_3p


def test_strobj_Y_GATA_H4():
    marker = STRMarkerObject(
        "Y-GATA-H4",
        "CTATCTATCTATCTATTCATCCATCTAATCTATCCATTCTATCTATCTATCTATCTATCTATCTATCTATCTATC"
        "TATCTACCTACCTACCTATCTATCTATAGATCTATCTATCTATCT",
        "straitrazor",
        kit="forenseq",
    )
    assert marker.convert == (
        "C [TATC]3 TATT [CATC]2 TAAT CTAT CCAT [TCTA]11 [CCTA]3 [TCTA]2 TAGA [TCTA]3 TCT"
    )
    assert str(marker.canonical) == "11"
    assert marker.designation == ("11", None, None)


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert, flank_5p, kit",
    [
        (
            "AGAAATGGATGACAGTAAAATGAAAACATTGCAATGTGTATACTCAGAAACAAGGAAAGATAGATAGATGATAGATAGATAGATAGAC"
            "AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACAGACAGACAGACAGATAGAT",
            "[TAGA]4 CAGA [TAGA]11 [CAGA]8 TAGA T",
            "24",
            "11",
            "8",
            "0",
            "AG AAAT GGAT GACA GTAA "
            "AATG AAAA CATT GCAA TGTG TATA CTCA GAAA CAAG GAAA [GATA]2 GATG A",
            "forenseq",
        ),
        (
            "AACAAGGAAAGATAGATAGATGATAGATAGATAGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAC"
            "AGACAGACAGACAGACAGACAGACAGACAGATAGATAGAATATATTATGGGGTACCAAAATGCAGGGCCCAAAAATGTGTAAAATATA"
            "TGTGT",
            "[TAGA]4 CAGA [TAGA]11 [CAGA]8 [TAGA]2",
            "24",
            "11",
            "8",
            "0",
            "AA CAAG GAAA [GATA]2 GATG A",
            "powerseq",
        ),
    ],
)
def test_strobj_DYS390(sequence, bracketed, conc, lus, sec, tert, flank_5p, kit):
    marker = STRMarkerObject("DYS390", sequence, "straitrazor", kit=kit)
    assert marker.convert == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5p


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert, kit",
    [
        (
            "TTTCTTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTATTTCTTTCCCTTCCTTCCTT"
            "CCTTCCTTCCTTTCTTTCTCTTTCCTCTTTCTCTTTCTTCTCTTTCTTTCTTTTTCTCTTTTTCTCTTTCTTTCTTTTTTACTTTCTT"
            "TCTCCTTCCTTCCTTCCTTTCTGAATTTCATTTCTTTTCTTT",
            "TTTC TTTT TCTC [TTTC]13 TTTA [TTTC]2 "
            "[CCTT]6 TCTT TCTC TTTC CTCT TTCT CTTT CTTC [TCTT]3 TTTC TCTT TTTC [TCTT]3 TTTT ACTT TCTT"
            " TCTC [CTTC]3 CTT",
            "16",
            "13",
            None,
            None,
            "forenseq",
        ),
        (
            "TTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCCCTTCCTTCCTTCCTTCCTTCCTTTCTT"
            "TCTCTTTCCTCTTTCTCTTTCTTCTCTTTCTTTCTTTTTCTCTTTTTCTCTTTCTTTCTTTTTTACTTTCTTTCTCCTTCCTTCCTTC"
            "CTTTCTGAATTTCATTTCTTTTCTTT",
            "TTTT TCTC [TTTC]13 [CCTT]6 TCTT TCTC TTTC CTCT TTCT CTTT "
            "CTTC [TCTT]3 TTTC TCTT TTTC [TCTT]3 TTTT ACTT TCTT TCTC [CTTC]3 CTT",
            "13",
            "13",
            None,
            None,
            "powerseq",
        ),
    ],
)
def test_strobj_DYS385(sequence, bracketed, conc, lus, sec, tert, kit):
    marker = STRMarkerObject("DYS385A-B", sequence, "straitrazor", kit=kit)
    assert marker.convert == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert",
    [
        (
            "GAGATAGAGACATGGATAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATATA"
            "GAGATAGAGAGATAGAGATAGAGATAGATAGATAGAGAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATA"
            "GAGA",
            "[AGAGAT]11 [ATAGAG]2 [AGATAG]3 ATAGAT AGAGAA [AGAGAT]8",
            19,
            11,
            8,
            None,
        ),
        (
            "GAGATAGAGACATGGATAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGA"
            "GATATAGAGATAGAGAGATAGAGATAGAGATAGATAGATAGAGAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGAGAGAGATA"
            "GAGATAGAGA",
            "[AGAGAT]12 [ATAGAG]2 [AGATAG]3 ATAGAT AGAGAA [AGAGAT]5 AGAGAG [AGAGAT]2",
            20,
            12,
            2,
            None,
        ),
    ],
)
def test_strobj_DYS448(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject("DYS448", sequence, "straitrazor", kit="forenseq")
    assert marker.convert == bracketed
    assert marker.canonical == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert",
    [
        (
            "TGCATGCACATACACATAACTAGATAGACTGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGA"
            "CAGATAGATA",
            "[TAGA]2 CTGA CAGA [TAGA]10 [CAGA]4 TAGA",
            18,
            10,
            2,
            None,
        ),
        (
            "TGCATGCACATACACATATCTAGATCGACTGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGA"
            "CAGATAGATA",
            "TAGA TCGA CTGA CAGA [TAGA]10 [CAGA]4 TAGA",
            18,
            10,
            1,
            None,
        ),
    ],
)
def test_strobj_DXS10103(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject("DXS10103", sequence, "straitrazor", kit="forenseq")
    assert marker.convert == bracketed
    assert marker.canonical == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize(
    "sequence, bracketed, conc, lus, sec, tert",
    [
        (
            "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGATACATAGATAATACAGATGAGAGTTGGATACA"
            "GAAGTAGGTATAATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACACA"
            "CACATAGATAATACAGAT",
            "[TAGA]11 [CAGA]3 TACA TAGA TAAT ACAG ATGA GAGT TGGA TACA GAAG TAGG"
            " TATA ATGA [TAGA]13 [CAGA]4",
            31,
            13,
            11,
            None,
        )
    ],
)
def test_strobj_DYS389II(sequence, bracketed, conc, lus, sec, tert):
    marker = STRMarkerObject("DYS389II", sequence, "straitrazor", kit="forenseq")
    assert marker.convert == bracketed
    assert marker.canonical == conc
    assert marker.designation == (lus, sec, tert)


@pytest.mark.parametrize(
    "sequence, bracketed, lus, sec, tert, flank_5, flank_3, kit",
    [
        (
            "GTCTCAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGAGAGAGGAAAGAAAGAGAAAAAG"
            "AAAAGAAATAGTAGCAACTGTTATTGTAAGA",
            "[AGAA]14 AAAG AGAG AG",
            "14",
            "1",
            None,
            "G TCTC",
            "GA [AAGA]2 GAAA AAGA AAAG AAAT AGTA GCA A CTGT TATT GT AAGA",
            "forenseq",
        ),
        (
            "AGGCTGCAGTGAGCCATGTTCATGCCACTGCACTTCACTCTGAGTGACAAATTGAGACCTTGTCTCAGAAAGAAAGAAAGAAAGAAAG"
            "AAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGAGAGAGGAAAGAAAGAGAAAAAGAAAAGAAATAG"
            "TAGCAACTGTTATTGTA",
            "[AGAA]18 AAAG AGAG AG",
            "18",
            "1",
            None,
            "A GGCT GCAG TGAG CCAT GTTC "
            "ATGC CACT GCAC TTCA CTCT GAGT GACA AATT GAGA CCTT G TCTC",
            "GA [AAGA]2 GAAA AAGA AAAG " "AAAT AGTA GCA A ACTG TTAT TGTA",
            "powerseq",
        ),
    ],
)
def test_strobj_D18S51(sequence, bracketed, lus, sec, tert, flank_5, flank_3, kit):
    marker = STRMarkerObject("D18S51", sequence, "straitrazor", kit=kit)
    assert marker.convert == bracketed
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5
    assert marker.flank_3p == flank_3


@pytest.mark.parametrize(
    "locus, sequence, bracketed, conc, lus, sec, tert, flank_5, flank_3",
    [
        (
            "CSF1PO",
            "CTAAGTACTTCCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAATCTATCTATCTTCTA"
            "TCTATGAAGGCAGTTACTGTTAATATCTTCATTTTACAGGTAGGAAAACTGAGACACAGGGTGGTTAGCAACCTGCTAGTCCTTGGCA"
            "GACTCAG",
            "[ATCT]12",
            "12",
            "12",
            "0",
            None,
            "CTA AGTA CT TCCT",
            "A [ATCT]3 T [CTAT]2 "
            "GAAG GCAG TTAC TGTT AATA TCTT CATT TTAC AGGT AGGA AAAC TGAG ACAC AGGG TGGT TAG CA ACCT "
            "GCTA GTCC TTGG CAGA CTCA G",
        ),
        (
            "D10S1248",
            "CCCCAGGACCAATCTGGTCACAAACATATTAATGAATTGAACAAATGAGTGAGTGGAAGGAAGGAAGGAAGGAAGG"
            "AAGGAAGGAAGGAAGGAAGGAAGGAAGGAA",
            "[GGAA]13",
            "13",
            "13",
            None,
            None,
            "CCCC AGGA CCAA " "TCTG GTCA CAAA CATA TTAA TGAA TT GAAC AAAT [GAGT]2",
            "",
        ),
        (
            "D12S391",
            "CAATGGATGCATAGGTAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGAC"
            "AGATGAGAGGG",
            "[AGAT]12 [AGAC]5 AGAT",
            "18",
            "12",
            "5",
            "0",
            "CAAT GGAT GCAT AGGT",
            "GAGA GGG",
        ),
        (
            "D13S317",
            "TTCTTTAGTGGGCATCCGTGACTCTCTGGACTCTGACCCATCTAACGCCTATCTGTATTTACAAATACATTATCTA"
            "TCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTCTTTTTGGGCTGCCTATGGCTCAACCCA"
            "AGTTGAAGGAGGAGATTT",
            "[TATC]9 [AATC]2 [ATCT]3 TTCT GTCT GTC",
            "9",
            "9",
            "3",
            "1",
            "TT "
            "CTT TAGT GGGC ATCC G TGAC TCTCT GGAC TC TGAC CCAT CTAA C G CCT ATCT GTAT TTAC AAAT ACAT",
            "TTTT TGGG CTGC CTAT GGCT CAAC CCAA GTTG AAGG AGGA GATT T",
        ),
        (
            "D16S539",
            "GTGCACAAATCTAAATGCAGAAAAGCACTGAAAGAAGAATCCAGAAAACCACAGTTCCCATTTTTATATGGGAGCAA"
            "ACAAAGGCAGATCCCAAGCTCTTCCTCTTCCCTAGATCAATACAGACAGACAGACAGGTGGATAGATAGATAGATAGATAGATAGATA"
            "GATAGATATCAT",
            "[GATA]9",
            "9",
            "9",
            "0",
            None,
            "GT GCAC AAAT CTAA ATGC AGAA AAGC ACTG "
            "AAAG AAGA ATCC AG AAAA CCAC AGTT CCCA TTTT TATA TGGG AG [CAAA]2 GGCA GATC CCAA G CTCT TC"
            " CTCT T CCCT AGAT CAAT [ACAG]4 GTG",
            "TCAT",
        ),
        (
            "D18S51",
            "AGGCTGCAGTGAGCCATGTTCATGCCACTGCACTTCACTCTGAGTGACAAATTGAGACCTTGTCTCAGAAAGAAAGAAAGAAAG"
            "AAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGAGAGAGGAAAGAAAGAGAAAAAGAAAAGAAATAGTAGCAACTGTTA"
            "TTGTA",
            "[AGAA]13 AAAG AGAG AG",
            "13",
            "13",
            "1",
            None,
            "A GGCT GCAG TGAG CCAT GTTC ATGC CACT GCAC TTCA CTCT GAGT GACA AATT GAGA CCTT G TCTC",
            "GA [AAGA]2 GAAA AAGA AAAG AAAT AGTA GCA A ACTG TTAT TGTA",
        ),
        (
            "D19S433",
            "AAGTTCTTTAGCAGTGATTTCTGATATTTTGGTGCACCCATTACCCGAATAAAAATCTTCTCTCTTTCTTCCTCTCT"
            "CCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTACCTTCTTTCCTTCAACAGAATCTTATTCTGTTGCCCAGGC"
            "TGGAGTGCA",
            "CT CTCT TTCT TCCT CTCT [CCTT]11 CCTA CCTT CTTT CCTT",
            "13",
            "11",
            "1",
            "0",
            "AAG TTCT TTAG CAGT GATT TCTG ATAT TTTG GTGC ACCC ATTA CCCG AATA AAAA TCTT",
            "CAAC AGAA " "TCTT ATTC TGTT GCCC AGGC TGGA GTGC A",
        ),
        (
            "D1S1656",
            "GAAATAGAATCACTAGGGAACCAAATATATATACATACAATTAAACACACACACATCTATCTATCTATCTATCTATC"
            "TATCTATCTATCTATCTATCTATCTATCTACATCACACAGTTGACCCTTGA",
            "CA [CACA]2 [TCTA]13",
            "13",
            "13",
            "0",
            "0",
            "G AAAT AGAA TCAC TAGG GAAC CAAA [TATA]2 CATA CAAT TAAA",
            "CATC AC ACA GTTG " "ACCC TTGA",
        ),
        (
            "D21S11",
            "TGAATTGCCTTCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATC"
            "TATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCTATCGTCTATCTATCCAGTCT"
            "ATCTACCTCCTATTAGTCT",
            "[TCTA]5 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]11 " "TA TCTA",
            "31.2",
            "11",
            "5",
            "6",
            "TG AATT GCCT",
            "TC G [TCTA]2 T CCAG [TCTA]2 CCTC " "CTAT T AGTC T",
        ),
        (
            "D22S1045",
            "CCTTCTTATAGCTGCTATGGGGGCTAGATTTTCCCCGATGATAGTAGTCTCATTATTATTATTATTATTATTATTA"
            "TTATTATTATTATTACTATTATTGTTATAAAAATATTG",
            "[ATT]13 ACT [ATT]2",
            "16",
            "13",
            None,
            None,
            "CCT TCT TAT AGC TGC TAT GGG GGC TAG ATT TTC CCC [GAT]2 [AGT]2 CTC",
            "GTT ATA AAA ATA TTG",
        ),
        (
            "D2S1338",
            "CTAGCATGGTACCTGCAGGTGGCCCATAATCATGAGTTATTCAGTAAGTTAAAGGATTGCAGGAGGGAAGGAAGGAC"
            "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGCAGGCAGGCAGGCAGGCAGGCAGGCAAGGCCAAGCCAT"
            "TTCTGTTTCCAAATCCACTGGCTCCCTCCCACAGCT",
            "[GGAA]2 GGAC [GGAA]12 [GGCA]7",
            "22",
            "12",
            "1",
            "7",
            "CT AGCA TGGT ACCT GCAG GTGG CCCA TAAT C ATGA GTTA TTCA GTAA GTTA AAGG ATTG CAG GAG",
            "AGGC CAAG CCAT TT CTGT TTCC AAAT CCAC TGGC [TCCC]2 ACAG CT",
        ),
        (
            "D2S441",
            "TGCACCCAACATTCTAACAAAAGGCTGTAACAAGGGCTACAGGAATCATGAGCCAGGAACTGTGGCTCATCTATGAAA"
            "ACTTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATATCA",
            "[TCTA]10",
            "10",
            "10",
            "0",
            None,
            "TG CACC CAAC ATTC TAAC AAAA GGCT GTAA CAAG GGCT ACAG GAA T CATG AG CCAG G AACT GTGG CTCA"
            " TCTA TGAA AACT",
            "TATC A",
        ),
        (
            "D3S1358",
            "TGCCCACTTCTGCCCAGGGATCTATTTTTCTGTGGTGTGTATTCCCTGTGCCTTTGGGGGCATCTCTTATACTCATG"
            "AAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAT"
            "CTATGAGACAGGGTCTTGCTCTGTC",
            "TCTA [TCTG]3 [TCTA]13",
            "17",
            "13",
            "3",
            None,
            "TGCC CACT "
            "TCTG CCCA GGGA TCTA TTTT TCTG TGGT GTGT ATTC CCTG TGCC TTTG GGGG CATC TCTT ATAC TCAT "
            "GAAA TCAA CAGA GGCT TGCA TGTA",
            "TGAG ACAG GGTC TTGC TC TGTC",
        ),
        (
            "D5S818",
            "AACATTTGTATCTTTATCTGTATCCTTATTTATACCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC"
            "TATCTATCTTCAAAATATTACGTAAGGATACCAAAGAGGAAAATCACCCTTGTCACATACTTGCTATTAAAATATACTTTTATTAGTA"
            "CA",
            "CTCT [ATCT]12",
            "12",
            "12",
            None,
            None,
            "AA CATT TGTA TCTT TATC TGTA TCCT T ATTT " "ATAC",
            "TCAA AAT ATTA CG TAAG GATA CCAA AGAG GAAA ATCA CCCT TGTC ACAT ACTT GCTA T TAAA "
            "ATAT ACTT TTAT TAGT ACA",
        ),
        (
            "D7S820",
            "AGAATTGCACCAAATATTGGTAATTAAATGTTTACTATAGACTATTTAGTGAGATAAAAAAAAACTATCAATCTGTCT"
            "ATCTATCTATCTATCTATCTATCTATCTATCGTTAGTTCGTTCTAAACTATGACAAGTGTTCTATCATACCCTTTATATATATTAACC"
            "TTAAAATAACTC",
            "AAAC TATC AATC TGTC [TATC]8",
            "8",
            "8",
            "1",
            0,
            "AGAA TTGC ACCA A " "ATAT TGGT AATT AAAT GTTT ACTA T AGAC T ATTT AGTG AGAT AAAAAA",
            "GTTA [GTTC]2 TAAA CTAT "
            "GACA AGTG TTCT ATCA TACC CTTT [ATAT]2 TAAC CTTA AAAT AACT C",
        ),
        (
            "FGA",
            "GTCTGAAATCGAAAATATGGTTATTGAAGTAGCTGCTGAGTGATTTGTCTGTAATTGCCAGCAAAAAAGAAAGGAAGAAAG"
            "GAAGGAAGGAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAGAAAAAAGAAAGAAAGAAA",
            "[GGAA]2 GGAG [AAAG]12 AGAA AAAA [GAAA]3",
            "20",
            "12",
            "3",
            "0",
            "G TCTG AAAT CGAA AATA "
            "TGGT TATT GAAG TAGC TGCT GAGT GATT TGTC TGTA ATTG CCA GCAA AAAA GAAA GGAA GAAA",
            "",
        ),
        (
            "PENTA D",
            "GAGCCATGATCACACCACTACACTCCAGCCTAGGTGACAGAGCAAGACACCATCTCAAGAAAGAAAAAAAAGAAAGAA"
            "AAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAAACGAAGGGGAAAAAAAGAGAATCATAAACA"
            "TAAATGTAAAATTTCTCAAAAAAATCGTTA",
            "AAAAG [AAAGA]12",
            "12",
            "12",
            None,
            None,
            "GA GCCAT GATCA CACCA CTACA CTCCA GCCTA GGTGA CAGAG CAAGA CACCA TCTCA AGAAA GAAAA",
            "AAAAA CGAA GGGGA AAAAA AGAGA ATCAT AAACA TAAAT GTAAA ATTTC TCAAA AAAAT CGTTA",
        ),
        (
            "PENTA E",
            "TAATGATTACATAACATACATGTGTGTAAAGTGCTTAGTATCATGATTGATACATGGAAAGAATTCTCTTATTTGGGT"
            "TATTAATTGAGAAAACTCCTTACAATTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTGAGAC",
            "[TCTTT]7",
            "7",
            "7",
            None,
            None,
            "TA ATGAT TACAT AACAT ACATG TGTGT AAAGT GCTTA GTATC ATGAT TGATA "
            "CATGG AAAGA ATTCT CTTAT TTGGG TTATT AATTG AGAA AACTC CTTAC AATTT",
            "GAGAC",
        ),
        (
            "TH01",
            "CTCCATGGTGAATGAATGAATGAATGAATGAATGAGGGAAATAAGGGAGGAACAGGCCAATGGGAATCACCCCAGAGCCC"
            "AGATACCCTTTGAATTTTGCCCCCTATTTGCCCAGGACCCCCCACCATGAGCTGCTGCTAGAGCCTGGGAAGGGCCTTGGGGCTGCCT"
            "CCCCAAGCAGGCAGGCTGGTTGGGGTGC",
            "[AATG]6",
            "6",
            "6",
            None,
            None,
            "CT CCAT GGTG",
            "AGGG AAAT AAGG GAGG AACA GGCC AATG GGAA TCAC CCCA GAGC CCAG ATAC CCTT TGAA TTTT GCCC "
            "CCTA TTTG CCCA GGAC CCCC CACC ATGA GCTG CTGC TAGA GCCT GGGA AGGG CCTT GGGG C TGCC TCCC "
            "CAAG [CAGG]2 CTGG TTGG GGTG C",
        ),
        (
            "TPOX",
            "CACTGGCCTGTGGGTCCCCCCATAGATCGTAAGCCCAGGAGGAAGGGCTGTGTTTCAGGGCTGTGATCACTAGCACCCAG"
            "AACCGTCGACTGGCACAGAACAGGCACTTAGGGAACCCTCACTGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATG"
            "TTTGGGCAAATAAACGCT",
            "[AATG]11",
            "11",
            "11",
            None,
            None,
            "CAC TGGC CTGT GGGT C CCCC CAT AGAT"
            " CG TAAG CCCA GGAG GAAG GGCT GTGT TTCA GGGC TGTG ATCA CTA GC ACCC AGAA CCGT CGAC TGGC "
            "ACAG AACA GGCA CTTA GGGA ACCC TCAC TG",
            "TTTG G GCAA ATAA ACGC T",
        ),
        (
            "VWA",
            "GGATAGATGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGATAGATCAATC"
            "CAAGTCACATACTGATTATTCTTATCATCCACTAGGGCTTTCACATCTCAGCCAAGTCAACTTGGATCCTCTAGACCTGTTTCTTCTT"
            "CTGGAA",
            "TAGA TGGA [TAGA]12 [CAGA]3 TAGA",
            "16",
            "12",
            "3",
            "1",
            "GGA",
            "TCAA T CCAA GTCA CATA CTGA TTAT TCTT ATCA TCCA CTAG GGCT TTCA CATC TCAG CCAA GTCA ACTT "
            "GGAT CCTC TAGA CCTG TTTC TTCT TCTG GAA",
        ),
        (
            "DYS643",
            "ACCTCATGCTCTGTGATTTTTGCAGGTGTTCACTGCAAGCCATGCCTGGTTAAACTACTGTGCCTTTTCTTTTCTTTT"
            "CTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTCTTTTTAAAACTTTTTACTTCAGTAGAATTTTGGGGGG",
            "[CTTTT]10 CTTTC TTTT",
            "10",
            "10",
            None,
            None,
            "ACC TCATG CTCTG TGATT TTTGC AGGTG " "TTCAC TGCAA GCCAT GCCTG GTTAA ACTAC TGTGC",
            "TAAAA CTT TTTAC TTCAG TAGAA TTTTG GGGGG",
        ),
        (
            "DYS635",
            "CCCAAATATCCATCAATCAATGAATGGATAAAGAAAATGTGATAGATAGATAGATAGATAGATAGATAGATAGATAGA"
            "TAGATAGATACATACATAGATAGATACATACATAGATAGATAGATAGAGATT",
            "[TAGA]11 [TACA]2 [TAGA]2 [TACA]2 " "[TAGA]4",
            "21",
            "11",
            None,
            None,
            "CCC AAAT ATCC A [TCAA]2 TGAA TGGA TAAA GAAA ATGT GA",
            "GATT",
        ),
        (
            "DYS576",
            "AAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAA"
            "GCCAAGACAAATACGCTTATTACTCCCATCTCCTCCTTCATCTCCAGGAAATGAGAC",
            "[AAAG]18",
            "18",
            "18",
            None,
            None,
            "A",
            "AAAA AGCC AAGA CAAA TACG CTTA TTAC TCCC ATCT CCT CCTT CATC TCCA GGAA ATGA GAC",
        ),
        (
            "DYS570",
            "TAAAATGAATGATGACTAGGTAGAAATCCTGGCTGTGTCCTCCAAGTTCCTTTTCTTTCTTTCTTTCTTTCTTTCTTT"
            "CTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTT",
            "[TTTC]17",
            "17",
            "17",
            None,
            None,
            "TAA [AATG]2 ATGA CTAG GTAG AAAT CCTG GCTG TGTC CTCC AAGT TCCT",
            "TTTT T",
        ),
        (
            "DYS549",
            "GTAAAGAACTATAAAAAGATTAATACAACAAAAATTTGGTAATCTGAAATAATAAGGTAGACATAGCAATTAGGTAGG"
            "TAAAGAGGAAGATGATAGATGATTAGAAAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAA"
            "AAATC",
            "[GATA]13",
            "13",
            "13",
            None,
            None,
            "G TAAA GAAC TATA AAAA GATT AATA CAAC AAAA "
            "ATTT GGTA ATCT GAAA TA ATAA GGTA GACA TAGC AATT [AGGT]2 AAAG AGGA AGAT GATA GATG ATTA "
            "GAAA GAT",
            "GAAA AAAT C",
        ),
        (
            "DYS533",
            "CTAATATTTATCTATATCATTCTAATTATGTCTCTTCTAACTATATAACTATGTATTATCTATCAATCTTCTACCTAT"
            "CATCTTTCTAGCTAGCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCATCTTCTATT"
            "GTTTGGTTGAGTTAAGAACTGATCATGAATAAATACATTTCATTGGT",
            "[TATC]12",
            "12",
            "12",
            None,
            None,
            "C TAAT ATTT ATCT ATAT CATT CTAA TTAT GTCT CTTC TAAC TATA TAAC TATG TATT ATCT ATCA ATCT "
            "TCTA CCTA TCAT CTTT [CTAG]2 CTAT CATC",
            "ATCT ATCA TCTT CTAT TGTT T GGTT GAGT TAAG AACT " "GATC ATGA ATAA ATAC ATTT CATT GGT",
        ),
        (
            "DYS481",
            "TAAAAGGAATGTGGCTAACGCTGTTCAGCATGCTGCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTC"
            "TTCTTCTTCTTCTTCTTCTTCTTTTTTGAGTCTTG",
            "[CTT]22",
            "22",
            "22",
            None,
            None,
            "TA AAA GGA ATG" " TGG CTA ACG CTG TTC AGC ATG CTG",
            "TTT TGA GTC TTG",
        ),
        (
            "DYS458",
            "GAAAGAAAGAAAAGGAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA"
            "GAAAGAAAGGAGGGTGGGCGTGGTGGCTCATGCTTGTAATGCCAGAACTTTGGGAGGCCGAGGTGG",
            "[GAAA]3 AG GAAG "
            "[GAAA]17 GGAG GGTG GGCG TGGT GGCT CATG CTTG TAAT GCCA GAAC TTTG GGAG GCCG AGGT GG",
            "17",
            "17",
            None,
            None,
            None,
            None,
        ),
        (
            "DYS456",
            "AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATTCCATTAG"
            "TTCTGTCCCTCTAGAGAACCCTAATACATCAGTTTAAGAA",
            "[AGAT]17 ATTC CATT AGTT CTGT CCCT CTAG AGAA " "CCCT AATA CATC AGTT TAAG AA",
            "17",
            "17",
            None,
            None,
            None,
            None,
        ),
        (
            "DYS448",
            "AGACATGGATAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGAT"
            "ATAGAGATAGAGAGATAGAGATAGAGATAGATAGATAGAGAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAG"
            "ATAGAGAGGTAAAGATAGA",
            "[AGAGAT]11 [ATAGAG]2 [AGATAG]3 ATAGAT AGAGAA [AGAGAT]8",
            "19",
            11,
            8,
            None,
            "AGACAT GGATAA",
            "AGAGA GGTAAA GATAGA",
        ),
        (
            "DYS439",
            "AAATAGAAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAGTATAAGTAAAGAGATGA" "TGG",
            "AAAT AGAA [GATA]12",
            "12",
            "12",
            None,
            None,
            None,
            "GAAA GTAT AAGT AAAG AGAT GAT " "GG",
        ),
        (
            "DYS438",
            "CAGTATATTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTATTTGAAATGGAGTTTCACTCTTGT"
            "TGCCCAGG",
            "[TTTTC]9",
            "9",
            "9",
            None,
            None,
            "CA GTATA",
            "TATTT GAAAT GGAGT TTCAC TCTTG " "TTGCC CAGG",
        ),
        (
            "DYS437",
            "GCCCATCCGGTCTATCTATCTATCTATCTATCTATCTATCTATCTGTCTGTCTATCTATCTATCTATCATCTATCATC"
            "TGTGAATGATGTCTATCTACTTATCTATGAATGATATTTATCTGTGGTTATCTATCTATCTATA",
            "[TCTA]8 [TCTG]2 [TCTA]4",
            "14",
            "8",
            None,
            None,
            "GC CCAT CCGG",
            "TCA TCTA TCAT CTGT "
            "GAAT GATG [TCTA]2 CTTA TCTA TGAA TGAT ATTT ATCT GTGG TTAT [CTAT]3 A",
        ),
        (
            "DYS393",
            "AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATGTATGTCTTTTCTATGAGACATA",
            "[AGAT]13 [ATGT]2 CTTT TCTA TGAG ACAT A",
            "13",
            "13",
            None,
            None,
            None,
            None,
        ),
        (
            "DYS392",
            "TAAATAATAATAATAATAATAATAATAATAATAATAATAATAAATAAATGGTGATACAAGAAAAAAATTTGTTTTCCT"
            "TCTTGGCTTTTAAATAACAAACACTTGAAATCAAATTAG",
            "[ATA]13",
            "13",
            "13",
            None,
            None,
            "TAA",
            "AAT AAA TGG TGA TAC AAG [AAA]2 ATT TGT TTT CCT TCT TGG CTT TTA AAT AAC AAA CAC TTG AAA "
            "TCA AAT TAG",
        ),
        (
            "DYS391",
            "TGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGCCTATCTGCCTGCCTACCTATCCCTCTA"
            "TGGCAATTGCTTGCAACCAGGGAGATTTTA",
            "TG TCTG [TCTA]10 TCTG",
            "10",
            "10",
            None,
            None,
            None,
            "CCTA TCT [GCCT]2 ACCT ATCC CTCT AT GGCA ATTG CTTG CAAC CAGG GAGA TTTT A",
        ),
        (
            "DYS390",
            "AACAAGGAAAGATAGATAGATGATAGATAGATAGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAG"
            "ATAGATAGACAGACAGACAGACAGACAGACAGACAGACAGATAGATAGAATATATTATGGGGTACCAAAATGCAGGGCCCAAAAATGT"
            "GTAAAATATATGTGT",
            "[TAGA]4 CAGA [TAGA]11 [CAGA]8 [TAGA]2",
            "24",
            "11",
            "8",
            "0",
            "AA CAAG GAAA [GATA]2 GATG A",
            "ATAT ATTA TGGG GTAC CAAA ATGC AGGG CCCA AAAA TGTG TAAA " "ATAT ATGT GT",
        ),
        (
            "DYS389II",
            "TCATAGATAGATGATGGACTGCTAGATAAATAGATAGATTGATAGAGGGAGGGATAGATAGATAGATAGATAGATA"
            "GATAGATAGATAGATAGATAGACAGACAGACAGATACATAGATAATACAGATGAGAGTTGGATACAGAAGTAGGTATAATGATAGATA"
            "GATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACA",
            "[TAGA]11 [CAGA]3 "
            "TACA TAGA TAAT ACAG ATGA GAGT TGGA TACA GAAG TAGG TATA ATGA [TAGA]13 [CAGA]4",
            "31",
            13,
            11,
            None,
            "TCA [TAGA]2 TGAT GGAC TGCT AGAT AAAT [AGAT]2 TGAT AGAG GGAG GGA",
            "CA",
        ),
        (
            "DYS385A-B",
            "TTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCCCTTCCTTCCTTCCT"
            "TCCTTCCTTTCTTTCTCTTTCCTCTTTCTCTTTCTTCTCTTTCTTTCTTTTTCTCTTTTTCTCTTTCTTTCTTTTTTACTTTCTTTCT"
            "CCTTCCTTCCTTCCTTTCTGAATTTCATTTCTTTTCTTT",
            "TTTT TCTC [TTTC]13 [CCTT]6 TCTT TCTC TTTC "
            "CTCT TTCT CTTT CTTC [TCTT]3 TTTC TCTT TTTC [TCTT]3 TTTT ACTT TCTT TCTC [CTTC]3 CTT",
            "13",
            "13",
            None,
            None,
            None,
            "TCTG AATT TCAT TTCT TTTC TTT",
        ),
        (
            "DYS19",
            "TCTGGGTTAAGGAGAGTGTCACTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTATCTAT"
            "CTATCTAAAACACTATATATATATAACACTATATATATAATACTATATATATATTA",
            "[TCTA]11 CCTA [TCTA]3",
            "14",
            "11",
            None,
            None,
            "TC TGGG TTAA GGAG AGTG TCAC TATA",
            "AA ACAC [TATA]3 ACAC [TATA]2 TA " "ATAC [TATA]2 TATT A",
        ),
        (
            "Y-GATA-H4",
            "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTACCTACCTATCTATCTATAGATCTATCT"
            "ATCTATCTTAAATTTGGAAATTCTCCTCAGCATAACATTTTAATGATGATTCCTAGGATACAAGTGATGTGCTGAAAGTATCAATGTG"
            "TATCAGAAAACCAACATCTCTGCTTAGGTCTCT",
            "[TCTA]11 [CCTA]3 [TCTA]2 TAGA [TCTA]3 TCT",
            "11",
            "11",
            None,
            None,
            None,
            "TAAA TTTG GAAA TTCT CCTC AGCA TAAC ATTT TAAT GATG ATTC CTAG GATA"
            " CAAG TGAT GTGC TGAA AGTA TCAA TGTG TATC AGAA AACC AACA TCTC TGCT TAGG TCTC T",
        ),
    ],
)
def test_new_power_config(locus, sequence, bracketed, conc, lus, sec, tert, flank_5, flank_3):
    marker = STRMarkerObject(locus, sequence, "straitrazor", kit="powerseq")
    assert marker.convert == bracketed
    assert str(marker.canonical) == conc
    assert marker.designation == (lus, sec, tert)
    assert marker.flank_5p == flank_5
    assert marker.flank_3p == flank_3


@pytest.mark.parametrize(
    "locus, sequence, cust_seq, bracketed",
    [
        (
            "VWA",
            "GGATAGATGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGATAGATCAAT"
            "CCAAGTCACATACTGATTATTCTTATCATCCACTAGGGCTTTCACATCTCAGCCAAGTCAACTTGGATCCTCTAGACCTGTTTCTTCT"
            "TCTGGAA",
            "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGATAGA",
            "[TAGA]11 [CAGA]4 TAGA",
        ),
        (
            "TPOX",
            "CACTGGCCTGTGGGTCCCCCCATAGATCGTAAGCCCAGGAGGAAGGGCTGTGTTTCAGGGCTGTGATCACTAGCAC"
            "CCAGAACCGTCGACTGGCACAGAACAGGCACTTAGGGAACCCTCACTGAATGAATGAATGAATGAATGAATGAATGAATGTTTG"
            "GGCAAATAAACGCT",
            "AATGAATGAATGAATGAATGAATGAATGAATG",
            "[AATG]8",
        ),
        (
            "TH01",
            "CTCCATGGTGAATGAATGAATGAATGAATGAATGAATGAGGGAAATAAGGGAGGAACAGGCCAATGGGAATCACCC"
            "CAGAGCCCAGATACCCTTTGAATTTTGCCCCCTATTTGCCCAGGACCCCCCACCATGAGCTGCTGCTAGAGCCTGGGAAGGGCC"
            "TTGGGGCTGCCTCCCCAAGCAGGCAGGCTGGTTGGGGTGC",
            "ATGGTGAATGAATGAATGAATGAATGAATGAATGAGGGA",
            "ATGG TG [AATG]7 AGGG A",
        ),
        (
            "TH01",
            "CTCCATGGTGAATGAATGAATGAATGAATGAATGAATGAATGAATGAGGGAAATAAGGGAGGAACAGGCCAATGGG"
            "AATCACCCCAGAGCCCAGATACCCTTTGAATTTTGCCCCTATTTGCCCAGGACCCCCCACCATGAGCTGCTGCTAGAGCCTGGG"
            "AAGGGCCTTGGGGCTGCCTCCCCAAGCAGGCAGGCTGGTTGGGGTGC",
            "ATGGTGAATGAATGAATGAATGAATGAATGAAT" "GAATGAATGAGGG",
            "ATGG TG [AATG]8 AAT GAGG G",
        ),
        (
            "PENTA E",
            "TAATGATTACATAACATACATGTGTGTAAAGTGCTTAGTATCATGATTGATACATGGAAAGAATTCTCTTATT"
            "TGGGTTATTAATTGAGAAAACTCCTTACAATTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTGAGAC",
            "TCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTT",
            "[TCTTT]7",
        ),
        (
            "PENTA D",
            "GAGCCATGATCACACCACTACACTCCAGCCTAGGTGACAGAGCAAGACACCATCTCAAGAAAGAAAAAAAAGA"
            "AAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAAACGAAGGGGAAAAAAAGAGAATCATAAACATAAATG"
            "TAAAATTTCTCAAAAAAATCGTTA",
            "AAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGAAAAGA",
            "[AAAGA]9",
        ),
        (
            "FGA",
            "GTCTGAAATCGAAAATATGGTTATTGAAGTAGCTGCTGAGTGATTTGTCTGTAATTGCCAGCAAAAAAGAAAGGAAG"
            "AAAGGAAGGAAGGAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAGAAAAAAGAAAGAAAGAAA",
            "GGAAGGAAGGAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAGAAAAAAGAAAGAAAGAAA",
            "[GGAA]2 GGAG [AAAG]12 AGAA AAAA [GAAA]3",
        ),
        (
            "D8S1179",
            "TTTCATGTGTACATTCGTATCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC"
            "TATTCCCCACAGTGAAAATAATCTACAGGATAGGTAAATAAATTAAGGCATATTCACGCAATGGGATACGATACAGTGATGAAA"
            "ATGAACTAATTATAGCTACGTGAAAC",
            "TCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATC" "TA",
            "TCTA TCTG [TCTA]12",
        ),
        (
            "D7S820",
            "AGAATTGCACCAAATATTGGTAATTAAATGTTTACTATAGACTATTTAGTGAGATAAAAAAAAACTATCAATCT"
            "GTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCGTTAGTTCGTTCTAAACTATGACAAGTGTTCTA"
            "TCATACCCTTTATATATATTAACCTTAAAATAACTC",
            "AGATAAAAAAAAACTATCAATCTGTCTATCTATCTATCTATCTA" "TCTATCTATCTATCTATCTATCTATCTATCGTTA",
            "AGAT [AAAA]2 AC TATC AATC TGTC [TATC]12 GTTA",
        ),
        (
            "D5S818",
            "AACATTTGTATCTTTATCTGTATCCTTATTTATACCTCTATCTATCTATCTATCTATCTATCTATCTATCTATC"
            "TATCTATCTTCAAAATATTACGTAAGGATACCAAAGAGGAAAATCACCCTTGTCACATACTTGCTATTAAAATATACTTTTATT"
            "AGTACA",
            "ATACCTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTTCAA",
            "ATAC CTCT [ATCT]11 TCAA",
        ),
        (
            "D3S1358",
            "TGCCCACTTCTGCCCAGGGATCTATTTTTCTGTGGTGTGTATTCCCTGTGCCTTTGGGGGCATCTCTTATACT"
            "CATGAAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAT"
            "CTATGAGACAGGGTCTTGCTCTGTC",
            "CATGTATCTATCTGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTAT" "CTATCTATCTATGAGA",
            "CATG TA TCTA [TCTG]2 [TCTA]12 TGAG A",
        ),
        (
            "D2S441",
            "TGCACCCAACATTCTAACAAAAGGCTGTAACAAGGGCTACAGGAATCATGAGCCAGGAACTGTGGCTCATCTAT"
            "GAAAACTTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTTATCTATCTATATCA",
            "AACTTCTATCTA" "TCTATCTATCTATCTATCTATCTATCTATCTATCTATTTATCTATCTATATCA",
            "AACT [TCTA]11 TTTA [TCTA]2 TATC A",
        ),
        (
            "D2S1338",
            "CTAGCATGGTACCTGCAGGTGGCCCATAATAATGAGTTATTCAGTAAGTTAAAGGATTGCAGGAGGGAAGGAA"
            "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGCAGGCAGGCAGGCAGGCAGGCAAGGCCAAGCCATTTCTGTTTCCAA"
            "ATCCACTGGCTCCCTCCCACAGCT",
            "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGCAGGCAGGCA" "GGCAGGCAGGCA",
            "[GGAA]11 [GGCA]6",
        ),
        (
            "D22S1045",
            "CCTTCTTATAGCTGCTATGGGGGCTAGATTTTCCCCGATGATAGTAGTCTCATTATTATTATTATTATTATT"
            "ATTATTATTATTATTACTATTATTGTTATAAAAATATTG",
            "TCTCATTATTATTATTATTATTATTATTATTATTATTATTA" "CTATTATTGTTA",
            "TCT C [ATT]12 ACT [ATT]2 GTT A",
        ),
        (
            "D21S11",
            "TGAATTGCCTTCTATCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATC"
            "TATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCGTCTATCTATCCAGT"
            "CTATCTACCTCCTATTAGTCT",
            "GCCTTCTATCTATCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTATCTATCT"
            "ATATCTATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCGTCT",
            "GCCT [TCTA]6 [TCTG]5 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCA TA [TCTA]11 TCGT CT",
        ),
        (
            "D1S1656",
            "GAAATAGAATCACTAGGGAACCAAATATATATACATACAATTAAACACACACACACCTATCTATCTATCTATC"
            "TATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATCTATCTACATCATACAGTTGACCCTTGA",
            "CCTATCTATC" "TATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCTATCTATCTA",
            "CCTA [TCTA]11 TCA [TCTA]4",
        ),
        (
            "D19S433",
            "AAGTTCTTTAGCAGTGATTTCTGATATTTTGGTGCACCCATTACCCGAATAAAAATCTTCTCTCTTTCTTCCT"
            "CTCTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTACCTTCTTTCCTTCAACAGAATCTTATTC"
            "TGTTGCCCAGGCTGGAGTGCA",
            "ATCTTCTCTCTTTCTTCCTCTCTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTTCCTT"
            "CCTTCCTTCCTTCCTACCTTCTTTCCTTCAACA",
            "ATC TTCT CTCT TTCT TCCT CTCT [CCTT]12 CCTA CCTT CTTT CCTT CAAC A",
        ),
        (
            "D18S51",
            "AGGCTGCAGTGAGCCATGTTCATGCCACTGCACTTCACTCTGAGTGACAAATTGAGACCTTGTCTCAGAAAGAA"
            "AGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGAGAGAGGAAAGAAAGAGAAAAAGAAAAGAAATAGTAGCAA"
            "CTGTTATTGTA",
            "TCTCAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAAGAGAGAGGAAAGAAA",
            "TCTC [AGAA]12 AAAG AGAG AGGA AAGA AA",
        ),
        (
            "D16S539",
            "GTGCACAAATCTAAATGCAGAAAAGCACTGAAAGAAGAATCCAGAAAACCACAGTTCCCATTTTTATATGGGA"
            "GCAAACAAAGGCAGATCCCAAGCTCTTCCTCTTCCCTAGATCAATACAGACAGACAGACAGGTGGATAGATAGATAGATAGATA"
            "GATAGATAGATAGATATCAT",
            "GATAGATAGATAGATAGATAGATAGATAGATAGATA",
            "[GATA]9",
        ),
        (
            "D13S317",
            "TTCTTTAGTGGGCATCCGTGACTCTCTGGACTCTGACCCATCTAACGCCTATCTGTATTTACAAATACATTAT"
            "CTATCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTCTTTTTGGGCTGCCTATGGCT"
            "CAACCCAAGTTGAAGGAGGAGATTT",
            "ACATTATCTATCTATCTATCTATCTATCTATCTATCTATCAATCAATCATCTATC" "TATCTTTCTGTCTGTCTTTTT",
            "ACAT [TATC]9 [AATC]2 [ATCT]3 TTCT [GTCT]2 TTTT",
        ),
        (
            "D12S391",
            "CAATGGATGCATAGGTAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACA"
            "GACAGACAGACAGACAGATGAGAGGG",
            "TGCATAGGTAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATA"
            "GATAGACAGACAGACAGACAGACAGACAGATGAGA",
            "TGCA TAGG T [AGAT]12 [AGAC]6 AGAT GAGA",
        ),
        (
            "D10S1248",
            "CCCCAGGACCAATCTGGTCACAAACATATTAATGAATTGAACAAATGAGTGAGTGGAAGGAAGGAAGGAAGG"
            "AAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA",
            "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA",
            "[GGAA]14",
        ),
        (
            "CSF1PO",
            "CTAAGTACTTCCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAATCTATCTATCTTCTATCTA"
            "TGAAGGCAGTTACTGTTAATATCTTCATTTTACAGGTAGGAAAACTGAGACACAGGGTGGTTAGCAACCTGCTAGTCCTTGGCA"
            "GACTCAG",
            "CTTCCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTAATCT",
            "CTTC CT [ATCT]10 A ATCT",
        ),
        (
            "Y-GATA-H4",
            "TCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTACCTACCTATCTATCTATAG"
            "ATCTATCTATCTATCTTAAATTTGGAAATTCTCCTCAGCATAACATTTTAATGATGATTCCTAGGATACAAGTGATGTGCTGAA"
            "AGTATCAATGTGTATCAGAAAACCAACATCTCTGCTTAGGTCTCT",
            "TCTATCTATCTATCTATCTATCTATCTATCTATCT" "ATCTATCTATCTA",
            "[TCTA]12",
        ),
        (
            "DYS643",
            "ACCTCATGCTCTGTGATTTTTGCAGGTGTTCACTGCAAGCCATGCCTGGTTAAACTACTGTGCCTTTTCTTTTC"
            "TTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTCTTTTTAAAACTTTTTACTTCAGTAGAATTTTGGGGGG",
            "GTGCCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTCTTTT",
            "GTGC [CTTTT]10 CTTT CTTTT",
        ),
        (
            "DYS635",
            "CCCAAATATCCATCAATCAATGAATGGATAAAGAAAATGTGATAGATAGATAGATAGATAGATAGATAGATAGA"
            "TAGATAGATACATACATAGATAGATACATACATAGATAGATACATACATAGATAGATAGATAGAGATT",
            "ATGTGATAGATA"
            "GATAGATAGATAGATAGATAGATAGATAGATAGATACATACATAGATAGATACATACATAGATAGATACATACATAGATAGATA"
            "GATAGAGATT",
            "ATGT GA [TAGA]10 [TACA]2 [TAGA]2 [TACA]2 [TAGA]2 [TACA]2 [TAGA]4 GATT",
        ),
        (
            "DYS576",
            "AAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAAA"
            "GCCAAGACAAATACGCTTATTACTCCCATCTCCTCCTTCATCTCCAGGAAATGAGAC",
            "AAAGAAAGAAAGAAAGAAAGAAA" "GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAG",
            "[AAAG]17",
        ),
        (
            "DYS570",
            "TAAAATGAATGATGACTAGGTAGAAATCCTGGCTGTGTCCTCCAAGTTCCTTTTCTTTCTTTCTTTCTTTCTTT"
            "CTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTTT",
            "TTTCTTTCTTTCTTTCTTTCTTTCTTTCTT" "TCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTC",
            "[TTTC]17",
        ),
        (
            "DYS549",
            "GTAAAGAACTATAAAAAGATTAATACAACAAAAATTTGGTAATCTGAAATAATAAGGTAGACATAGCAATTAGG"
            "TAGGTAAAGAGGAAGATGATAGATGATTAGAAAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATA"
            "GATAGAAAAAATC",
            "AGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAAA",
            "AGAT [GATA]13 GAAA AA",
        ),
        (
            "DYS533",
            "CTAATATTTATCTATATCATTCTAATTATGTCTCTTCTAACTATATAACTATGTATTATCTATCAATCTTCTAC"
            "CTATCATCTTTCTAGCTAGCTATCATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCATCTATCAT"
            "CTTCTATTGTTTGGTTGAGTTAAGAACTGATCATGAATAAATACATTTCATTGGT",
            "TATCATCTATCTATCTATCTATCTA" "TCTATCTATCTATCTATCTATCTATCTATCATCT",
            "TATC ATC [TATC]12 ATCT",
        ),
        (
            "DYS481",
            "TAAAAGGAATGTGGCTAACGCTGTTCAGCATGCTGCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTT"
            "CTTCTTCTTCTTCTTCTTCTTCTTCTTTTTTGAGTCTTG",
            "CTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCTTCT" "TCTTCTTCTTCTTCTTCTTCTTCTT",
            "[CTT]22",
        ),
        (
            "DYS458",
            "GAAAGAAAGAAAAGGAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA"
            "GAAAGAAAGAAAGGAGGGTGGGCGTGGTGGCTCATGCTTGTAATGCCAGAACTTTGGGAGGCCGAGGTGG",
            "GAAAGAAAGA"
            "AAAGGAAGGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA",
            "[GAAA]3 AG GAAG [GAAA]17",
        ),
        (
            "DYS456",
            "AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATTCCATTAGTTCT"
            "GTCCCTCTAGAGAACCCTAATACATCAGTTTAAGAA",
            "AGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT" "AGATAGATAGATAGATATTCC",
            "[AGAT]15 ATTC C",
        ),
        (
            "DYS448",
            "AGACATGGATAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAG"
            "AGATATAGAGATAGAGAGATAGAGATAGAGATAGATAGATAGAGAAAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAG"
            "AGATAGAGATAGAGAGGTAAAGATAGA",
            "AGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGATAGAGA"
            "TAGAGATAGAGATATAGAGATAGAGAGATAGAGATAGAGATAGATAGATAGAGAAAGAGATAGAGATAGAGATAGAGATAGAGA"
            "TAGAGATAGAGATAGAGAT",
            "[AGAGAT]11 [ATAGAG]2 [AGATAG]3 ATAGAT AGAGAA [AGAGAT]8",
        ),
        (
            "DYS439",
            "AAATAGAAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAGTATAAGTAA"
            "AGAGATGATGG",
            "TAGAAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAAAGT",
            "TAGA A [GATA]13 GAAA GT",
        ),
        (
            "DYS438",
            "CAGTATATTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTATTTGA"
            "AATGGAGTTTCACTCTTGTTGCCCAGG",
            "TATATTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTTCTTTT" "CTTTTCTTTTCTATTT",
            "TATA [TTTTC]12 TATTT",
        ),
        (
            "DYS437",
            "GCCCATCCGGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGTCTGTCTATCTATCTATCTATCAT"
            "CTATCATCTGTGAATGATGTCTATCTACTTATCTATGAATGATATTTATCTGTGGTTATCTATCTATCTATA",
            "CCGGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGTCTGTCTATCTATCTATCTATCATCT",
            "CCGG [TCTA]9 [TCTG]2 [TCTA]4 TCAT CT",
        ),
        (
            "DYS393",
            "CGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATATGTATGTCTTTTCTATGAGAC" "ATA",
            "CGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAT",
            "CGAT [AGAT]12",
        ),
        (
            "DYS392",
            "TAAATAATAATAATAATAATAATAATAATAATAATAATAATAAATAAATGGTGATACAAGAAAAAAATTTGTTT"
            "TCCTTCTTGGCTTTTAAATAACAAACACTTGAAATCAAATTAG",
            "ATAATAATAATAATAATAATAATAATAATAATAATAA" "TA",
            "[ATA]13",
        ),
        (
            "DYS391",
            "TGTCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTGCCTATCTGCCTGCCTACCTA"
            "TCCCTCTATGGCAATTGCTTGCAACCAGGGAGATTTTA",
            "TCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATC" "TATCTATCTGCCTATC",
            "TCTG [TCTA]11 TCTG CCTA TC",
        ),
        (
            "DYS390",
            "AACAAGGAAAGATAGATAGATGATAGATAGATAGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAG"
            "ATAGATAGACAGACAGACAGACAGACAGACAGACAGACAGATAGATAGAATATATTATGGGGTACCAAAATGCAGGGCCCAAAA"
            "ATGTGTAAAATATATGTGT",
            "TAGATAGATAGATAGACAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGAC"
            "AGACAGACAGACAGACAGACAGACAGACAGA",
            "[TAGA]4 CAGA [TAGA]10 [CAGA]8",
        ),
        (
            "DYS389II",
            "TCATAGATAGATGATGGACTGCTAGATAAATAGATAGATTGATAGAGGGAGGGATAGATAGATAGATAGATA"
            "GATAGATAGATAGATAGATAGATAGACAGACAGACAGATACATAGATAATACAGATGAGAGTTGGATACAGAAGTAGGTATAAT"
            "GATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACAGACA",
            "TAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGATACATAGATAATACAGATGAGAGTTGGA"
            "TACAGAAGTAGGTATAATGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACAGA",
            "[TAGA]11 [CAGA]3 TACA TAGA TAAT ACAG ATGA GAGT TGGA TACA GAAG TAGG TATA ATGA [TAGA]11"
            " [CAGA]5",
        ),
        (
            "DYS385A-B",
            "TTTTTCTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCCCTTCCTTCCTTCCT"
            "TCCTTCCTTTCTTTCTCTTTCCTCTTTCTCTTTCTTCTCTTTCTTTCTTTTTCTCTTTTTCTCTTTCTTTCTTTTTTACTTTCT"
            "TTCTCCTTCCTTCCTTCCTTTCTGAATTTCATTTCTTTTCTTT",
            "TTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCTTTCT" "TTCTTTCTTTC",
            "[TTTC]12",
        ),
        (
            "DYS19",
            "TCTGGGTTAAGGAGAGTGTCACTATATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTACCTAT"
            "CTATCTATCTAAAACACTATATATATATAACACTATATATATAATACTATATATATATTA",
            "TCTATCTATCTATCTATCTA" "TCTATCTATCTATCTATCTATCTACCTATCTATCTATCTA",
            "[TCTA]11 CCTA [TCTA]3",
        ),
    ],
)
def test_custom_ranges(locus, sequence, cust_seq, bracketed):
    marker = STRMarkerObject(locus, sequence, "straitrazor", custom=True, kit="powerseq")
    assert marker.custom_sequence == cust_seq
    assert marker.custom_brack == bracketed
