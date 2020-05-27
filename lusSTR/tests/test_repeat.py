#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import lusSTR
from lusSTR.repeat import collapse_tandem_repeat, collapse_all_repeats
from lusSTR.repeat import split_by_n, get_blocks
from lusSTR.repeat import collapse_repeats_by_length, sequence_to_bracketed_form
import pytest


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
    final_output = collapse_all_repeats(sequence, repeat_list)
    assert final_output == output


def test_split_by_n():
    sequence = 'AGGTAGGTAGGTCGAACGAATTGG'
    blocks = list(split_by_n(sequence, n=4))
    assert blocks == [
        'AGGT', 'AGGT', 'AGGT', 'CGAA', 'CGAA', 'TTGG'
    ]


def test_sequence_to_bracketed_form():
    sequence = (
        'TCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATC'
        'TATCTATCTATCTATCTATCTATCTATCTATCTA'
    )
    repeats = ['TCTA', 'TCTG']
    final_output = sequence_to_bracketed_form(sequence, 6, repeats)
    assert final_output == '[TCTA]3 [TCTG]6 [TCTA]3 TA [TCTA]3 TCA [TCTA]2 TCCATA [TCTA]11'


def test_collapse_repeats_by_length():
    sequence = 'TCTATCTATCTATCTATCTATCTATCTATATATCTATCTATCTATCTA'
    assert collapse_repeats_by_length(sequence, 4) == '[TCTA]7 TATA [TCTA]4'
