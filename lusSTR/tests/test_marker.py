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
from lusSTR.tests import data_file
import re
from tempfile import NamedTemporaryFile


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


# @pytest.mark.parametrize('sequence, uas_seq, front, back', [
#     (
#         'CTATGCATCTATCTATCTATCTATCTATCTATCTATCTATCTAATGGTTA',
#         'ATCTATCTATCTATCTATCTATCTATCTATCTATCT',
#         6,
#         8,
#     ),
#     (
#         'TCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATTCCC',
#         'TCTATCTGTCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTA',
#         0,
#         5,
#     ),
# ])
# def test_full_seq_to_uas(sequence, uas_seq, front, back):
#     uas_sequence = lusSTR.annot.full_seq_to_uas(sequence, front, back)
#     assert uas_sequence == uas_seq
