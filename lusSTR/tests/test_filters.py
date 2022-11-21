#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
import lusSTR
from lusSTR.filter_settings import get_filter_metadata_file
from lusSTR.tests import data_file
import pytest


with open(get_filter_metadata_file(), 'r') as fh:
    filter_marker_data = json.load(fh)


@pytest.mark.parametrize('filter, locus, total_reads, allele_reads, final_reads, pass_filt', [
    ('Analytical', 'VWA', 318, 19, 318, False),
    ('Analytical', 'VWA', 318, 240, 318, True),
    ('Detection', 'VWA', 600, 9, 591, False),
    ('Detection', 'VWA', 600, 123, 600, True)
])
def test_thresholds(filter, locus, total_reads, allele_reads, final_reads, pass_filt):
    metadata = filter_marker_data[locus]
    test_total_reads, test_passfilt = lusSTR.filter_settings.thresholds(
        filter, metadata, total_reads, allele_reads
    )
    assert test_total_reads == final_reads
    assert test_passfilt == pass_filt
