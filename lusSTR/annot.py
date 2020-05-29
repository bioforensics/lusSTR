#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import csv
import json
import os
import pandas as pd
import re

import lusSTR
from lusSTR.repeat import collapse_all_repeats, collapse_repeats_by_length
from lusSTR.repeat import sequence_to_bracketed_form, split_by_n
from lusSTR.repeat import reverse_complement, reverse_complement_bracketed


def split_sequence_into_two_strings(sequence, repeat_for_split):
    '''
    Function to split a sequence into two separate strings at a specified repeat unit.
    '''
    last = 0
    prev = 0
    for m in re.finditer(repeat_for_split, sequence):
        if m.start() == prev or m.start() == last or prev == 0:
            prev = m.end()
        else:
            last = m.end()
    first_string = sequence[:prev]
    second_string = sequence[prev:]
    return first_string, second_string


def main(args):
    data = pd.read_csv(args.input)
    list_of_lists = []
    flanks_list = []
    for i, row in data.iterrows():
        locus = data.iloc[i, 0]
        reads = data.iloc[i, 1]
        sequence = data.iloc[i, 2]
        sampleid = data.iloc[i, 3]
        try:
            project = data.iloc[i, 4]
            analysis = data.iloc[i, 5]
        except IndexError:
            project = 'NA'
            analysis = 'NA'

        marker = lusSTR.marker.STRMarkerObject(locus, sequence, uas=args.uas, kit=args.kit)
        summary = [sampleid, project, analysis, locus] + marker.summary + [reads]
        list_of_lists.append(summary)

        if not args.uas and args.kit == 'forenseq':
            flank_summary = [
                sampleid, project, analysis, locus, reads, marker.canonical, marker.sequence,
                marker.flank_5p, marker.annotation, marker.flank_3p
            ]
            flanks_list.append(flank_summary)

    columns = [
        'SampleID', 'Project', 'Analysis', 'Locus', 'UAS_Output_Sequence',
        'Forward_Strand_Sequence', 'Traditional_STR_Allele', 'Forward_Strand_Bracketed_form',
        'UAS_Output_Bracketed_Form', 'LUS', 'LUS_Plus', 'Reads'
    ]
    final_output = pd.DataFrame(list_of_lists, columns=columns)
    name = os.path.splitext(args.out)[0]
    if not args.uas:
        flanks_columns = [
            'SampleID', 'Project', 'Analysis', 'Locus', 'Reads', 'Length_Allele', 'Full_Sequence',
            '5_Flank_Anno', 'UAS_Region_Anno', '3_Flank_Anno'
        ]
        final_flank_output = pd.DataFrame(flanks_list, columns=flanks_columns)
        final_flank_output.to_csv(f'{name}_flanks_anno.txt', sep='\t', index=False)
        if args.combine:
            final_output = final_output.groupby(columns[:-1], as_index=False)['Reads'].sum()
            final_output.to_csv(args.out, sep='\t', index=False)
        else:
            final_output.to_csv(f'{name}_no_combined_reads.txt', sep='\t', index=False)
    else:
        final_output.to_csv(args.out, sep='\t', index=False)
