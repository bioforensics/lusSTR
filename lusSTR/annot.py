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
from pkg_resources import resource_filename


def get_str_metadata_file():
    return resource_filename('lusSTR', 'str_markers.json')


with open(get_str_metadata_file(), 'r') as fh:
    str_marker_data = json.load(fh)


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


def format_table(input, uas=False, kit='forenseq'):
    '''
    Function to format final output table and the flanking report (if necessary).
    '''
    data = pd.read_csv(input, keep_default_na=False)
    data.iloc[:, 3] = data.iloc[:, 3].astype(str)
    list_of_lists = []
    flanks_list = []
    for i, row in data.iterrows():
        locus = data.iloc[i, 0].upper()
        reads = data.iloc[i, 1]
        sequence = data.iloc[i, 2]
        sampleid = re.sub(" ", "_", data.iloc[i, 3])
        try:
            project = re.sub(" ", "_", data.iloc[i, 4])
            analysis = re.sub(" ", "_", data.iloc[i, 5])
        except IndexError:
            project = 'NA'
            analysis = 'NA'
        if locus == 'PENTAD' or locus == 'PENTA_D':
            locus = 'PENTA D'
        if locus == 'PENTAE' or locus == 'PENTA_E':
            locus = 'PENTA E'
        if locus == 'DYS385A-B' or locus == 'DYS385':
            locus = 'DYS385A-B'
        if locus == 'AMELOGENIN':
            continue
        metadata = str_marker_data[locus]
        if kit == 'forenseq':
            remove_5p = metadata['Foren_5']
            remove_3p = metadata['Foren_3']
        else:
            remove_5p = metadata['Power_5']
            remove_3p = metadata['Power_3']
        if len(sequence) <= (remove_5p + remove_3p) and not uas:
            flank_summary = [
                sampleid, project, analysis, locus, reads, 'NA', sequence, 'NA', 'NA', 'NA',
                'Partial sequence'
            ]
            flanks_list.append(flank_summary)
            continue
        elif 'N' in sequence:
            flank_summary = [
                sampleid, project, analysis, locus, reads, 'NA', sequence, 'NA', 'NA', 'NA',
                'Sequence contains Ns'
            ]
            flanks_list.append(flank_summary)
            continue

        marker = lusSTR.marker.STRMarkerObject(locus, sequence, uas=uas, kit=kit)
        summary = [sampleid, project, analysis, locus] + marker.summary + [reads]
        list_of_lists.append(summary)

        if not uas:
            flank_summary = [
                sampleid, project, analysis, locus, reads, marker.canonical, marker.sequence,
                marker.flank_5p, marker.annotation, marker.flank_3p, marker.indel_flag
            ]
            flanks_list.append(flank_summary)

    columns = [
        'SampleID', 'Project', 'Analysis', 'Locus', 'UAS_Output_Sequence',
        'Forward_Strand_Sequence', 'RU_Allele', 'Forward_Strand_Bracketed_Notation',
        'UAS_Output_Bracketed_Notation', 'LUS', 'LUS_Plus', 'Reads'
    ]
    if not list_of_lists:
        final_output = pd.DataFrame(list_of_lists, columns=columns)
    else:
        final_output = sort_table(pd.DataFrame(list_of_lists, columns=columns))
    if not uas:
        flanks_columns = [
            'SampleID', 'Project', 'Analysis', 'Locus', 'Reads', 'RU_Allele',
            'Full_Sequence', '5_Flank_Bracketed_Notation', 'UAS_Region_Bracketed_Notation',
            '3_Flank_Bracketed_Notation', 'Potential_Issues'
        ]
        if not flanks_list:
            final_flank_output = pd.DataFrame(flanks_list, columns=flanks_columns)
        else:
            final_flank_output = sort_table(pd.DataFrame(flanks_list, columns=flanks_columns))
    else:
        final_flank_output = ''
    return final_output, final_flank_output, columns


def combine_reads(table, columns):
    comb_table = table.groupby(columns[:-1], as_index=False)['Reads'].sum()
    sorted = sort_table(comb_table)
    return sorted


def sort_table(table):
    sorted_table = table.sort_values(
                by=['SampleID', 'Project', 'Analysis', 'Locus', 'Reads', 'RU_Allele'],
                ascending=False
            )
    return sorted_table


def main(args):
    output_name = os.path.splitext(args.out)[0]
    input_name = os.path.splitext(args.input)[0]
    autosomal_final_table, autosomal_flank_table, columns = format_table(
        args.input, args.uas, args.kit
    )
    if args.sex:
        sex_final_table, sex_flank_table, columns = format_table(
            f'{input_name}_sexloci.csv', args.uas, args.kit
        )
        if not args.uas:
            sex_flank_table.to_csv(
                f'{output_name}_sexloci_flanks_anno.txt', sep='\t', index=False
            )
            if args.combine:
                if not sex_final_table.empty:
                    sex_final_table = combine_reads(sex_final_table, columns)
                sex_final_table.to_csv(f'{output_name}_sexloci.txt', sep='\t', index=False)
            else:
                sex_final_table.to_csv(
                    f'{output_name}_sexloci_no_combined_reads.txt', index=False
                )
        else:
            sex_final_table.to_csv(f'{output_name}_sexloci.txt', sep='\t', index=False)
    if not args.uas:
        autosomal_flank_table.to_csv(f'{output_name}_flanks_anno.txt', sep='\t', index=False)
        if args.combine:
            if not autosomal_final_table.empty:
                autosomal_final_table = combine_reads(autosomal_final_table, columns)
            autosomal_final_table.to_csv(args.out, sep='\t', index=False)
        else:
            autosomal_final_table.to_csv(
                f'{output_name}_no_combined_reads.txt', sep='\t', index=False
            )
    else:
        autosomal_final_table.to_csv(args.out, sep='\t', index=False)
