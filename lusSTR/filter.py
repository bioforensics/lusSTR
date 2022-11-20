#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import json
from re import L
import numpy as np
import os
import pandas as pd
import sys
from pathlib import Path
from pkg_resources import resource_filename

import lusSTR
from lusSTR.filter_settings import filters


strs = [
    'CSF1PO', 'D10S1248', 'D12S391', 'D13S317', 'D16S539', 'D17S1301', 'D18S51', 'D19S433',
    'D1S1656', 'D20S482', 'D21S11', 'D22S1045', 'D2S1338', 'D2S441', 'D3S1358', 'D4S2408',
    'D5S818', 'D6S1043', 'D7S820', 'D8S1179', 'D9S1122', 'FGA', 'PENTA D', 'PENTA E', 'TH01',
    'TPOX', 'VWA'
]


def get_filter_metadata_file():
    return resource_filename('lusSTR', 'filters.json')


with open(get_filter_metadata_file(), 'r') as fh:
    filter_marker_data = json.load(fh)


def process_strs(dict_loc, allele_des):
    final_df = pd.DataFrame()
    for key, value in dict_loc.items():
        data = dict_loc[key].reset_index(drop=True)
        if allele_des == 'ru':
            data_combine = data.groupby(
                ['SampleID', 'Locus', 'RU_Allele'], as_index=False
            )['Reads'].sum()
            data_order = data_combine.sort_values(by=['RU_Allele'], ascending=False)
        total_reads = data_order['Reads'].sum()
        metadata = filter_marker_data[key[1]]
        data_order = data_order.reindex(columns=[
            *data_order.columns.tolist(),
            'allele_type', 'stuttering_allele1', 'stuttering_allele2', 'allele1_ref_reads',
            'allele2_ref_reads', 'perc_noise', 'perc_stutter'
        ], fill_value=None)
        filtered_df = filters(data_order, metadata, total_reads, allele_des)
        final_df = final_df.append(filtered_df)
    return final_df


def EFM_output(df, outfile, separate=False):
    infile = df[df.allele_type != 'noise']
    infile_sort = infile.sort_values(by=['SampleID', 'Locus', 'RU_Allele'], ascending=True)
    infile_sort['merged'] = (
        infile_sort['RU_Allele'].astype(str)+', '+infile_sort['Reads'].astype(str)
    )
    alleles_typed_col = infile_sort.groupby(
        ['SampleID', 'Locus']
    ).agg({'merged': ','.join}).reset_index()
    alleles_typed_final = alleles_typed_col['merged'].str.split(',', expand=True)
    new_df = pd.DataFrame()
    al_num = 0
    for name, values in alleles_typed_final.iteritems():
        if int(name) % 2:
            new_df[name] = values
        else:
            new_df.insert(al_num, name, values)
            al_num += 1
    col_len = int(len(new_df.columns)/2)
    for k in range(col_len):
        all_num = k+1
        read_num = k+col_len
        new_df.rename(columns={new_df.columns[k]: f'Allele{all_num}'}, inplace=True)
        new_df.rename(columns={new_df.columns[read_num]: f'Height{all_num}'}, inplace=True)
    sampl_snp_ids = alleles_typed_col.loc[:, ['SampleID', 'Locus']]
    final_df2 = pd.concat([sampl_snp_ids, new_df], axis=1).reset_index(drop=True)
    final_df2 = final_df2.rename({'SampleID': 'SampleName', 'Locus': 'Marker'}, axis=1)
    ids_list = final_df2['SampleName'].unique()
    df_complete = pd.DataFrame()
    for id in ids_list:
        df_sub = final_df2[final_df2['SampleName'] == id]
        for strloc in strs:
            if strloc in df_sub['Marker'].values:
                continue
            else:
                new_row = pd.DataFrame({'SampleName': [id], 'Marker': [strloc]})
                df_sub = df_sub.append(new_row)
        df_order = df_sub.sort_values(by=['Marker'])
        if separate:
            Path('Separated_EFM_Files').mkdir(exist_ok=True)
            df_order = df_order.dropna(axis=1, how='all', inplace=False)
            df_order.to_csv(f'Separated_EFM_Files/{id}.csv', index=False)
        else:
            df_complete = df_complete.append(df_order)
            df_complete.to_csv(outfile, index=False)


def STRmix_output(df):
    data_combine = df.groupby(['SampleID', 'Locus', 'RU_Allele'], as_index=False)['Reads'].sum()
    dict_loc = {k: v for k, v in data_combine.groupby(['SampleID', 'Locus'])}
    final_df = pd.DataFrame()
    for key, value in dict_loc.items():
        data = dict_loc[key].reset_index(drop=True)
        metadata = filter_marker_data[key[1]]
        slope = metadata['Slope']
        intercept = metadata['Intercept']
        data['Size'] = data['RU_Allele']*intercept + slope
        final_df = final_df.append(data)
    final_df = final_df.rename({'RU_Allele': 'CE Allele'})
    id_list = final_df['SampleID'].unique()
    Path('STRmix_Files').mkdir(exist_ok=True)
    for id in id_list:
        df_sub = final_df[final_df['SampleID'] == id]
        df_sub.iloc[:, 1:].to_csv(f'STRmix_Files/{id}.csv', index=False)


def main(args):
    full_df = pd.read_csv(args.input, sep='\t')
    output_type = args.output.lower()
    allele_des = args.allele.lower()
    if args.out is None:
        args.out = sys.stdout
    if output_type != 'efm' and output_type != 'strmix':
        raise ValueError('Incorrect output type specified. Please use EFM or STRmix only!')
    if args.nofilters:
        if output_type == 'strmix':
            STRmix_output(full_df)
        else:
            full_df['allele_type'] = 'real_allele'
            EFM_output(full_df, args.out, args.separate)
    else:
        dict_loc = {k: v for k, v in full_df.groupby(['SampleID', 'Locus'])}
        final_df = process_strs(dict_loc, args.allele)
        if args.info:
            if args.out != sys.stdout:
                name = args.out.replace('.csv', '').split(os.sep)[-1]
                final_df.to_csv(f'{name}_sequence_info.csv', index=False)
            else:
                raise ValueError('No outfile provided. Please specify --out to create info file.')
        if output_type == 'efm':
            EFM_output(final_df, args.out, args.separate)
        else:
            STRmix_output(final_df)
