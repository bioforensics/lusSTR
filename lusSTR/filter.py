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
import lusSTR
from lusSTR.filter_settings import allele_counts, allele_imbalance_check, filters
import numpy as np
import os
import pandas as pd
from pathlib import Path
from pkg_resources import resource_filename
import re
import sys


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


def process_strs(dict_loc, datatype):
    final_df = pd.DataFrame()
    flags_df = pd.DataFrame()
    for key, value in dict_loc.items():
        data = dict_loc[key].reset_index(drop=True)
        if datatype == 'ce':
            data_combine = data.groupby(
                ['SampleID', 'Locus', 'RU_Allele'], as_index=False
            )['Reads'].sum()
            data_order = data_combine.sort_values(
                by=['RU_Allele'], ascending=False
            ).reset_index(drop=True)
        else:
            data_combine = data[['SampleID', 'Locus', 'UAS_Output_Sequence', 'RU_Allele', 'Reads']]
            data_order = data_combine.sort_values(
                by=['RU_Allele'], ascending=True
            ).reset_index(drop=True)
        total_reads = data_order['Reads'].sum()
        locus = key[1]
        data_order = data_order.reindex(columns=[
            *data_order.columns.tolist(),
            'allele_type', 'stuttering_allele1', 'stuttering_allele2', 'allele1_ref_reads',
            'allele2_ref_reads', 'perc_noise', 'perc_stutter'
        ], fill_value=None)
        filtered_df = filters(data_order, locus, total_reads, datatype)
        final_df = final_df.append(filtered_df)
        flags_df = flags_df.append(allele_counts(filtered_df))
        flags_df = flags_df.append(allele_imbalance_check(filtered_df))
    final_df = final_df.astype({'RU_Allele': 'float64', 'Reads': 'int'})
    return final_df, flags_df


def EFM_output(df, outfile, profile, separate=False):
    if outfile is None:
        outfile = sys.stdout
    if profile == 'reference':
        infile = df[df.allele_type == 'real_allele']
    else:
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
            if profile == 'evidence':
                df_order.to_csv(f'Separated_EFM_Files/{id}.csv', index=False)
            else:
                df_order.iloc[:, :4].to_csv(f'Separated_EFM_Files/{id}.csv', index=False)
        else:
            df_complete = df_complete.append(df_order)
            if profile == 'evidence':
                df_complete.to_csv(outfile, index=False)
            else:
                for i in range(len(df_complete)):
                    if pd.isna(df_complete.loc[i, 'Allele2']):
                        df_complete.loc[i, 'Allele2'] = df_complete.loc[i, 'Allele1']
                name = outfile.replace('.csv', '')
                df_complete.iloc[:, :4].to_csv(f'{name}_reference.csv', index=False)


def STRmix_output(df, outdir, profile, datatype):
    df_filt = df[df['Locus'] != 'Amelogenin'].reset_index(drop=True)
    if profile == 'reference':
        infile = df_filt[df_filt.allele_type == 'real_allele']
    else:
        infile = df_filt[df_filt.allele_type != 'noise']
    if datatype == 'ce':
        data_combine = infile.groupby(
            ['SampleID', 'Locus', 'RU_Allele'], as_index=False
        )['Reads'].sum()
        dict_loc = {k: v for k, v in data_combine.groupby(['SampleID', 'Locus'])}
        final_df = pd.DataFrame()
        for key, value in dict_loc.items():
            data = dict_loc[key].reset_index(drop=True)
            metadata = filter_marker_data[key[1]]
            slope = metadata['Slope']
            intercept = metadata['Intercept']
            data['Size'] = data['RU_Allele']*intercept + slope
            final_df = final_df.append(data)
        final_df.rename(
            {'RU_Allele': 'Allele', 'Reads': 'Height'}, axis=1, inplace=True
        )
    else:
        final_df = infile[[
            'SampleID', 'Locus', 'RU_Allele', 'UAS_Output_Sequence', 'Reads'
        ]].copy()
        final_df.rename(
            {'RU_Allele': 'CE Allele', 'UAS_Output_Sequence': 'Allele Seq'}, axis=1, inplace=True
        )
    final_df.replace(
            {'Locus': {'VWA': 'vWA', 'PENTA D': 'PentaD', 'PENTA E': 'PentaE'}}, inplace=True
        )
    id_list = final_df['SampleID'].unique()
    if outdir is None:
        outdir = 'STRmix_Files'
    Path(outdir).mkdir(exist_ok=True)
    for id in id_list:
        df_sub = final_df[final_df['SampleID'] == id].reset_index(drop=True)
        if profile == 'evidence':
            df_sub.iloc[:, 1:].to_csv(f'{outdir}/{id}_{datatype}.csv', index=False)
        else:
            ref_df = reference_table(df_sub.iloc[:, 1:3])
            ref_df.to_csv(f'{outdir}/{id}_reference_{datatype}.csv', index=False)


def reference_table(data):
    new_rows = []
    for i, row in data.iterrows():
        locus = data.loc[i, 'Locus']
        try:
            next_col = data.loc[i+1, 'Locus']
        except KeyError:
            next_col = None
        try:
            prev_col = data.loc[i-1, 'Locus']
        except KeyError:
            prev_col = None
        if locus == next_col or locus == prev_col:
            continue
        else:
            new_rows.append(list(row))
    new_df = pd.DataFrame(new_rows, columns=['Locus', 'Allele'])
    concat_df = pd.concat([data, new_df]).reset_index(drop=True)
    sort_df = concat_df.sort_values(by=['Locus', 'Allele'])
    return sort_df


def main(args):
    full_df = pd.read_csv(args.input, sep='\t')
    if args.out is None:
        outpath = sys.stdout
    else:
        outpath = args.out
    if args.nofilters:
        full_df['allele_type'] = 'real_allele'
        if args.output == 'efm':
            EFM_output(full_df, outpath, args.profile, args.separate)
        else:
            STRmix_output(full_df, outpath, args.profile, args.data)
    else:
        dict_loc = {k: v for k, v in full_df.groupby(['SampleID', 'Locus'])}
        final_df, flags_df = process_strs(dict_loc, args.data)
        if args.output == 'efm':
            EFM_output(final_df, outpath, args.profile, args.separate)
        else:
            STRmix_output(final_df, outpath, args.profile, args.data)
        if args.info:
            if outpath != sys.stdout:
                if args.output == 'efm':
                    name = args.out.replace('.csv', '')
                    final_df.to_csv(f'{name}_sequence_info.csv', index=False)
                    if not flags_df.empty:
                        flags_df.to_csv(f'{name}_Flagged_Loci.csv', index=False)
                else:
                    if outpath == sys.stdout:
                        outpath = 'STRmix_Files'
                    final_df.to_csv(f'{outpath}/STRmix_Files_sequence_info.csv', index=False)
                    if not flags_df.empty:
                        flags_df.to_csv(f'{outpath}/Flagged_Loci.csv', index=False)
            else:
                raise ValueError('No outfile provided. Please specify --out to create info file.')
