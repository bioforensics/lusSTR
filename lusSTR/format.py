#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import lusSTR
import argparse
import os
import pandas as pd
import re
import sys


def uas_load(input_file, sex=False):
    '''Format a UAS Sample Details Report (.xlsx) for use with `lusSTR annotate`.'''
    if os.path.isdir(input_file):
        myfiles = os.listdir(input_file)
        autosomal_data = pd.DataFrame()
        sex_data = pd.DataFrame()
        for filename in sorted(myfiles):
            filepath = os.path.join(input_file, filename)
            sex_final, results = uas_format(filepath, sex)
            autosomal_data = autosomal_data.append(results)
            if sex is True:
                sex_data = sex_data.append(sex_final)
    else:
        sex_data, autosomal_data = uas_format(filepath, sex)
    return sex_data, autosomal_data


def uas_format(input_file, sex=False):
    data = pd.read_excel(io=input_file, sheet_name=0)
    sampleID = data.iloc[1, 1]
    projectID = data.iloc[2, 1]
    analysisID = data.iloc[3, 1]
    offset = data[data.iloc[:, 0] == "Coverage Information"].index.tolist()[0]
    results = data.iloc[offset + 2:]
    results.columns = data.iloc[offset + 1]
    results = results[results.Locus != 'Amelogenin'][['Locus', 'Reads', 'Repeat Sequence']]
    results['SampleID'] = sampleID
    results['Project'] = projectID
    results['Analysis'] = analysisID
    if sex is True:
        y_data = pd.read_excel(io=input_file, sheet_name='Y STRs')
        y_offset = y_data[y_data.iloc[:, 0] == "Coverage Information"].index.tolist()[0]
        y_strs = y_data.iloc[y_offset + 2:]
        y_strs.columns = y_data.iloc[y_offset + 1]
        y_strs = y_strs[['Locus', 'Reads', 'Repeat Sequence']]
        x_data = pd.read_excel(io=input_file, sheet_name='X STRs')
        x_offset = x_data[x_data.iloc[:, 0] == "Coverage Information"].index.tolist()[0]
        x_strs = x_data.iloc[x_offset + 2:]
        x_strs.columns = x_data.iloc[x_offset + 1]
        x_strs = x_strs[['Locus', 'Reads', 'Repeat Sequence']]
        sex_final = pd.concat([y_strs, x_strs], ignore_index=True)
        sex_final['SampleID'] = sampleID
        sex_final['Project'] = projectID
        sex_final['Analysis'] = analysisID
    else:
        sex_final = ''
    return sex_final, results


def strait_razor_concat(input_dir, sex=False):
    '''Format a directory of STRait Razor output files for use with `lusSTR annotate`.'''
    locus_list = [
        'CSF1PO', 'D10S1248', 'D12S391', 'D13S317', 'D16S539', 'D17S1301', 'D18S51', 'D19S433',
        'D1S1656', 'D20S482', 'D21S11', 'D22S1045', 'D2S1338', 'D2S441', 'D3S1358', 'D4S2408',
        'D5S818', 'D6S1043', 'D7S820', 'D8S1179', 'D9S1122', 'FGA', 'PentaD', 'PENTAD',
        'Penta D', 'PentaE', 'PENTAE', 'Penta E', 'TH01', 'TPOX', 'vWA', 'VWA'
        ]
    sex_locus_list = [
        'Y-GATA-H4', 'DYS643', 'DYS635', 'DYS612', 'DYS576', 'DYS570', 'DYS549', 'DYS533',
        'DYS522', 'DYS505', 'DYS481', 'DYS460', 'DYS448', 'DYS439', 'DYS438', 'DYS437', 'DYS392',
        'DYS391', 'DYS390', 'DYS389I', 'DYS389II', 'DYS385a-b', 'DYS19', 'DYF387S1', 'DYS393',
        'DYS458', 'DYS456', 'HPRTB', 'DXS8378', 'DXS7423', 'DXS7132', 'DXS10135', 'DXS10074',
        'DXS10103', 'DYS385'
    ]
    myfiles = os.listdir(input_dir)
    autosomal_data = pd.DataFrame()
    xydata = pd.DataFrame()
    for filename in sorted(myfiles):
        name = re.sub('_STRaitRazor.txt', '', filename)
        filepath = os.path.join(input_dir, filename)
        data = pd.read_csv(
            filepath, sep='\t', header=None,
            names=['Locus_allele', 'Length', 'Sequence', 'Forward_Reads', 'Reverse_Reads']
        )
        data[['Locus', 'Allele']] = data.Locus_allele.str.split(":", expand=True)
        data['Total_Reads'] = data['Forward_Reads'] + data['Reverse_Reads']
        data['SampleID'] = name
        data = data[['Locus', 'Total_Reads', 'Sequence', 'SampleID']]
        auto_only_data = data[data.Locus.isin(locus_list)]
        autosomal_data = autosomal_data.append(auto_only_data)
        if sex is True:
            sex_only_data = data[data.Locus.isin(sex_locus_list)]
            xydata = xydata.append(sex_only_data)
    analysisID = input_dir.rstrip(os.sep)
    analysisID_final = os.path.basename(analysisID)
    autosomal_data['Project'] = analysisID_final
    autosomal_data['Analysis'] = analysisID_final
    if sex is True:
        xydata['Project'] = analysisID_final
        xydata['Analysis'] = analysisID_final
    return xydata, autosomal_data


def main(args):
    if args.uas:
        sex_results, results = uas_load(args.input, args.sex)
    else:
        sex_results, results = strait_razor_concat(args.input, args.sex)
    if args.out is None:
        args.out = sys.stdout
    results.to_csv(args.out, index=False)
    if args.sex:
        name = os.path.splitext(args.out)[0]
        sex_results.to_csv(f'{name}_sexloci.csv', index=False)
