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


def uas_load(inpath, sexloci=False):
    '''Format a UAS Sample Details Report (.xlsx) for use with `lusSTR annotate`.

    The `inpath` argument can refer to a report file or a directory of report files. Any files
    without the `.xlsx` file extension are ignored. The `sexloci` argument determines whether X and
    Y chromosome STRs are included along with autosomal STRs.
    '''
    if os.path.isdir(inpath):
        files = os.listdir(inpath)
        auto_strs = pd.DataFrame()
        sex_strs = pd.DataFrame() if sexloci is True else None
        for filename in sorted(files):
            if not filename.endswith('.xlsx'):
                continue
            filepath = os.path.join(inpath, filename)
            autodata, sexdata = uas_format(filepath, sexloci)
            auto_strs = auto_strs.append(autodata)
            if sexloci is True:
                sex_strs = sex_strs.append(sexdata)
    else:
        auto_strs, sex_strs = uas_format(inpath, sexloci)
    return auto_strs, sex_strs


def parse_str_table_from_sheet(infile, sheet, exclude=None):
    table = pd.read_excel(io=infile, sheet_name=sheet)
    offset = table[table.iloc[:, 0] == "Coverage Information"].index.tolist()[0]
    data = table.iloc[offset + 2:]
    data.columns = table.iloc[offset + 1]
    if exclude is not None:
        data = data[~data.Locus.isin(exclude)]
    data = data[['Locus', 'Reads', 'Repeat Sequence']]
    data['SampleID'] = table.iloc[1, 1]
    data['Project'] = table.iloc[2, 1]
    data['Analysis'] = table.iloc[3, 1]
    return data


def uas_format(infile, sexloci=False):
    auto_strs = parse_str_table_from_sheet(infile, sheet=0, exclude=['Amelogenin'])
    sex_strs = None
    if sexloci is True:
        y_strs = parse_str_table_from_sheet(infile, 'Y STRs')
        x_strs = parse_str_table_from_sheet(infile, 'X STRs')
        sex_strs = pd.concat([y_strs, x_strs], ignore_index=True)
    return auto_strs, sex_strs


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
        results, sex_results = uas_load(args.input, args.sex)
    else:
        sex_results, results = strait_razor_concat(args.input, args.sex)
    if args.out is None:
        args.out = sys.stdout
    results.to_csv(args.out, index=False)
    if args.sex:
        name = os.path.splitext(args.out)[0]
        sex_results.to_csv(f'{name}_sexloci.csv', index=False)
