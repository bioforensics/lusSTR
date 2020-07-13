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


def uas_load(input_file):
    '''Format a UAS Sample Details Report (.xlsx) for use with `lusSTR annotate`.'''
    data = pd.read_excel(io=input_file, sheet_name=0)
    offset = data[data.iloc[:, 0] == "Coverage Information"].index.tolist()[0]
    results = data.iloc[offset + 2:]
    results.columns = data.iloc[offset + 1]
    results = results[results.Locus != 'Amelogenin'][['Locus', 'Reads', 'Repeat Sequence']]
    results['SampleID'] = data.iloc[1, 1]
    results['Project'] = data.iloc[2, 1]
    results['Analysis'] = data.iloc[3, 1]
    return results


def strait_razor_concat(input_dir):
    '''Format a directory of STRait Razor output files for use with `lusSTR annotate`.'''
    locus_list = [
        'CSF1PO', 'D10S1248', 'D12S391', 'D13S317', 'D16S539', 'D17S1301', 'D18S51', 'D19S433',
        'D1S1656', 'D20S482', 'D21S11', 'D22S1045', 'D2S1338', 'D2S441', 'D3S1358', 'D4S2408',
        'D5S818', 'D6S1043', 'D7S820', 'D8S1179', 'D9S1122', 'FGA', 'PentaD', 'PENTAD',
        'Penta D', 'PentaE', 'PENTAE', 'Penta E', 'TH01', 'TPOX', 'vWA', 'VWA'
        ]
    myfiles = os.listdir(input_dir)
    alldata = pd.DataFrame()
    for filename in sorted(myfiles):
        name = re.sub('_STRaitRazor.txt', '', filename)
        filepath = os.path.join(input_dir, filename)
        data = pd.read_csv(
            filepath, sep='\t', header=None,
            names=['Locus_allele', 'Length', 'Sequence', 'Forward_Reads', 'Reverse_Reads']
        )
        data[['Locus', 'Allele']] = data.Locus_allele.str.split(":", expand=True)
        data = data[data.Locus.isin(locus_list)]
        data['Total_Reads'] = data['Forward_Reads'] + data['Reverse_Reads']
        data['SampleID'] = name
        data = data[['Locus', 'Total_Reads', 'Sequence', 'SampleID']]
        alldata = alldata.append(data)
    analysisID = input_dir.rstrip(os.sep)
    analysisID_final = os.path.basename(analysisID)
    alldata['Project'] = analysisID_final
    alldata['Analysis'] = analysisID_final
    return alldata


def main(args):
    if args.uas:
        results = uas_load(args.input)
    else:
        results = strait_razor_concat(args.input)
    if args.out is None:
        args.out = sys.stdout
    results.to_csv(args.out, index=False)
