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


def strait_razor_concat(input_dir):
    '''
    Function to prepare STRait Razor output for use in the 'annotate' lusSTR command.

    STRait Razor outputs individual files per sample. The function formats each file
    appropriately ('Locus', 'Total Reads', 'Sequence', 'SampleID') and then concatenates
    all samples into one large file.
    '''
    loci_list = [
        'CSF1PO', 'D10S1248', 'D12S391', 'D13S317', 'D16S539', 'D17S1301', 'D18S51', 'D19S433',
        'D1S1656', 'D20S482', 'D21S11', 'D22S1045', 'D2S1338', 'D2S441', 'D3S1358', 'D4S2408',
        'D5S818', 'D6S1043', 'D7S820', 'D8S1179', 'D9S1122', 'FGA', 'PentaD', 'PentaE', 'TH01',
        'TPOX', 'vWA'
        ]
    myfiles = os.listdir(input_dir)
    straitrazorcomp = pd.DataFrame()
    for filename in myfiles:
        name = re.sub("_STRaitRazor.txt", "", filename)
        print(filename)
        file = pd.read_table(input_dir + filename, sep="\t", header=None)
        file.columns = ['Locus_allele', 'Length', 'Sequence', 'Forward_Reads', 'Reverse_Reads']
        file[['Locus', 'Allele']] = file.Locus_allele.str.split(":", expand=True)
        filtered_file = file[file['Locus'].isin(loci_list)]
        filtered_file['Total_Reads'] = (
            filtered_file['Forward_Reads'] + filtered_file['Reverse_Reads']
        )
        filtered_file['SampleID'] = name
        final_file = filtered_file.loc[:, ['Locus', 'Total_Reads', 'Sequence', 'SampleID']]
        straitrazorcomp = straitrazorcomp.append(final_file)
    straitrazorcomp.columns = ['Locus', 'Total_Reads', 'Sequence', 'SampleID']
    return straitrazorcomp


def main(args):
    '''
    Script to convert either the UAS Sample Details Report (.xlsx format using the --uas flag)
    or STRait Razor output to a more user-friendly format. Also removes the Amelogenin locus
    and extract relevant information (e.g. Sample ID, Project ID and Analysis ID).
    '''
    if args.uas:
        file = pd.read_excel(io=args.input, sheet_name=0)
        well_index = file[
            file["Sample Autosomal STR Report"] == "Coverage Information"].index.tolist()
        results_newdf = file[(well_index[0] + 2):]
        results_newdf.columns = file.iloc[(well_index[0] + 1)]
        results_filt = results_newdf[results_newdf.Locus != "Amelogenin"]
        results_final = results_filt[['Locus', 'Reads', 'Repeat Sequence']]
        results_final['SampleID'] = file.iloc[1, 1]
        results_final['Project'] = file.iloc[2, 1]
        results_final['Analysis'] = file.iloc[3, 1]
    else:
        results_final = strait_razor_concat(args.input)
        print(args.input)
        analysisID = re.sub("_FASTQ/", "", args.input)
        results_final['Project'] = "NA"
        results_final['Analysis'] = re.sub("_FASTQ/", "", args.input)

    output_file = sys.stdout
    if args.out is not None:
        results_final.to_csv(args.out, index=False)
