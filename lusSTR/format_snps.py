#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import glob
import json
import lusSTR
import pandas as pd
import os
from pkg_resources import resource_filename


def get_snp_metadata_file():
    return resource_filename('lusSTR', 'snp_data.json')


with open(get_snp_metadata_file(), 'r') as fh:
    snp_marker_data = json.load(fh)


snp_type_dict = {
    'a': 'Ancestry',
    'i': 'Identity',
    'p': 'Phenotype',
    'p/a': 'Phenotype;Ancestry'
}


def uas_load(indir, type='i'):
    '''Format SNP data from UAS output files for use with `lusSTR annotate_snps`'''
    snp_final_output = pd.DataFrame()
    files = glob.glob(os.path.join(indir, '*.xlsx'))
    for filename in sorted(files):
        if 'Phenotype' in filename or 'Sample Details' in filename:
            snps = uas_format(filename, type)
            snp_final_output = snp_final_output.append(snps)
        else:
            continue
    return snp_final_output


def parse_snp_table_from_sheet(infile, sheet, snp_type_arg):
    table = pd.read_excel(io=infile, sheet_name=sheet)
    offset = table[table.iloc[:, 0] == "Coverage Information"].index.tolist()[0]
    data = table.iloc[offset + 2:]
    data.columns = table.iloc[offset + 1]
    data = data[['Locus', 'Reads', 'Allele Name']]
    final_df = pd.DataFrame()
    for type in snp_type_arg:
        filtered_dict = {k: v for k, v in snp_marker_data.items() if type in v['Type']}
        filtered_data = data[data['Locus'].isin(filtered_dict)].reset_index()
        if type == 'i':
            filtered_data['Type'] = 'Identity'
        elif type == 'p':
            filtered_data['Type'] = 'Phenotype'
        elif type == 'a':
            filtered_data['Type'] = 'Ancestry'
        else:
            filtered_data['Type'] = 'Phenotype;Ancestry'
        final_df = final_df.append(filtered_data)
    final_df['SampleID'] = table.iloc[1, 1]
    final_df['Project'] = table.iloc[2, 1]
    final_df['Analysis'] = table.iloc[3, 1]
    return final_df


def uas_format(infile, snp_type_arg):
    if 'all' in snp_type_arg:
        type = ['i', 'a', 'p']
    else:
        type = snp_type_arg
    if 'Phenotype' in infile and 'i' in type:
        snp_data = parse_snp_table_from_sheet(infile, 'SNP Data', type)
    elif 'Sample Details' in infile and ('a' or 'p' in type):
        snp_data = parse_snp_table_from_sheet(infile, 'iSNPs', type)
    else:
        snp_data = ''
    return snp_data


def strait_razor_concat(indir, snp_type_arg):
    '''Format a directory of STRait Razor output files for use with `lusSTR annotate_snps`.'''
    snps_final = pd.DataFrame()
    analysisID = os.path.basename(indir.rstrip(os.sep))
    files = glob.glob(os.path.join(indir, '*.txt'))
    for filename in sorted(files):
        name = filename.replace('.txt', '').split(os.sep)[-1]
        table = pd.read_csv(
            filename, sep='\t', header=None,
            names=['Locus_allele', 'Length', 'Sequence', 'Total_Reads', 'Other_Reads']
        )
        try:
            table[['SNP', 'Bases_off']] = table.Locus_allele.str.split(":", expand=True)
        except ValueError:
            print(
                f'Error found with {filename}. Will bypass and continue. Please check file'
                f' and rerun the command, if necessary.'
            )
            continue
        snps = []
        if 'all' in snp_type_arg:
            print('yes')
            snps_only = pd.DataFrame(
                table[table['SNP'].str.contains('rs')]
            ).reset_index(drop=True)
            for j, row in snps_only.iterrows():
                snpid = snps_only.iloc[j, 5]
                try:
                    metadata = snp_marker_data[snpid]
                except KeyError:
                    continue
                snp_type = metadata['Type']
                seq = snps_only.iloc[j, 2]
                snp_call = seq[metadata['Coord']]
                row_tmp = [
                    snpid, snps_only.iloc[j, 3], snp_call, snp_type_dict[snp_type], name,
                    analysisID, analysisID, seq, snps_only.iloc[j, 6]
                ]
                snps.append(row_tmp)
        else:
            for type in snp_type_arg:
                for j, row in table.iterrows():
                    snpid = table.iloc[j, 5]
                    try:
                        metadata = snp_marker_data[snpid]
                    except KeyError:
                        continue
                    snp_type = metadata['Type']
                    seq = snps_only.iloc[j, 3]
                    snp_call = seq[metadata['Coord']]
                    if type in snp_type:
                        row_tmp = [
                            snpid, snps_only.iloc[j, 3], snp_call, snp_type_dict[snp_type], name,
                            analysisID, analysisID, seq, snps_only.iloc[j, 6]
                        ]
                        snps.append(row_tmp)
        snps_tmp = pd.DataFrame(snps)
        snps_final = snps_final.append(snps_tmp)
        snps_final.columns = [
            'SNP', 'Reads', 'Allele', 'Type', 'SampleID', 'Project', 'Analysis', 'Sequence',
            'Bases_off'
        ]
    return snps_final


def main(args):
    if args.uas:
        results = uas_load(args.input, args.type)
    else:
        results = strait_razor_concat(args.input, args.type)
    if args.out is None:
        args.out = sys.stdout
    results.to_csv(args.out, index=False)
