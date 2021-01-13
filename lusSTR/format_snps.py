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


def strait_razor_concat(indir, type):
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
        if type == 'all':
            snps_only = pd.DataFrame(table[table['SNP'].str.contains('rs')]).reset_index(drop=True)
            for j, row in snps_only.iterrows():
                snpid = snps_only.iloc[j, 5]
                try:
                    metadata = snp_marker_data[snpid]
                except KeyError:
                    continue
                snp_type = metadata['Type']
                row_tmp = [
                    snpid, snps_only.iloc[j, 2], snps_only.iloc[j, 3], snps_only.iloc[j, 6],
                    snp_type_dict[snp_type]
                ]
                snps.append(row_tmp)
        else:
            for type in ['i']:
                for j, row in table.iterrows():
                    snpid = table.iloc[j, 5]
                    try:
                        metadata = snp_marker_data[snpid]
                    except KeyError:
                        continue
                    snp_type = metadata['Type']
                    if type in snp_type:
                        row_tmp = [
                            snpid, table.iloc[j, 2], table.iloc[j, 3], table.iloc[j, 6],
                            snp_type_dict[snp_type]
                        ]
                        snps.append(row_tmp)
        snps_tmp = pd.DataFrame(snps)
        snps_tmp['SampleID'] = name
        snps_tmp['Project'] = analysisID
        snps_tmp['Analysis'] = analysisID
        snps_final = snps_final.append(snps_tmp)
        snps_final.columns = [
            'SNP', 'Sequence', 'Reads', 'Bases_off', 'Type', 'SampleID', 'Project', 'Analysis'
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
