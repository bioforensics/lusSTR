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


snps_within_loci = {
    'mh16-MC1RB': {
        'SNPs': ['rs1805005', 'rs1805006', 'rs2228479']
    },
    'mh16-MC1RC': {
        'SNPs': [
            'rs11547464',
            'rs1805007',
            'rs201326893_Y152OCH',
            'rs1110400',
            'rs1805008',
            'rs885479'
            ]
    }
}


def complement_base(base):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    comp_base = complement[base]
    return comp_base


def uas_load(indir, type='i'):
    '''Format SNP data from UAS output files'''
    snp_final_output = pd.DataFrame()
    files = glob.glob(os.path.join(indir, '[!~]*.xlsx'))
    for filename in sorted(files):
        if 'Phenotype' in filename or 'Sample Details' in filename:
            snps = uas_types(filename, type)
            if snps is not None:
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
    if snp_type_arg == 'all':
        final_df = data
    elif snp_type_arg == 'i':
        filtered_dict = {k: v for k, v in snp_marker_data.items() if 'i' in v['Type']}
        filtered_data = data[data['Locus'].isin(filtered_dict)].reset_index()
        final_df = final_df.append(filtered_data)
    else:
        filtered_dict = {k: v for k, v in snp_marker_data.items() if 'i' not in v['Type']}
        filtered_data = data[data['Locus'].isin(filtered_dict)].reset_index()
        final_df = final_df.append(filtered_data)
    final_df['SampleID'] = table.iloc[1, 1]
    final_df['Project'] = table.iloc[2, 1]
    final_df['Analysis'] = table.iloc[3, 1]
    return final_df


def uas_types(infile, snp_type_arg):
    if 'Sample Details' in infile and (snp_type_arg == 'all' or snp_type_arg == 'i'):
        snp_data = parse_snp_table_from_sheet(infile, 'iSNPs', snp_type_arg)
    elif 'Phenotype' in infile and (snp_type_arg == 'all' or snp_type_arg == 'p'):
        snp_data = parse_snp_table_from_sheet(infile, 'SNP Data', snp_type_arg)
    else:
        snp_data = None
    return snp_data


def uas_format(infile, snp_type_arg):
    data = uas_load(infile, snp_type_arg)
    data_filt = data.loc[data['Reads'] != 0].reset_index()
    data_df = []
    for j, row in data_filt.iterrows():
        snpid = data_filt.iloc[j, 2]
        metadata = snp_marker_data[snpid]
        type = metadata['Type']
        uas_allele = data_filt.iloc[j, 4]
        if metadata['ReverseCompNeeded'] == 'Yes':
            forward_strand_allele = complement_base(uas_allele)
        else:
            forward_strand_allele = uas_allele
        row_tmp = [
            data_filt.iloc[j, 5], data_filt.iloc[j, 6], data_filt.iloc[j, 7], snpid,
            data_filt.iloc[j, 3], forward_strand_allele, uas_allele, snp_type_dict[type]
        ]
        data_df.append(row_tmp)
    data_final = pd.DataFrame(data_df, columns=[
        'SampleID', 'Project', 'Analysis', 'SNP', 'Reads', 'Forward_Strand_Allele', 'UAS_Allele',
        'Type'
    ])
    return data_final


def compile_row_of_snp_data(infile, snp, table_loc, type, name, analysis):
    snp_df = []
    if 'mh16' in snp:
        locus_data = snps_within_loci[snp]
        for k in range(0, len(locus_data['SNPs'])):
            snp_id = locus_data['SNPs'][k]
            row_tmp = collect_snp_info(infile, snp_id, table_loc, type, name, analysis)
            if row_tmp is not None:
                snp_df.append(row_tmp)
    else:
        row_tmp = collect_snp_info(infile, snp, table_loc, type, name, analysis)
        if row_tmp is not None:
            snp_df.append(row_tmp)
    final_snp_df = pd.DataFrame(snp_df)
    return final_snp_df


def snp_call_exception(seq, expected_size, metadata, base):
    new_size = len(seq) + expected_size
    new_base_call = seq[new_size]
    if new_base_call in metadata['Alleles']:
        flag = (
            'Sequence length different than expected (check for indels); allele position adjusted'
        )
        return new_base_call, flag
    else:
        flag = (
            'Allele call does not match expected allele! Check for indels '
            '(does not match expected sequence length)'
        )
        return base, flag


def collect_snp_info(infile, snpid, j, type, name, analysis):
    if snpid == 'N29insA':
        snpid = 'rs312262906_N29insA'
    metadata = snp_marker_data[snpid]
    snp_type = metadata['Type']
    seq = infile.iloc[j, 2]
    expected_alleles = metadata['Alleles']
    snp_loc = metadata['Coord']
    snp_call = seq[snp_loc]
    if snpid == 'rs312262906_N29insA' and snp_call == 'A':
        snp_call = 'insA'
    if metadata['ReverseCompNeeded'] == 'Yes':
        snp_call_uas = complement_base(snp_call)
    else:
        snp_call_uas = snp_call
    if snpid == 'rs2402130':
        expected_size = len(seq) - 73
    else:
        expected_size = int(infile.iloc[j, 6])
    if snp_call not in expected_alleles and expected_size != '0':
        if snpid == 'rs1821380':
            snp_call, allele_flag = snp_call_exception(seq, expected_size, metadata, snp_call)
            snp_call_uas = complement_base(snp_call)
        else:
            allele_flag = (
                'Allele call does not match expected allele! Check for indels '
                '(does not match expected sequence length)'
            )
    elif snp_call not in expected_alleles:
        allele_flag = 'Allele call does not match expected allele!'
    elif expected_size != 0:
        allele_flag = 'Check for indels (does not match expected sequence length)'
    else:
        allele_flag = ''
    if (
        (type == 'p' and (snp_type == 'p' or snp_type == 'a' or snp_type == 'p/a')) or
        (type == 'i' and snp_type == 'i') or (type == 'all')
    ):
        row_tmp = [
            name, analysis, analysis, snpid, seq, infile.iloc[j, 7], snp_call, snp_call_uas,
            snp_type_dict[snp_type], allele_flag
        ]
    else:
        row_tmp = None
    return row_tmp


def strait_razor_concat(indir, snp_type_arg):
    '''Format a directory of STRait Razor output files.'''
    snps = pd.DataFrame()
    analysisID = os.path.basename(indir.rstrip(os.sep))
    files = glob.glob(os.path.join(indir, '[!~]*.txt'))
    for filename in sorted(files):
        name = filename.replace('.txt', '').split(os.sep)[-1]
        table = pd.read_csv(
            filename, sep='\t', header=None,
            names=['Locus_allele', 'Length', 'Sequence', 'Forward_Reads', 'Reverse_Reads']
        )
        try:
            table[['SNP', 'Bases_off']] = table.Locus_allele.str.split(":", expand=True)
        except ValueError:
            print(
                f'Error found with {filename}. Will bypass and continue. Please check file'
                f' and rerun the command, if necessary.'
            )
            continue
        table['Total_Reads'] = table['Forward_Reads'] + table['Reverse_Reads']
        snps_only = pd.DataFrame(
            table[table['SNP'].str.contains('rs|mh16|insA')]
        ).reset_index(drop=True)
        for j, row in snps_only.iterrows():
            snpid = snps_only.iloc[j, 5]
            try:
                row = compile_row_of_snp_data(snps_only, snpid, j, snp_type_arg, name, analysisID)
            except KeyError:
                continue
            if row is not None:
                snps = snps.append(row)
        snps.columns = [
            'SampleID', 'Project', 'Analysis', 'SNP', 'Sequence', 'Reads', 'Forward_Strand_Allele',
            'UAS_Allele', 'Type', 'Potential_Issues'
        ]
    return snps


def strait_razor_format(infile, snp_type_arg):
    '''
    This function formats STRait Razor input data for two separate reports. The Reads are summed
     for identical allele calls per SNP.
    '''
    results = strait_razor_concat(infile, snp_type_arg)
    results_combine = results.groupby(
        [
            'SNP', 'Forward_Strand_Allele', 'UAS_Allele', 'Type', 'SampleID', 'Project',
            'Analysis'
        ],
        as_index=False
    )['Reads'].sum()
    results_combine = results_combine[[
        'SampleID', 'Project', 'Analysis', 'SNP', 'Reads', 'Forward_Strand_Allele',
        'UAS_Allele', 'Type'
    ]]
    return results, results_combine


def main(args):
    if args.uas:
        results = uas_format(args.input, args.type)
        results.to_csv(args.out, index=False, sep='\t')
    else:
        results, results_combined = strait_razor_format(args.input, args.type)
        output_name = os.path.splitext(args.out)[0]
        results_combined.to_csv(args.out, index=False, sep='\t')
        results.to_csv(f'{output_name}_full_output.txt', index=False, sep='\t')
