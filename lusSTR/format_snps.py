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
    files = glob.glob(os.path.join(indir, '*.xlsx'))
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
    for type in snp_type_arg:
        filtered_dict = {k: v for k, v in snp_marker_data.items() if type in v['Type']}
        filtered_data = data[data['Locus'].isin(filtered_dict)].reset_index()
        final_df = final_df.append(filtered_data)
    final_df['SampleID'] = table.iloc[1, 1]
    final_df['Project'] = table.iloc[2, 1]
    final_df['Analysis'] = table.iloc[3, 1]
    return final_df


def uas_types(infile, snp_type_arg):
    if 'all' in snp_type_arg:
        type = ['i', 'a', 'p']
    else:
        type = snp_type_arg
    if 'Sample Details' in infile and 'i' in type:
        snp_data = parse_snp_table_from_sheet(infile, 'iSNPs', type)
    elif 'Phenotype' in infile and ('a' or 'p' in type):
        snp_data = parse_snp_table_from_sheet(infile, 'SNP Data', type)
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
            snpid, data_filt.iloc[j, 3], forward_strand_allele, uas_allele, snp_type_dict[type],
            data_filt.iloc[j, 5], data_filt.iloc[j, 6], data_filt.iloc[j, 7]
        ]
        if snpid == 'rs16891982' or snpid == 'rs12913832':
            match_on = [
                snpid, data_filt.iloc[j, 3], forward_strand_allele, uas_allele,
                data_filt.iloc[j, 5], data_filt.iloc[j, 6], data_filt.iloc[j, 7]
            ]
            if [s for s in data_df if all(xs in s for xs in match_on)]:
                continue
        data_df.append(row_tmp)
    data_final = pd.DataFrame(data_df, columns=[
        'SNP', 'Reads', 'Forward_Strand_Allele', 'UAS_Allele', 'Type', 'SampleID', 'Project',
        'Analysis'
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
    expected_size = infile.iloc[j, 6]
    if snp_call not in expected_alleles and expected_size != '0':
        allele_flag = (
            'Allele call does not match expected allele! Check for indels '
            '(does not match expected sequence length)'
        )
    elif snp_call not in expected_alleles:
        allele_flag = 'Allele call does not match expected allele!'
    elif expected_size != '0':
        allele_flag = 'Check for indels (does not match expected sequence length)'
    else:
        allele_flag = ''
    if type in snp_type or 'all' in type:
        row_tmp = [
            snpid, infile.iloc[j, 3], snp_call, snp_call_uas, snp_type_dict[snp_type], name,
            analysis, analysis, allele_flag, seq, infile.iloc[j, 6]
        ]
    else:
        row_tmp = None
    return row_tmp


def strait_razor_concat(indir, snp_type_arg):
    '''Format a directory of STRait Razor output files for use with `lusSTR annotate_snps`.'''
    snps = pd.DataFrame()
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
        if 'all' in snp_type_arg:
            snps_only = pd.DataFrame(
                table[table['SNP'].str.contains('rs|mh16|insA')]
            ).reset_index(drop=True)
            for j, row in snps_only.iterrows():
                snpid = snps_only.iloc[j, 5]
                try:
                    row = compile_row_of_snp_data(snps_only, snpid, j, 'all', name, analysisID)
                except KeyError:
                    continue
                if row is not None:
                    snps = snps.append(row)
        else:
            for type in snp_type_arg:
                for j, row in table.iterrows():
                    snpid = table.iloc[j, 5]
                    try:
                        row = compile_row_of_snp_data(table, snpid, j, type, name, analysisID)
                    except KeyError:
                        continue
                    if row is not None:
                        snps = snps.append(row)
        snps.columns = [
            'SNP', 'Reads', 'Forward_Strand_Allele', 'UAS_Allele', 'Type', 'SampleID', 'Project',
            'Analysis', 'Potential_Issues', 'Sequence', 'Bases_off'
        ]
    return snps


def strait_razor_format(infile, snp_type_arg):
    '''
    This function formats STRait Razor input data for two separate reports. The Reads are summed
     for identical allele calls per SNP. Because rs16891982 and rs12913832 are listed as both
     ancestry and phenotype SNPs, they are double counted when both 'a' and 'p' are specified as
     arguments; therefore, their read totals are divided by two.
    '''
    results = strait_razor_concat(infile, snp_type_arg)
    results_combine = results.groupby(
            [
                'SNP', 'Forward_Strand_Allele', 'UAS_Allele', 'Type', 'SampleID', 'Project',
                'Analysis'
            ],
            as_index=False
        )['Reads'].sum()
    if 'a' in snp_type_arg and 'p' in snp_type_arg:
        results_combine.loc[
            (results_combine['SNP'] == 'rs16891982') |
            (results_combine['SNP'] == 'rs12913832'), 'Reads'
        ] = results_combine['Reads']/2
        results_combine = results_combine.astype({'Reads': int})
    return results, results_combine


def main(args):
    if args.uas:
        results = uas_format(args.input, args.type)
    else:
        results, results_combined = strait_razor_format(args.input, args.type)
        output_name = os.path.splitext(args.out)[0]
        results_combined.to_csv(f'{output_name}_full_output.txt', index=False, sep='\t')
    if args.out is None:
        args.out = sys.stdout
    results.to_csv(args.out, index=False, sep='\t')
