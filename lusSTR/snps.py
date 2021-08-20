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
import openpyxl
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


def uas_format(infile, snp_type_arg):
    '''
    This function begins with the compiled data from all files within the specified directory.
    It removes any allele with Reads of 0; identifies whether the allele call needs to be reverse
    complemented to be reported on the forward strand; and checks that the called allele is one of
    two expected alleles for the SNP (and flags any SNP call which is unexpected).
    '''
    data_filt = uas_load(infile, snp_type_arg).reset_index(drop=True)
    data_df = []
    for j, row in data_filt.iterrows():
        snpid = data_filt.iloc[j, 0]
        metadata = snp_marker_data[snpid]
        type = metadata['Type']
        uas_allele = data_filt.iloc[j, 2]
        if metadata['ReverseCompNeeded'] == 'Yes':
            forward_strand_allele = complement_base(uas_allele)
        else:
            forward_strand_allele = uas_allele
        if data_filt.loc[j, 'Typed Allele?'] == 'No':
            flag = 'Contains untyped allele'
        elif forward_strand_allele in metadata['Alleles']:
            flag = ''
        else:
            flag = 'Allele call does not match expected allele!'
        row_tmp = [
            data_filt.loc[j, 'SampleID'], data_filt.loc[j, 'Project'],
            data_filt.loc[j, 'Analysis'], snpid, data_filt.loc[j, 'Reads'], forward_strand_allele,
            uas_allele, snp_type_dict[type], flag
        ]
        data_df.append(row_tmp)
    data_final = pd.DataFrame(data_df, columns=[
        'SampleID', 'Project', 'Analysis', 'SNP', 'Reads', 'Forward_Strand_Allele', 'UAS_Allele',
        'Type', 'Issues'
    ])
    data_final_sort = data_final.sort_values(
        by=['SampleID', 'Project', 'Analysis', 'SNP', 'Reads'], ascending=False
    )
    return data_final_sort


def uas_load(indir, type='i'):
    '''
    This function lists input .xlsx files within the specified directory and performs a check to
    ensure the correct file is processed (must contain either "Phenotype" or "Sample Details").
    This also compiles the SNP data for each file within the directory.
    '''
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


def uas_types(infile, snp_type_arg):
    '''
    This function determines which tab within the specified file is required to extract the SNP
    data from based on the name of the file.
    '''
    if 'Sample Details' in infile and (snp_type_arg == 'all' or snp_type_arg == 'i'):
        snp_data = parse_snp_table_from_sheet(infile, 'iSNPs', snp_type_arg)
    elif 'Phenotype' in infile and (snp_type_arg == 'all' or snp_type_arg == 'p'):
        snp_data = parse_snp_table_from_sheet(infile, 'SNP Data', snp_type_arg)
    else:
        snp_data = None
    return snp_data


def parse_snp_table_from_sheet(infile, sheet, snp_type_arg):
    '''
    This function formats the SNP data from the original file and filters the SNPs based on the
    indicated SNP type.
    '''
    file = openpyxl.load_workbook(infile)
    file_sheet = file[sheet]
    table = pd.DataFrame(file_sheet.values)
    offset = table[table.iloc[:, 0] == 'Coverage Information'].index.tolist()[0]
    data = table.iloc[offset + 2:]
    data.columns = table.iloc[offset + 1]
    data = data[['Locus', 'Reads', 'Allele Name', 'Typed Allele?']]
    final_df = pd.DataFrame()
    if snp_type_arg == 'all':
        final_df = data
    elif snp_type_arg == 'i':
        filtered_dict = {k: v for k, v in snp_marker_data.items() if 'i' in v['Type']}
        filtered_data = data[data['Locus'].isin(filtered_dict)].reset_index(drop=True)
        final_df = final_df.append(filtered_data)
    else:
        filtered_dict = {k: v for k, v in snp_marker_data.items() if 'i' not in v['Type']}
        filtered_data = data[data['Locus'].isin(filtered_dict)].reset_index(drop=True)
        final_df = final_df.append(filtered_data)
    final_df['SampleID'] = table.iloc[2, 1]
    final_df['Project'] = table.iloc[3, 1]
    final_df['Analysis'] = table.iloc[4, 1]
    return final_df


def strait_razor_format(infile, snp_type_arg):
    '''
    This function formats STRait Razor input data for two separate reports. The full output
    includes all reads, the SNP allele calls and any results flags. In the main report, the reads
    are summed for identical allele calls per SNP. This function also checks that the allele call
    is one of two expected alleles for the SNP (and flags the allele if not).
    '''
    results = strait_razor_concat(infile, snp_type_arg)
    results_sort = results.sort_values(
        by=['SampleID', 'Project', 'Analysis', 'SNP', 'Reads'], ascending=False
    )
    results_combine = results_sort.groupby(
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
    results_combine['Issues'] = ''
    for j, row in results_combine.iterrows():
        snpid = results_combine.iloc[j, 3]
        metadata = snp_marker_data[snpid]
        if results_combine.iloc[j, 5] not in metadata['Alleles']:
            results_combine.iloc[j, 8] = 'Allele call does not match expected allele!'
    results_combine_sort = results_combine.sort_values(
        by=['SampleID', 'Project', 'Analysis', 'SNP', 'Reads'], ascending=False
    )
    return results_sort, results_combine_sort


def strait_razor_concat(indir, snp_type_arg):
    '''
    This function reads in all .txt files within the specified directory. For each file, the
    forward and reverse reads are summed and each sequence is processed and compiled into one
    final dataframe.
    '''
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
        'SampleID', 'Project', 'Analysis', 'SNP', 'Sequence', 'Reads',
        'Forward_Strand_Allele', 'UAS_Allele', 'Type', 'Potential_Issues'
    ]
    return snps


def compile_row_of_snp_data(infile, snp, table_loc, type, name, analysis):
    '''
    This function is necessary to account for the two sets of SNPs reported from the same
    sequence amplicon. Sequences labeled as mh16-MC1RB and mh16-MC1RC contain 3 and 6 SNPs,
    respectively. This function reports out each SNP from the sequence amplicon as individual
    rows and calls another function to compile data on each SNP.
    '''
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
    '''
    This function compiles allele calls, reads, reverse complements allele call if necessary to
    match how the UAS reports the allele, and any flags associated with the allele call. The flags
    indicate potential issues, including an unexpected allele call (not one of two expected
    alleles for the SNP) or unexpected length of the sequence amplicon which could result in an
    incorrect allele call. This function also determines if the SNP should be included in the
    final table based on the specified SNP type from the CLI.
    '''
    if snpid == 'N29insA':
        snpid = 'rs312262906_N29insA'
    metadata = snp_marker_data[snpid]
    snp_type = metadata['Type']
    seq = infile.iloc[j, 2]
    expected_alleles = metadata['Alleles']
    snp_loc = metadata['Coord']
    if len(seq) > snp_loc:
        snp_call = seq[snp_loc]
        if snpid == 'rs312262906_N29insA' and snp_call == 'A':
            snp_call = 'insA'
        if metadata['ReverseCompNeeded'] == 'Yes':
            snp_call_uas = complement_base(snp_call)
        else:
            snp_call_uas = snp_call
        if snpid == 'rs2402130':
            differ_length = len(seq) - 73
        else:
            differ_length = int(infile.iloc[j, 6])
        if snp_call not in expected_alleles and differ_length != 0:
            if snpid == 'rs1821380':
                snp_call, allele_flag = snp_call_exception(seq, differ_length, metadata, snp_call)
                snp_call_uas = complement_base(snp_call)
            else:
                allele_flag = (
                    'Allele call does not match expected allele! Check for indels '
                    '(does not match expected sequence length)'
                )
        elif snp_call not in expected_alleles:
            allele_flag = 'Allele call does not match expected allele!'
        elif differ_length != 0:
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
    else:
        row_tmp = None
    return row_tmp


def snp_call_exception(seq, expected_size, metadata, base):
    '''
    This function accounts for insertions and deletions in sequences to identify the correct base
    coordinate for the SNP. If the identified allele is still not one of the expected alleles, the
    sequence will be flagged appropriately.
    '''
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


def indiv_files(table, input_dir, ext):
    output_dir = f'Separated_lusstr_Files/{input_dir}'
    os.makedirs(output_dir, exist_ok=True)
    for samp in table['SampleID'].unique():
        new_df = table[table['SampleID'] == samp]
        new_df.to_csv(f'{output_dir}/{samp}{ext}', sep='\t', index=False)


def main(args):
    output_name = os.path.splitext(args.out)[0]
    if args.uas:
        results = uas_format(args.input, args.type)
        if args.separate:
            indiv_files(results, output_name, '.txt')
        else:
            results.to_csv(args.out, index=False, sep='\t')
    else:
        results, results_combined = strait_razor_format(args.input, args.type)
        if args.separate:
            indiv_files(results_combined, output_name, '.txt')
        else:
            results_combined.to_csv(args.out, index=False, sep='\t')
        results.to_csv(f'{output_name}_full_output.txt', index=False, sep='\t')
