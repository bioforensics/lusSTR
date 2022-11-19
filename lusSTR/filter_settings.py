#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
import lusSTR
import numpy as np
import pandas as pd
from pkg_resources import resource_filename


def get_filter_metadata_file():
    return resource_filename('lusSTR', 'filters.json')


with open(get_filter_metadata_file(), 'r') as fh:
    filter_marker_data = json.load(fh)


strs = [
    'CSF1PO', 'D10S1248', 'D12S391', 'D13S317', 'D16S539', 'D17S1301', 'D18S51', 'D19S433',
    'D1S1656', 'D20S482', 'D21S11', 'D22S1045', 'D2S1338', 'D2S441', 'D3S1358', 'D4S2408',
    'D5S818', 'D6S1043', 'D7S820', 'D8S1179', 'D9S1122', 'FGA', 'PENTA D', 'PENTA E', 'TH01',
    'TPOX', 'VWA'
]


def filters(data_order, metadata, total_reads, allele_des):
    if len(data_order) == 1:
        if thresholds('Detection', metadata, total_reads, data_order['Reads'][0])[1] is False:
            data_order = pd.DataFrame()
        elif thresholds('Analytical', metadata, total_reads, data_order['Reads'][0])[1] is False:
            data_order[['allele_type', 'perc_noise']] = ['noise', 1]
        elif thresholds('Analytical', metadata, total_reads, data_order['Reads'][0])[1] is True:
            data_order['allele_type'] = 'real_allele'
    else:
        for i in range(len(data_order)):  # check for alleles below detection threshold
            al_reads = data_order.loc[i, 'Reads']
            if thresholds('Detection', metadata, total_reads, al_reads)[1] is False:
                data_order = data_order.drop(data_order.index[i])
                total_reads = total_reads - al_reads
        data_order = data_order.reset_index(drop=True)
        for i in range(len(data_order)):  # check for alleles below AT threshold
            ref_allele_reads = data_order.loc[i, 'Reads']
            if thresholds('Analytical', metadata, total_reads, ref_allele_reads)[1] is False:
                data_order.loc[i, ['allele_type', 'perc_noise']] = [
                    'noise', round(ref_allele_reads/total_reads, 3)
                ]
            else:
                data_order.loc[i, 'allele_type'] = 'real_allele'
        for i in range(len(data_order)):  # check for stutter alleles
            if data_order.loc[i, 'allele_type'] != 'real_allele':
                continue
            else:
                ref_allele_reads = data_order.loc[i, 'Reads']
                for j in range(len(data_order)):
                    if j == i:
                        continue
                    init_type_all = data_order.loc[j, 'allele_type']
                    if init_type_all == 'noise':
                        continue
                    al_reads = data_order.loc[j, 'Reads']
                    ref_allele = float(data_order.loc[i, 'RU_Allele'])
                    question_allele = float(data_order.loc[j, 'RU_Allele'])
                    data_order.loc[j, ['allele_type', 'perc_stutter']] = allele_type_ru(
                        ref_allele, question_allele, init_type_all, metadata, al_reads,
                        ref_allele_reads, data_order.loc[j, 'allele1_ref_reads'], data_order
                    )
                    if 'stutter' in data_order.loc[j, 'allele_type']:
                        if (
                            '/' in data_order.loc[j, 'allele_type'] and
                            pd.isnull(data_order.loc[j, 'stuttering_allele2'])
                        ):
                            data_order.loc[j, ['stuttering_allele2', 'allele2_ref_reads']] = [
                                ref_allele, ref_allele_reads
                            ]
                        elif pd.isnull(data_order.loc[j, 'stuttering_allele1']):
                            data_order.loc[j, ['stuttering_allele1', 'allele1_ref_reads']] = [
                                ref_allele, ref_allele_reads
                            ]
            for j in range(len(data_order)):
                type_all = data_order.loc[j, 'allele_type']
                if 'stutter' in type_all:
                    if '/' not in type_all:
                        if pd.isnull(data_order.loc[j, 'perc_stutter']):
                            data_order.loc[j, 'perc_stutter'] = round(
                                data_order.loc[j, 'Reads'] /
                                data_order.loc[j, 'allele1_ref_reads'], 3
                            )
                    else:
                        data_order.loc[j, 'perc_stutter'] = ''
                    data_order.loc[j, 'perc_noise'] = ''
                elif 'noise' in data_order.loc[j, 'allele_type']:
                    data_order.loc[j, 'perc_noise'] = round(
                            data_order.loc[j, 'Reads']/total_reads, 3
                        )
    return data_order


def thresholds(filter, metadata, total_reads, al_reads):
    use = metadata[f'{filter}ThresholdUse']
    count = metadata[f'{filter}ThresholdStaticCount']
    perc = metadata[f'{filter}ThresholdDynamicPercent']
    thresh_perc = round(perc * total_reads, 1)
    if (
        use.lower() == 'dynamic' and
        total_reads < metadata['MinimumNumberReadsForDynamicThresholds']
    ):
        use = 'static'
    if use.lower() == 'both':
        if thresh_perc >= count:
            thresh = thresh_perc
        else:
            thresh = count
    elif use.lower() == 'static':
        thresh = count
    elif use.lower() == 'dynamic':
        thresh = thresh_perc
    else:
        print(
            f'Specified {filter} Threshold Use Not Acceptable!'
            f'Check filters.json file and re-run.'
        )
    if al_reads < thresh:
        if filter == 'dynamic':
            total_reads = total_reads - al_reads
        return total_reads, False
    else:
        return total_reads, True


def minus1_stutter(
    all_type, stutter_thresh, forward_thresh, stutter_thresh_reads,
    ref_reads, al1_ref_reads, al_reads
):
    stut_perc = None
    if all_type == '+1_stutter':
        all_thresh = (
            stutter_thresh_reads +
            forward_stut_thresh(forward_thresh, stutter_thresh, al1_ref_reads)
        )
        all_type = output_allele_call(
            al_reads, all_thresh, '-1_stutter/+1_stutter'
        )
    elif all_type == '-2_stutter':
        all_thresh = stutter_thresh_reads + np.ceil(stutter_thresh * al1_ref_reads)
        all_type = output_allele_call(
            al_reads, all_thresh, '-1_stutter/-2_stutter'
        )
    elif al_reads <= stutter_thresh_reads:
        all_type = '-1_stutter'
        stut_perc = round(al_reads/ref_reads, 3)
    else:
        all_type = 'real_allele'
    return all_type, stut_perc


def minus2_stutter(
    all_type, stutter_thresh, forward_thresh, stutter_thresh_reads,
    al1_ref_reads, al_reads
):
    stut_perc = None
    if all_type == '-1_stutter':
        all_thresh = (
            stutter_thresh_reads +
            np.ceil(stutter_thresh * al1_ref_reads)
        )
        all_type = output_allele_call(
            al_reads, all_thresh, '-1_stutter/-2_stutter'
        )
    elif all_type == '+1_stutter':
        all_thresh = (
            stutter_thresh_reads +
            forward_stut_thresh(forward_thresh, stutter_thresh, al1_ref_reads)
        )
        all_type = output_allele_call(
            al_reads, all_thresh, '+1_stutter/-2_stutter'
        )
    elif al_reads <= stutter_thresh_reads:
        all_type = '-2_stutter'
        stut_perc = al_reads / al1_ref_reads
    else:
        all_type = 'real_allele'
    return all_type, stut_perc


def plus1_stutter(all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads, al_reads):
    stut_perc = None
    if all_type == '-1_stutter':
        all_thresh = (
            np.ceil(stutter_thresh * al1_ref_reads) +
            forward_stut_thresh(forward_thresh, stutter_thresh, ref_reads)
        )
        all_type = output_allele_call(
            al_reads, all_thresh, '-1_stutter/+1_stutter'
        )
    elif all_type == '-2_stutter':
        all_thresh = (
            np.ceil(stutter_thresh * al1_ref_reads) +
            forward_stut_thresh(forward_thresh, stutter_thresh, ref_reads)
        )
        all_type = output_allele_call(
            al_reads, all_thresh, '+1_stutter/-2_stutter'
        )
    elif al_reads <= forward_stut_thresh(forward_thresh, stutter_thresh, ref_reads):
        all_type = '+1_stutter'
        stut_perc = round(al_reads/ref_reads, 3)
    else:
        all_type = 'real_allele'
    return all_type, stut_perc


def allele_type_ru(ref, ru, all_type, metadata, al_reads, ref_reads, al1_ref_reads, data):
    stutter_thresh = metadata['StutterThresholdDynamicPercent']
    forward_thresh = metadata['StutterForwardThresholdDynamicPercent']
    stutter_thresh_reads = np.ceil(stutter_thresh * ref_reads)
    stut_perc = None
    if ref - ru == 1:  # -1 stutter
        all_type, stut_perc = minus1_stutter(
            all_type, stutter_thresh, forward_thresh, stutter_thresh_reads, ref_reads,
            al1_ref_reads, al_reads
        )
    elif ref - ru == 2:  # -2 stutter
        if check_2stutter(data, 'ru', ru)[0] is True:
            ref_reads = check_2stutter(data, 'ru', ru)[1]
            stutter_thresh_reads = stutter_thresh * ref_reads
            all_type, stut_perc = minus2_stutter(
                all_type, stutter_thresh_reads, forward_thresh, stutter_thresh_reads, ref_reads,
                al1_ref_reads, al_reads
            )
    elif ref - ru == -1:  # +1 stutter
        all_type, stut_perc = plus1_stutter(
            all_type, stutter_thresh_reads, forward_thresh, ref_reads,
            al1_ref_reads, al_reads
        )
    elif pd.isnull(all_type):
        all_type = 'noise'
    return all_type, stut_perc


def output_allele_call(al_reads, all_thresh, orig_type):
    if al_reads <= all_thresh:
        all_type = orig_type
    else:
        all_type = 'real_allele'
    return all_type


def forward_stut_thresh(perc, perc_stut, reads):
    if perc == 0:
        forward_thresh = np.ceil(perc_stut**2 * reads)
    else:
        forward_thresh = np.ceil(perc * reads)
    return forward_thresh


def check_2stutter(data, allele_des, allele):
    is_true, reads = False, None
    if '-1_stutter' in data.loc[:, 'allele_type'].values:
        if allele_des == 'ru':
            for k, row in data.iterrows():
                ru_test = data.loc[k, 'RU_Allele']
                if ru_test - allele == 1 and data.loc[k, 'allele_type'] == '-1_stutter':
                    is_true, reads = True, data.loc[k, 'Reads']
                    break
    return is_true, reads
