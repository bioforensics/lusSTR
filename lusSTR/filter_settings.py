# -------------------------------------------------------------------------------------------------
# Copyright (c) 2020, DHS.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from collections import defaultdict
import json
import numpy as np
import pandas as pd
from pkg_resources import resource_filename
import re


def get_filter_metadata_file():
    return resource_filename("lusSTR", "filters.json")


with open(get_filter_metadata_file(), "r") as fh:
    filter_marker_data = json.load(fh)


def filters(locus_allele_info, locus, locus_reads, datatype):
    metadata = filter_marker_data[locus]
    if len(locus_allele_info) == 1:
        locus_allele_info = single_allele_thresholds(metadata, locus_reads, locus_allele_info)
    else:
        locus_allele_info, locus_reads = multiple_allele_thresholds(
            metadata, locus_reads, locus_allele_info
        )
        locus_allele_info = ce_filtering(locus_allele_info, locus_reads, metadata, datatype)
        if datatype == "ngs":
            locus_allele_info = same_size_filter(locus_allele_info, metadata)
    return locus_allele_info


def single_allele_thresholds(metadata, locus_reads, single_all_df):
    if thresholds("Detection", metadata, locus_reads, single_all_df["Reads"][0])[1] is False:
        single_all_df = pd.DataFrame()
    elif thresholds("Analytical", metadata, locus_reads, single_all_df["Reads"][0])[1] is False:
        single_all_df[["allele_type", "perc_noise"]] = ["noise", 1.0]
    elif thresholds("Analytical", metadata, locus_reads, single_all_df["Reads"][0])[1] is True:
        single_all_df["allele_type"] = "real_allele"
    return single_all_df


def multiple_allele_thresholds(metadata, locus_reads, locus_allele_info):
    for i in range(len(locus_allele_info)):  # check for alleles below detection threshold
        quest_al_reads = locus_allele_info.loc[i, "Reads"]
        if thresholds("Detection", metadata, locus_reads, quest_al_reads)[1] is False:
            locus_allele_info = locus_allele_info.drop(locus_allele_info.index[i])
            locus_reads = thresholds("Detection", metadata, locus_reads, quest_al_reads)[0]
    locus_allele_info = locus_allele_info.reset_index(drop=True)
    for i in range(len(locus_allele_info)):  # check for alleles below AT threshold
        quest_allele_reads = locus_allele_info.loc[i, "Reads"]
        if thresholds("Analytical", metadata, locus_reads, quest_allele_reads)[1] is False:
            locus_allele_info.loc[i, ["allele_type", "perc_noise"]] = [
                "noise",
                round(quest_allele_reads / locus_reads, 3),
            ]
        else:
            locus_allele_info.loc[i, "allele_type"] = "real_allele"
    return locus_allele_info, locus_reads


def thresholds(filter, metadata, locus_reads, quest_al_reads):
    use = metadata[f"{filter}ThresholdUse"]
    count = metadata[f"{filter}ThresholdStaticCount"]
    perc = metadata[f"{filter}ThresholdDynamicPercent"]
    thresh_perc = round(perc * locus_reads, 1)
    if (
        use.lower() == "dynamic"
        and locus_reads < metadata["MinimumNumberReadsForDynamicThresholds"]
    ):
        use = "static"
    if use.lower() == "both":
        if thresh_perc >= count:
            thresh = thresh_perc
        else:
            thresh = count
    elif use.lower() == "static":
        thresh = count
    elif use.lower() == "dynamic":
        thresh = thresh_perc
    else:
        print(
            f"Specified {filter} Threshold Use Not Acceptable!"
            f"Check filters.json file and re-run."
        )
    if quest_al_reads < thresh:
        if filter == "Detection":
            locus_reads = locus_reads - quest_al_reads
        return locus_reads, False
    else:
        return locus_reads, True


def ce_filtering(locus_allele_info, locus_reads, metadata, datatype):
    for i in range(len(locus_allele_info)):  # check for stutter alleles
        if locus_allele_info.loc[i, "allele_type"] != "real_allele":
            continue
        else:
            ref_allele_reads = locus_allele_info.loc[i, "Reads"]
            for j in range(len(locus_allele_info)):
                if j == i:
                    continue
                init_type_all = locus_allele_info.loc[j, "allele_type"]
                if init_type_all == "noise":
                    continue
                locus_allele_info = allele_ident(
                    locus_allele_info, init_type_all, metadata, ref_allele_reads, i, j, datatype
                )
        for j in range(len(locus_allele_info)):
            type_all = locus_allele_info.loc[j, "allele_type"]
            if "stutter" in type_all:
                if "/" not in type_all:
                    if pd.isnull(locus_allele_info.loc[j, "perc_stutter"]):
                        locus_allele_info.loc[j, "perc_stutter"] = round(
                            locus_allele_info.loc[j, "Reads"]
                            / locus_allele_info.loc[j, "allele1_ref_reads"],
                            3,
                        )
                else:
                    locus_allele_info.loc[j, "perc_stutter"] = ""
                locus_allele_info.loc[j, "perc_noise"] = ""
            elif "noise" in locus_allele_info.loc[j, "allele_type"]:
                locus_allele_info.loc[j, "perc_noise"] = round(
                    locus_allele_info.loc[j, "Reads"] / locus_reads, 3
                )
    return locus_allele_info


def allele_ident(locus_allele_info, init_type_all, metadata, ref_allele_reads, i, j, datatype):
    quest_al_reads = locus_allele_info.loc[j, "Reads"]
    ref_allele = float(locus_allele_info.loc[i, "RU_Allele"])
    question_allele = float(locus_allele_info.loc[j, "RU_Allele"])
    if datatype == "ngs":
        ref_bracket = locus_allele_info.loc[i, "UAS_Output_Bracketed_Notation"]
        question_bracket = locus_allele_info.loc[j, "UAS_Output_Bracketed_Notation"]
    else:
        ref_bracket = None
        question_bracket = None
    locus_allele_info.loc[j, ["allele_type", "perc_stutter"]] = allele_type(
        ref_allele,
        question_allele,
        init_type_all,
        metadata,
        quest_al_reads,
        ref_allele_reads,
        locus_allele_info.loc[j, "allele1_ref_reads"],
        locus_allele_info,
        ref_bracket,
        question_bracket,
        datatype,
    )
    if "stutter" in locus_allele_info.loc[j, "allele_type"]:
        if "/" in locus_allele_info.loc[j, "allele_type"] and pd.isnull(
            locus_allele_info.loc[j, "stuttering_allele2"]
        ):
            stut_allele2 = ref_allele if datatype == "ce" else ref_bracket
            locus_allele_info.loc[j, ["stuttering_allele2", "allele2_ref_reads"]] = [
                stut_allele2,
                ref_allele_reads,
            ]
        elif pd.isnull(locus_allele_info.loc[j, "stuttering_allele1"]):
            stut_allele1 = ref_allele if datatype == "ce" else ref_bracket
            locus_allele_info.loc[j, ["stuttering_allele1", "allele1_ref_reads"]] = [
                stut_allele1,
                ref_allele_reads,
            ]
    return locus_allele_info


def minus1_stutter(
    all_type,
    stutter_thresh,
    forward_thresh,
    stutter_thresh_reads,
    ref_reads,
    al1_ref_reads,
    quest_al_reads,
):
    stut_perc = None
    if all_type == "+1_stutter":
        all_thresh = stutter_thresh_reads + forward_stut_thresh(
            forward_thresh, stutter_thresh, al1_ref_reads
        )
        all_type = output_allele_call(quest_al_reads, all_thresh, "-1_stutter/+1_stutter")
    elif all_type == "-2_stutter":
        all_thresh = stutter_thresh_reads + np.ceil(stutter_thresh * al1_ref_reads)
        all_type = output_allele_call(quest_al_reads, all_thresh, "-1_stutter/-2_stutter")
    elif quest_al_reads <= stutter_thresh_reads:
        all_type = "-1_stutter"
        stut_perc = round(quest_al_reads / ref_reads, 3)
    return all_type, stut_perc


def minus2_stutter(
    all_type,
    stutter_thresh,
    forward_thresh,
    stutter_thresh_reads,
    ref_reads,
    al1_ref_reads,
    quest_al_reads,
):
    stut_perc = None
    if all_type == "-1_stutter":
        all_thresh = stutter_thresh_reads + np.ceil(stutter_thresh * al1_ref_reads)
        all_type = output_allele_call(quest_al_reads, all_thresh, "-1_stutter/-2_stutter")
    elif all_type == "+1_stutter":
        all_thresh = stutter_thresh_reads + forward_stut_thresh(
            forward_thresh, stutter_thresh, al1_ref_reads
        )
        all_type = output_allele_call(quest_al_reads, all_thresh, "+1_stutter/-2_stutter")
    elif quest_al_reads <= stutter_thresh_reads:
        all_type = "-2_stutter"
        stut_perc = round(quest_al_reads / ref_reads, 3)
    return all_type, stut_perc


def plus1_stutter(
    all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads, quest_al_reads
):
    stut_perc = None
    if all_type == "-1_stutter":
        all_thresh = np.ceil(stutter_thresh * al1_ref_reads) + forward_stut_thresh(
            forward_thresh, stutter_thresh, ref_reads
        )
        all_type = output_allele_call(quest_al_reads, all_thresh, "-1_stutter/+1_stutter")
    elif all_type == "-2_stutter":
        all_thresh = np.ceil(stutter_thresh * al1_ref_reads) + forward_stut_thresh(
            forward_thresh, stutter_thresh, ref_reads
        )
        all_type = output_allele_call(quest_al_reads, all_thresh, "+1_stutter/-2_stutter")
    elif quest_al_reads <= forward_stut_thresh(forward_thresh, stutter_thresh, ref_reads):
        all_type = "+1_stutter"
        stut_perc = round(quest_al_reads / ref_reads, 3)
    return all_type, stut_perc


def allele_type(
    ref,
    ru,
    all_type,
    metadata,
    quest_al_reads,
    ref_reads,
    al1_ref_reads,
    all_type_df,
    ref_bracket,
    question_bracket,
    datatype,
):
    stutter_thresh = metadata["StutterThresholdDynamicPercent"]
    forward_thresh = metadata["StutterForwardThresholdDynamicPercent"]
    stutter_thresh_reads = np.ceil(stutter_thresh * ref_reads)
    stut_perc = None
    allele_diff = round(ref - ru, 1)
    if allele_diff == 1 and ref_reads > quest_al_reads:  # -1 stutter
        if (
            datatype == "ngs" and bracketed_stutter_id(ref_bracket, question_bracket, -1) == -1
        ) or datatype == "ce":
            all_type, stut_perc = minus1_stutter(
                all_type,
                stutter_thresh,
                forward_thresh,
                stutter_thresh_reads,
                ref_reads,
                al1_ref_reads,
                quest_al_reads,
            )
    elif allele_diff == 2 and ref_reads > quest_al_reads:  # -2 stutter
        allele = ru if datatype == "ce" else question_bracket
        if check_2stutter(all_type_df, datatype, allele)[0] is True:
            if (
                datatype == "ngs" and bracketed_stutter_id(ref_bracket, question_bracket, -2) == -2
            ) or datatype == "ce":
                ref_reads = check_2stutter(all_type_df, datatype, allele)[1]
                stutter_thresh_reads = stutter_thresh * ref_reads
                all_type, stut_perc = minus2_stutter(
                    all_type,
                    stutter_thresh_reads,
                    forward_thresh,
                    stutter_thresh_reads,
                    ref_reads,
                    al1_ref_reads,
                    quest_al_reads,
                )
    elif allele_diff == -1 and ref_reads > quest_al_reads:  # +1 stutter
        if (
            datatype == "ngs" and bracketed_stutter_id(ref_bracket, question_bracket, 1) == 1
        ) or datatype == "ce":
            all_type, stut_perc = plus1_stutter(
                all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads, quest_al_reads
            )
    elif pd.isnull(all_type):
        all_type = "noise"
    return all_type, stut_perc


def output_allele_call(quest_al_reads, all_thresh, orig_type):
    if quest_al_reads <= all_thresh:
        all_type = orig_type
    else:
        all_type = "real_allele"
    return all_type


def forward_stut_thresh(perc, perc_stut, reads):
    if perc == 0:
        forward_thresh = np.ceil(perc_stut**2 * reads)
    else:
        forward_thresh = np.ceil(perc * reads)
    return forward_thresh


def check_2stutter(stutter_df, allele_des, allele):
    is_true, reads = False, None
    if "-1_stutter" in stutter_df.loc[:, "allele_type"].values:
        if allele_des == "ce":
            for k, row in stutter_df.iterrows():
                ru_test = stutter_df.loc[k, "RU_Allele"]
                if ru_test - allele == 1 and stutter_df.loc[k, "allele_type"] == "-1_stutter":
                    is_true, reads = True, stutter_df.loc[k, "Reads"]
                    break
        else:
            for k, row in stutter_df.iterrows():
                bracket_test = stutter_df.loc[k, "UAS_Output_Bracketed_Notation"]
                if bracketed_stutter_id(bracket_test, allele, -1) == -1:
                    is_true, reads = True, stutter_df.loc[k, "Reads"]
    return is_true, reads


def allele_counts(allele_df):
    mix_df = pd.DataFrame(columns=["SampleID", "Locus", "Flags"])
    try:
        if allele_df.allele_type.value_counts()["real_allele"] > 2:
            mix_df.loc[len(mix_df.index)] = [
                allele_df.loc[1, "SampleID"],
                allele_df.loc[1, "Locus"],
                ">2Alleles",
            ]
    except (KeyError, AttributeError):
        pass
    return mix_df


def create_repeat_dict(bracket):
    repeats = defaultdict(lambda: defaultdict(dict))
    num = 0
    for unit in bracket.split(" "):
        num += 1
        if "[" in unit:
            match = re.match(r"\[([ACGT]+)\](\d+)", unit)
            if match:
                repeat_motif = match.group(1)
                count = match.group(2)
        else:
            repeat_motif = unit
            count = 1
        repeats[num]["motif"] = repeat_motif
        repeats[num]["counts"] = int(count)
    return repeats


def bracketed_stutter_id(ref_bracket, quest_bracket, stutter_id):
    ref_repeats = create_repeat_dict(ref_bracket)
    quest_repeats = create_repeat_dict(quest_bracket)
    stutter = None
    if len(quest_repeats) == len(ref_repeats):
        diffcount = 0
        for key, value in ref_repeats.items():
            ref_info = ref_repeats[key]
            quest_info = quest_repeats[key]
            if ref_info["motif"] != quest_info["motif"]:
                stutter = None
                break
            else:
                diff = quest_info["counts"] - ref_info["counts"]
                if diff == 0:
                    continue
                elif diff == stutter_id:
                    if diffcount == 0:
                        diffcount = diff
                        stutter = stutter_id
                    else:
                        stutter = None
                        break
                else:
                    stutter = None
                    break
    return stutter


def allele_imbalance_check(allele_df):
    imbalance_df = pd.DataFrame(columns=["SampleID", "Locus", "Flags"])
    try:
        locus = allele_df.loc[0, "Locus"]
        metadata = filter_marker_data[locus]
        het_perc = metadata["MinimumHeterozygousBalanceThresholdDynamicPercent"]
        if allele_df.allele_type.value_counts()["real_allele"] >= 2:
            real_df = allele_df[allele_df["allele_type"] == "real_allele"].reset_index(drop=True)
            max_reads = real_df["Reads"].max()
            min_reads = real_df["Reads"].min()
            if min_reads / max_reads < het_perc:
                imbalance_df.loc[len(imbalance_df.index)] = [
                    real_df.loc[0, "SampleID"],
                    locus,
                    "IntraLocusImbalance",
                ]
    except (KeyError, AttributeError):
        pass
    return imbalance_df


def check_D7(allele_df):
    D7_df = pd.DataFrame(columns=["SampleID", "Locus", "Flags"])
    for i in range(len(allele_df)):
        al = allele_df.loc[i, "RU_Allele"]
        try:
            if str(al).split(".")[1] == "1" and allele_df.loc[i, "allele_type"] != "noise":
                print("D7 microvariants detected! Check flagged file for details.")
                D7_df.loc[len(D7_df.index)] = [
                    allele_df.loc[0, "SampleID"],
                    allele_df.loc[0, "Locus"],
                    "Microvariant",
                ]
        except IndexError:
            continue
    return D7_df


def flags(allele_df):
    flags_df = pd.DataFrame(columns=["SampleID", "Locus", "Flags"])
    flags_df = flags_df.append(allele_counts(allele_df))
    flags_df = flags_df.append(allele_imbalance_check(allele_df))
    flags_df = flags_df.append(check_D7(allele_df))
    return flags_df


def same_size_filter(df, metadata):
    final_df = pd.DataFrame()
    al_list = df["RU_Allele"].unique()
    for ru_allele in al_list:
        df_filt = (
            df[df["RU_Allele"] == ru_allele]
            .sort_values(by=["Reads"], ascending=False)
            .reset_index(drop=True)
        )
        if len(df_filt) > 1:
            high_reads = df_filt.loc[0, "Reads"]
            final_df = final_df.append(df_filt.loc[0, :])
            for i in range(1, len(df_filt)):
                size_thresh = metadata["SameSizeThresholdDynamicPercent"]
                size_thresh_reads = high_reads * size_thresh
                if df_filt.loc[i, "Reads"] >= size_thresh_reads:
                    final_df = final_df.append(df_filt.loc[i, :])
        else:
            final_df = final_df.append(df_filt)
    final_df = final_df.reset_index(drop=True)
    return final_df
