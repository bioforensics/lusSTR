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
import importlib.resources
import json
import numpy as np
import pandas as pd
import re


def get_filter_metadata_file():
    return importlib.resources.files("lusSTR") / "data/filters.json"


with open(get_filter_metadata_file(), "r") as fh:
    filter_marker_data = json.load(fh)


def filters(locus_allele_info, locus, locus_reads, datatype, brack_col):
    metadata = filter_marker_data[locus]
    if locus == "AMELOGENIN":
        locus_allele_info = filter_amel(metadata, locus_allele_info, locus_reads)
    else:
        locus_allele_info["CE_Allele"] = locus_allele_info["CE_Allele"].astype(float)
        if len(locus_allele_info) == 1:
            locus_allele_info = single_allele_thresholds(metadata, locus_reads, locus_allele_info)
        else:
            locus_allele_info, locus_reads = multiple_allele_thresholds(
                metadata, locus_reads, locus_allele_info
            )
            locus_allele_info = ce_filtering(
                locus_allele_info, locus_reads, metadata, datatype, brack_col
            )
            if datatype != "ce":
                locus_allele_info = same_size_filter(locus_allele_info, metadata, datatype)
    return locus_allele_info


def filter_amel(metadata, amel_df, locus_reads):
    for filter in ["Detection", "Analytical"]:
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
            thresh = thresh_perc if thresh_perc >= count else count
        elif use.lower() == "static":
            thresh = count
        elif use.lower() == "dynamic":
            thresh = thresh_perc
        if filter == "Detection":
            amel_dt = amel_df[amel_df["Reads"] >= thresh].reset_index(drop=True)
            locus_reads = amel_df["Reads"].sum()
        else:
            for i in range(len(amel_dt)):
                al_reads = amel_dt.loc[i, "Reads"]
                if al_reads < thresh:
                    amel_dt.loc[i, ["allele_type", "perc_noise"]] = [
                        "BelowAT",
                        round(al_reads / locus_reads, 3),
                    ]
                else:
                    amel_dt.loc[i, "allele_type"] = "Typed"
    return amel_dt


def single_allele_thresholds(metadata, locus_reads, single_all_df):
    if thresholds("Detection", metadata, locus_reads, single_all_df["Reads"][0])[1] is False:
        single_all_df = pd.DataFrame()
    elif thresholds("Analytical", metadata, locus_reads, single_all_df["Reads"][0])[1] is False:
        single_all_df[["allele_type", "perc_noise"]] = ["BelowAT", 1.0]
    elif thresholds("Analytical", metadata, locus_reads, single_all_df["Reads"][0])[1] is True:
        single_all_df["allele_type"] = "Typed"
    return single_all_df


def multiple_allele_thresholds(metadata, locus_reads, locus_allele_info):
    for i in range(len(locus_allele_info)):  # check for alleles below detection threshold
        quest_al_reads = locus_allele_info.loc[i, "Reads"]
        if thresholds("Detection", metadata, locus_reads, quest_al_reads)[1] is False:
            locus_allele_info = locus_allele_info.drop(index=i)
            locus_reads = thresholds("Detection", metadata, locus_reads, quest_al_reads)[0]
    locus_allele_info = locus_allele_info.reset_index(drop=True)
    for i in range(len(locus_allele_info)):  # check for alleles below AT threshold
        quest_allele_reads = locus_allele_info.loc[i, "Reads"]
        if thresholds("Analytical", metadata, locus_reads, quest_allele_reads)[1] is False:
            locus_allele_info.loc[i, ["allele_type", "perc_noise"]] = [
                "BelowAT",
                round(quest_allele_reads / locus_reads, 3),
            ]
        else:
            locus_allele_info.loc[i, "allele_type"] = "Typed"
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
        thresh = thresh_perc if thresh_perc >= count else count
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


def ce_filtering(locus_allele_info, locus_reads, metadata, datatype, brack_col):
    for i in range(len(locus_allele_info)):  # check for stutter alleles
        if locus_allele_info.loc[i, "allele_type"] != "Typed":
            continue
        else:
            ref_allele_reads = locus_allele_info.loc[i, "Reads"]
            for j in range(len(locus_allele_info)):
                if j == i:
                    continue
                init_type_all = locus_allele_info.loc[j, "allele_type"]
                if init_type_all == "BelowAT":
                    continue
                locus_allele_info = allele_ident(
                    locus_allele_info,
                    init_type_all,
                    metadata,
                    ref_allele_reads,
                    i,
                    j,
                    datatype,
                    brack_col,
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
                    locus_allele_info.loc[j, "perc_stutter"] = None
                locus_allele_info.loc[j, "perc_noise"] = None
            elif "BelowAT" in locus_allele_info.loc[j, "allele_type"]:
                locus_allele_info.loc[j, "perc_noise"] = round(
                    locus_allele_info.loc[j, "Reads"] / locus_reads, 3
                )
    return locus_allele_info


def allele_ident(
    locus_allele_info, init_type_all, metadata, ref_allele_reads, i, j, datatype, brack_col
):
    quest_al_reads = locus_allele_info.loc[j, "Reads"]
    if datatype == "lusplus":
        ref_allele = float(locus_allele_info.loc[i, "LUS_Plus"].split("_")[0])
        question_allele = float(locus_allele_info.loc[j, "LUS_Plus"].split("_")[0])
        ref_altformat = locus_allele_info.loc[i, "LUS_Plus"]
        question_altformat = locus_allele_info.loc[j, "LUS_Plus"]
    else:
        ref_allele = float(locus_allele_info.loc[i, "CE_Allele"])
        question_allele = float(locus_allele_info.loc[j, "CE_Allele"])
        ref_altformat = None
        question_altformat = None
        if datatype == "ngs":
            ref_altformat = locus_allele_info.loc[i, brack_col]
            question_altformat = locus_allele_info.loc[j, brack_col]
    locus_allele_info.loc[j, ["allele_type", "perc_stutter"]] = allele_type(
        ref_allele,
        question_allele,
        init_type_all,
        metadata,
        quest_al_reads,
        ref_allele_reads,
        locus_allele_info.loc[j, "allele1_ref_reads"],
        locus_allele_info,
        ref_altformat,
        question_altformat,
        datatype,
        brack_col,
    )
    if "stutter" in locus_allele_info.loc[j, "allele_type"]:
        if (
            "/" in locus_allele_info.loc[j, "allele_type"]
            and locus_allele_info.loc[j, "parent_allele2"] == "nan"
        ):
            stut_allele2 = ref_allele if datatype == "ce" else ref_altformat
            locus_allele_info.loc[j, ["parent_allele2", "allele2_ref_reads"]] = [
                stut_allele2,
                ref_allele_reads,
            ]
        elif locus_allele_info.loc[j, "parent_allele1"] == "nan":
            stut_allele1 = ref_allele if datatype == "ce" else ref_altformat
            locus_allele_info.loc[j, ["parent_allele1", "allele1_ref_reads"]] = [
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
    ce,
    all_type,
    metadata,
    quest_al_reads,
    ref_reads,
    al1_ref_reads,
    all_type_df,
    ref_altformat,
    question_altformat,
    datatype,
    brack_col,
):
    stutter_thresh = metadata["StutterThresholdDynamicPercent"]
    forward_thresh = metadata["StutterForwardThresholdDynamicPercent"]
    stutter_thresh_reads = np.ceil(stutter_thresh * ref_reads)
    stut_perc = None
    allele_diff = round(ref - ce, 1)
    if allele_diff == 1 and ref_reads > quest_al_reads:  # -1 stutter
        if (
            (
                datatype == "ngs"
                and bracketed_stutter_id(ref_altformat, question_altformat, -1) == -1
            )
            or datatype == "ce"
            or (
                datatype == "lusplus"
                and lusplus_stutter_id(ref_altformat, question_altformat, -1) == -1
            )
        ):
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
        allele = ce if datatype == "ce" else question_altformat
        if check_2stutter(all_type_df, datatype, allele, brack_col)[0] is True:
            if (
                (
                    datatype == "ngs"
                    and bracketed_stutter_id(ref_altformat, question_altformat, -2) == -2
                )
                or datatype == "ce"
                or (
                    datatype == "lusplus"
                    and lusplus_stutter_id(ref_altformat, question_altformat, -2) == -2
                )
            ):
                ref_reads = check_2stutter(all_type_df, datatype, allele, brack_col)[1]
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
            (datatype == "ngs" and bracketed_stutter_id(ref_altformat, question_altformat, 1) == 1)
            or datatype == "ce"
            or (
                datatype == "lusplus"
                and lusplus_stutter_id(ref_altformat, question_altformat, 1) == 1
            )
        ):
            all_type, stut_perc = plus1_stutter(
                all_type, stutter_thresh, forward_thresh, ref_reads, al1_ref_reads, quest_al_reads
            )
    elif pd.isnull(all_type):
        all_type = "BelowAT"
    return all_type, stut_perc


def output_allele_call(quest_al_reads, all_thresh, orig_type):
    if quest_al_reads <= all_thresh:
        all_type = orig_type
    else:
        all_type = "Typed"
    return all_type


def forward_stut_thresh(perc, perc_stut, reads):
    if perc == 0:
        forward_thresh = np.ceil(perc_stut**2 * reads)
    else:
        forward_thresh = np.ceil(perc * reads)
    return forward_thresh


def check_2stutter(stutter_df, allele_des, allele, brack_col):
    is_true, reads = False, None
    if "-1_stutter" in stutter_df.loc[:, "allele_type"].values:
        if allele_des == "ce":
            for k, row in stutter_df.iterrows():
                ce_test = stutter_df.loc[k, "CE_Allele"]
                if ce_test - allele == 1 and stutter_df.loc[k, "allele_type"] == "-1_stutter":
                    is_true, reads = True, stutter_df.loc[k, "Reads"]
                    break
        else:
            for k, row in stutter_df.iterrows():
                if allele_des == "ngs":
                    bracket_test = stutter_df.loc[k, brack_col]
                    if bracketed_stutter_id(bracket_test, allele, -1) == -1:
                        is_true, reads = True, stutter_df.loc[k, "Reads"]
                else:
                    lusp_test = stutter_df.loc[k, "LUS_Plus"]
                    if lusplus_stutter_id(lusp_test, allele, -1) == -1:
                        is_true, reads = True, stutter_df.loc[k, "Reads"]
    return is_true, reads


def allele_counts(allele_df):
    mix_df = pd.DataFrame(columns=["SampleID", "Locus", "Flags"])
    try:
        if allele_df.allele_type.value_counts()["Typed"] > 2:
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


def lusplus_stutter_id(ref_lusp, question_lusp, stutter_id):
    diffcount = 0
    stutter = None
    for j in range(1, len(ref_lusp.split("_"))):
        diff = float(question_lusp.split("_")[j]) - round(float(ref_lusp.split("_")[j]))
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
        if allele_df.allele_type.value_counts()["Typed"] >= 2:
            real_df = allele_df[allele_df["allele_type"] == "Typed"].reset_index(drop=True)
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


def check_D7(allele_df, datatype):
    D7_df = pd.DataFrame(columns=["SampleID", "Locus", "Flags"])
    for i in range(len(allele_df)):
        if datatype == "lusplus":
            al = allele_df.loc[i, "LUS_Plus"].split("_")[0]
        else:
            al = allele_df.loc[i, "CE_Allele"]
        try:
            if str(al).split(".")[1] == "1" and allele_df.loc[i, "allele_type"] != "BelowAT":
                print("D7 microvariants detected! Check flagged file for details.")
                D7_df.loc[len(D7_df.index)] = [
                    allele_df.loc[0, "SampleID"],
                    allele_df.loc[0, "Locus"],
                    "Microvariant",
                ]
        except IndexError:
            continue
    return D7_df


def flags(allele_df, datatype):
    flags_df = pd.DataFrame(columns=["SampleID", "Locus", "Flags"])
    flags_df = pd.concat([flags_df, allele_counts(allele_df)])
    flags_df = pd.concat([flags_df, allele_imbalance_check(allele_df)])
    flags_df = pd.concat([flags_df, check_D7(allele_df, datatype)])
    return flags_df


def same_size_filter(df, metadata, datatype):
    final_df = pd.DataFrame()
    if datatype == "lusplus":
        df["CE_Allele"] = df["LUS_Plus"].apply(lambda x: x.split("_")[0])
    al_list = df["CE_Allele"].unique()
    for ce_allele in al_list:
        df_filt = (
            df[df["CE_Allele"] == ce_allele]
            .sort_values(by=["Reads"], ascending=False)
            .reset_index(drop=True)
        )
        if len(df_filt) > 1:
            high_reads = df_filt.loc[0, "Reads"]
            final_df = pd.concat([final_df, pd.DataFrame([df_filt.loc[0, :]])])
            for i in range(1, len(df_filt)):
                size_thresh = metadata["SameSizeThresholdDynamicPercent"]
                size_thresh_reads = high_reads * size_thresh
                if df_filt.loc[i, "Reads"] >= size_thresh_reads:
                    final_df = pd.concat([final_df, pd.DataFrame([df_filt.loc[i, :]])])
        else:
            final_df = pd.concat([final_df, df_filt])
    final_df = final_df.reset_index(drop=True)
    return final_df
