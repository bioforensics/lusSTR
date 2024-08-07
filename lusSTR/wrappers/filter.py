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

import argparse
from collections import defaultdict
import importlib.resources
import json
import lusSTR
from lusSTR.scripts.filter_settings import filters, flags
import math
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from pathlib import Path
import re
import sys


strs = [
    "CSF1PO",
    "D10S1248",
    "D12S391",
    "D13S317",
    "D16S539",
    "D17S1301",
    "D18S51",
    "D19S433",
    "D1S1656",
    "D20S482",
    "D21S11",
    "D22S1045",
    "D2S1338",
    "D2S441",
    "D3S1358",
    "D4S2408",
    "D5S818",
    "D6S1043",
    "D7S820",
    "D8S1179",
    "D9S1122",
    "FGA",
    "PENTA D",
    "PENTA E",
    "TH01",
    "TPOX",
    "VWA",
    "DYS19",
    "DYS385A-B",
    "DYS389II",
    "DYS390",
    "DYS391",
    "DYS392",
    "DYS393",
    "DYS437",
    "DYS438",
    "DYS439",
    "DYS448",
    "DYS456",
    "DYS458",
    "DYS481",
    "DYS533",
    "DYS549",
    "DYS570",
    "DYS576",
    "DYS635",
    "DYS643",
    "Y-GATA-H4",
]


def get_filter_metadata_file():
    return importlib.resources.files("lusSTR") / "data/filters.json"


with open(get_filter_metadata_file(), "r") as fh:
    filter_marker_data = json.load(fh)


def process_strs(dict_loc, datatype, seq_col, brack_col):
    final_df = pd.DataFrame()
    flags_df = pd.DataFrame()
    for key, value in dict_loc.items():
        data = dict_loc[key].reset_index(drop=True)
        if datatype == "ce":
            data_combine = data.groupby(["SampleID", "Locus", "CE_Allele"], as_index=False)[
                "Reads"
            ].sum()
        elif datatype == "lusplus":
            data_combine = data.groupby(
                ["SampleID", "Locus", "CE_Allele", "LUS_Plus"], as_index=False
            )["Reads"].sum()
        else:
            data_combine = data[
                [
                    "SampleID",
                    "Locus",
                    seq_col,
                    "CE_Allele",
                    brack_col,
                    "Reads",
                ]
            ]
        data_order = data_combine.sort_values(by=["CE_Allele"], ascending=False).reset_index(
            drop=True
        )
        total_reads = data_order["Reads"].sum()
        locus = key[1]
        data_order = data_order.reindex(
            columns=[
                *data_order.columns.tolist(),
                "allele_type",
                "parent_allele1",
                "parent_allele2",
                "allele1_ref_reads",
                "allele2_ref_reads",
                "perc_noise",
                "perc_stutter",
            ],
            fill_value=None,
        )
        filtered_df = filters(data_order, locus, total_reads, datatype, brack_col)
        final_df = pd.concat([final_df, filtered_df])
        flags_df = pd.concat([flags_df, flags(filtered_df, datatype)])
    if datatype == "ce" or datatype == "ngs":
        final_df = final_df.astype({"CE_Allele": "float64", "Reads": "int"})
    return final_df, flags_df


def EFM_output(profile, outfile, profile_type, data_type, col, separate=False):
    if profile_type == "reference":
        profile = profile[profile.allele_type == "real_allele"]
    else:
        profile = profile[profile.allele_type != "BelowAT"]
    efm_profile = populate_efm_profile(profile, data_type, col)
    if separate:
        write_sample_specific_efm_profiles(efm_profile, profile_type, data_type, outfile)
    else:
        write_aggregate_efm_profile(efm_profile, profile_type, data_type, outfile)


def populate_efm_profile(profile, data_type, colname):
    if data_type == "ce":
        prof_col = "CE_Allele"
    elif data_type == "lusplus":
        prof_col = "LUS_Plus"
    elif data_type == "ngs":
        prof_col = colname
    else:
        message = (
            f"Incorrect data type {data_type} specified for EFM. Please choose either "
            "'ce' or 'lusplus'."
        )
        raise ValueError(message)
    profile = profile.sort_values(by=["SampleID", "Locus", prof_col])
    profile = profile.rename(columns={prof_col: "Allele"})
    allele_heights = defaultdict(lambda: defaultdict(dict))
    for i, row in profile.iterrows():
        if data_type == "ce":
            allele_heights[row.SampleID][row.Locus][float(row.Allele)] = int(row.Reads)
        else:
            allele_heights[row.SampleID][row.Locus][row.Allele] = int(row.Reads)
    max_num_alleles = determine_max_num_alleles(allele_heights)
    reformatted_profile = list()
    for sampleid, loci in allele_heights.items():
        for locusid, alleles in loci.items():
            allele_list, height_list = list(), list()
            for allele, height in alleles.items():
                allele_list.append(allele)
                height_list.append(height)
            while len(allele_list) < max_num_alleles:
                allele_list.append(None)
                height_list.append(None)
            entry = [sampleid, locusid] + allele_list + height_list
            reformatted_profile.append(entry)
    for sampleid in allele_heights:
        for locusid in strs:
            if locusid not in allele_heights[sampleid]:
                entry = [sampleid, locusid] + ([None] * max_num_alleles * 2)
                reformatted_profile.append(entry)
    allele_columns = [f"Allele{n + 1}" for n in range(max_num_alleles)]
    height_columns = [f"Height{n + 1}" for n in range(max_num_alleles)]
    column_names = ["SampleName", "Marker"] + allele_columns + height_columns
    efm_profile = pd.DataFrame(reformatted_profile, columns=column_names)
    for col in height_columns:
        efm_profile[col] = efm_profile[col].astype("Int64")
    efm_profile = efm_profile.sort_values(by=["SampleName", "Marker"])
    return efm_profile


def write_sample_specific_efm_profiles(efm_profile, profile_type, data_type, outdir):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    for sample in efm_profile.SampleName:
        sample_profile = efm_profile[efm_profile.SampleName == sample].reset_index(drop=True)
        sample_profile.dropna(axis=1, how="all", inplace=True)
        if profile_type == "evidence":
            sample_profile.to_csv(f"{outdir}/{sample}_evidence_{data_type}.csv", index=False)
        else:
            num_alleles = (len(sample_profile.columns) - 2) / 2
            if num_alleles > 2:
                message = (
                    f"reference profile {sample} has at least one locus with {num_alleles} "
                    "alleles; stubbornly refusing to proceed with more than two alleles for "
                    "a reference profile"
                )
                raise ValueError(message)
            for i in range(len(sample_profile)):
                if pd.isna(sample_profile.loc[i, "Allele2"]):
                    sample_profile.loc[i, "Allele2"] = sample_profile.loc[i, "Allele1"]
            sample_profile.iloc[:, :4].to_csv(
                f"{outdir}/{sample}_reference_{data_type}.csv", index=False
            )


def write_aggregate_efm_profile(efm_profile, profile_type, data_type, outfile):
    Path(outfile).mkdir(parents=True, exist_ok=True)
    name = os.path.basename(outfile)
    if profile_type == "evidence":
        efm_profile.to_csv(f"{outfile}/{name}_evidence_{data_type}.csv", index=False)
    else:
        for i in range(len(efm_profile)):
            if pd.isna(efm_profile.loc[i, "Allele2"]):
                efm_profile.loc[i, "Allele2"] = efm_profile.loc[i, "Allele1"]
        efm_profile.iloc[:, :4].to_csv(f"{outfile}/{name}_reference_{data_type}.csv", index=False)


def determine_max_num_alleles(allele_heights):
    max_num_alleles = 0
    for sampleid, loci in allele_heights.items():
        for locusid, alleles in loci.items():
            if len(alleles) > max_num_alleles:
                max_num_alleles = len(alleles)
    return max_num_alleles


def STRmix_output(profile, outdir, profile_type, data_type, seq_col):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    if profile_type == "reference":
        filtered_df = profile[profile.allele_type == "real_allele"]
    else:
        filtered_df = profile[profile.allele_type != "BelowAT"]
    if data_type == "ce":
        strmix_profile = strmix_ce_processing(filtered_df)
    elif data_type == "lusplus":
        error_message = (
            "LUSPlus specified as the data type. STRmix does not accept LUSPlus! Please rerun with "
            "'ce' or 'ngs' specified."
        )
        raise ValueError(error_message)
    else:
        strmix_profile = filtered_df.loc[:, ["SampleID", "Locus", "CE_Allele", seq_col, "Reads"]]
        strmix_profile.rename(
            {"CE_Allele": "CE Allele", seq_col: "Allele Seq"}, axis=1, inplace=True
        )
        strmix_profile = strmix_profile.sort_values(by=["SampleID", "Locus", "CE Allele"])
    strmix_profile.replace(
        {"Locus": {"VWA": "vWA", "PENTA D": "PentaD", "PENTA E": "PentaE"}}, inplace=True
    )
    Path(outdir).mkdir(exist_ok=True)
    id_list = strmix_profile["SampleID"].unique()
    for id in id_list:
        sample_df = strmix_profile[strmix_profile["SampleID"] == id].reset_index(drop=True)
        if profile_type == "evidence":
            sample_df.iloc[:, 1:].to_csv(f"{outdir}/{id}_evidence_{data_type}.csv", index=False)
        else:
            reference_df = reference_table(sample_df, data_type)
            reference_df.to_csv(f"{outdir}/{id}_reference_{data_type}.csv", index=False)


def strmix_ce_processing(profile):
    data_combine = profile.groupby(["SampleID", "Locus", "CE_Allele"], as_index=False)[
        "Reads"
    ].sum()
    dict_loc = {k: v for k, v in data_combine.groupby(["SampleID", "Locus"])}
    locus_df = pd.DataFrame()
    for key, value in dict_loc.items():
        data = dict_loc[key].reset_index(drop=True)
        metadata = filter_marker_data[key[1]]
        slope = metadata["Slope"]
        intercept = metadata["Intercept"]
        data["Size"] = data["CE_Allele"] * slope + intercept
        locus_df = pd.concat([locus_df, data])
    locus_df.rename({"CE_Allele": "Allele", "Reads": "Height"}, axis=1, inplace=True)
    return locus_df


def reference_table(sample_data, datatype):
    new_rows = []
    for i, row in sample_data.iterrows():
        locus = sample_data.loc[i, "Locus"]
        try:
            next_col = sample_data.loc[i + 1, "Locus"]
        except KeyError:
            next_col = None
        try:
            prev_col = sample_data.loc[i - 1, "Locus"]
        except KeyError:
            prev_col = None
        if next_col == prev_col:
            message = (
                f"reference profile has at least one locus with > 2 alleles; "
                "stubbornly refusing to proceed with more than two alleles for "
                "a reference profile"
            )
            raise ValueError(message)
        elif locus == next_col or locus == prev_col:
            continue
        else:
            new_rows.append(list(row))
        final_reference = format_ref_table(new_rows, sample_data, datatype)
    return final_reference


def format_ref_table(new_rows, sample_data, datatype):
    if datatype == "ce":
        ref_filtered = pd.DataFrame(
            new_rows, columns=["SampleID", "Locus", "Allele", "Height", "Size"]
        )
        concat_df = pd.concat([sample_data.iloc[:, 1:3], ref_filtered.iloc[:, 1:3]]).reset_index(
            drop=True
        )
        sort_df = concat_df.sort_values(by=["Locus", "Allele"])
    else:
        ref_filtered = pd.DataFrame(
            new_rows, columns=["SampleID", "Locus", "CE Allele", "Allele Seq", "Reads"]
        )
        concat_df = pd.concat([sample_data.iloc[:, 1:4], ref_filtered.iloc[:, 1:4]]).reset_index(
            drop=True
        )
        sort_df = concat_df.sort_values(by=["Locus", "CE Allele"])
    return sort_df


def marker_plots(df, output_name, sex):
    Path("MarkerPlots").mkdir(parents=True, exist_ok=True)
    df["CE_Allele"] = df["CE_Allele"].astype(float)
    filt_df = df[df["allele_type"] == "real_allele"]
    for sample_id in df["SampleID"].unique():
        # sample_id = f"{id}_ystrs" if sex else id
        with PdfPages(f"MarkerPlots/{output_name}_{sample_id}_marker_plots.pdf") as pdf:
            make_plot(filt_df, sample_id, filters=True, at=False)
            pdf.savefig()
            make_plot(df, sample_id)
            pdf.savefig()
            make_plot(df, sample_id, sameyaxis=True)
            pdf.savefig()


def make_plot(df, sample_id, sameyaxis=False, filters=False, at=True):
    sample_df = df[df["SampleID"] == sample_id].copy()
    # sample_id = f"{id}_sexchr" if sex else id
    conditions = [
        sample_df["allele_type"].str.contains("real"),
        sample_df["allele_type"].str.contains("BelowAT"),
        sample_df["allele_type"].str.contains("stutter"),
    ]
    values = ["Typed", "BelowAT", "Stutter"]
    sample_df.loc[:, "Type"] = np.select(conditions, values)
    max_reads = max(sample_df["Reads"])
    n = 100 if max_reads > 1000 else 10
    max_yvalue = int(math.ceil(max_reads / n)) * n
    increase_value = int(math.ceil((max_yvalue / 5)) / n) * n
    fig = plt.figure(figsize=(30, 30))
    n = 0
    for marker in sample_df["Locus"].unique():
        n += 1
        colors = {"Typed": "g", "Stutter": "b", "BelowAT": "r"}
        marker_df = sample_df[sample_df["Locus"] == marker].sort_values(by="CE_Allele")
        ax = fig.add_subplot(6, 5, n)
        p = ax.bar(
            marker_df["CE_Allele"],
            marker_df["Reads"],
            edgecolor="k",
            color=[colors[i] for i in marker_df["Type"]],
        )
        if at:
            at = get_at(marker_df, marker)
            ax.axhline(at, linestyle="--", color="k")
            ax.text(round(min(marker_df["CE_Allele"])) - 0.9, at + (at * 0.1), f"AT", size=12)
        labels = marker_df["Type"].unique()
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors[l]) for l in labels]
        if not filters:
            plt.legend(handles, labels, title="Allele Type")
        else:
            for i, row in marker_df.iterrows():
                marker_df.loc[i, "Label"] = (
                    str(int(marker_df.loc[i, "CE_Allele"]))
                    if ".0" in str(marker_df.loc[i, "CE_Allele"])
                    else str(marker_df.loc[i, "CE_Allele"])
                )
            ax.bar_label(p, labels=marker_df["Label"])
        if sameyaxis:
            ax.set_yticks(np.arange(0, max_yvalue, increase_value))
        ax.set_xticks(
            np.arange(
                round(min(marker_df["CE_Allele"]) - 1), round(max(marker_df["CE_Allele"])) + 2, 1.0
            )
        )
        ax.title.set_text(marker)
    if sameyaxis:
        title = "Marker Plots for All Alleles With Same Y-Axis Scale"
    elif filters:
        title = "Marker Plots for Typed Alleles With Custom Y-Axis Scale"
    else:
        title = "Marker Plots for All Alleles With Custom Y-Axis Scale"
    plt.text(0.4, 0.95, title, transform=fig.transFigure, size=24)


def get_at(df, locus):
    metadata = filter_marker_data[locus]
    thresh_use = metadata["AnalyticalThresholdUse"]
    at_st = float(metadata["AnalyticalThresholdStaticCount"])
    at_dy = metadata["AnalyticalThresholdDynamicPercent"]
    at_dy_num = df["Reads"].sum() * float(at_dy)
    if thresh_use.lower() == "both":
        at = at_st if at_st > at_dy_num else at_dy_num
    elif thresh_use.lower() == "static":
        at = at_st
    elif thresh_use.lower() == "dynamic":
        at = at_dy_num
    else:
        raise ValueError("Incorrect AT specified in filters.json. Please check and re-run.")
    return at


def process_input(
    input_name,
    outpath,
    profile_type,
    data_type,
    output_type,
    nofilters,
    separate,
    custom,
    sex,
    info,
):
    full_df = pd.read_csv(f"{input_name}.txt", sep="\t")
    if custom:
        seq_col = "Custom_Range_Sequence"
        brack_col = "Custom_Bracketed_Notation"
    else:
        seq_col = "UAS_Output_Sequence" if strand == "uas" else "Forward_Strand_Sequence"
        brack_col = (
            "UAS_Output_Bracketed_Notation"
            if strand == "uas"
            else "Forward_Strand_Bracketed_Notation"
        )
    if nofilters:
        full_df["allele_type"] = "real_allele"
        marker_plots(full_df, input_name, sex)
        if output_type == "efm" or output_type == "mpsproto":
            EFM_output(full_df, outpath, profile_type, data_type, brack_col, separate)
        else:
            STRmix_output(full_df, outpath, profile_type, data_type, seq_col)
    else:
        dict_loc = {k: v for k, v in full_df.groupby(["SampleID", "Locus"])}
        final_df, flags_df = process_strs(dict_loc, data_type, seq_col, brack_col)
        marker_plots(final_df, input_name, sex)
        if output_type == "efm" or output_type == "mpsproto":
            EFM_output(final_df, outpath, profile_type, data_type, brack_col, separate)
        else:
            STRmix_output(final_df, outpath, profile_type, data_type, seq_col)
        if info:
            name = os.path.basename(outpath)
            final_df.to_csv(f"{outpath}/{input_name}_sequence_info.csv", index=False)
            if not flags_df.empty:
                flags_df.to_csv(f"{outpath}/{input_name}_Flagged_Loci.csv", index=False)


def main(
    input,
    output_type,
    profile_type,
    data_type,
    output_dir,
    info,
    separate,
    nofilters,
    strand,
    custom,
    sex,
):
    input = str(input)
    if profile_type not in ("evidence", "reference"):
        raise ValueError(f"unknown profile type '{profile_type}'")
    if data_type not in ("ce", "ngs", "lusplus"):
        raise ValueError(f"unknown data type '{data_type}'")
    if output_type not in ("efm", "strmix", "mpsproto"):
        raise ValueError(f"unknown output type '{output_type}'")
    if output_dir is None:
        raise ValueError("No output specified using --out.")
    if sex:
        outpath_sex = f"{output_dir}/ystrs/"
        input_name_sex = f"{os.path.splitext(input)[0]}_sexloci"
        process_input(
            input_name_sex,
            outpath_sex,
            profile_type,
            data_type,
            output_type,
            nofilters,
            separate,
            custom,
            sex,
            info,
        )
    input_name = os.path.splitext(input)[0]
    outpath = output_dir
    process_input(
        input_name,
        outpath,
        profile_type,
        data_type,
        output_type,
        separate,
        nofilters,
        custom,
        sex,
        info,
    )


if __name__ == "__main__":
    main(
        snakemake.input,
        output_type=snakemake.params.output_type,
        profile_type=snakemake.params.profile_type,
        data_type=snakemake.params.data_type,
        output_dir=snakemake.params.output_dir,
        info=snakemake.params.info,
        separate=snakemake.params.separate,
        nofilters=snakemake.params.filters,
        strand=snakemake.params.strand,
        custom=snakemake.params.custom,
        sex=snakemake.params.sex,
    )
