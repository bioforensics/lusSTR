# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import importlib.resources
import json
import os
import pandas as pd
from pathlib import Path


def get_snp_metadata_file():
    return importlib.resources.files("lusSTR") / "data/snp_data.json"


with open(get_snp_metadata_file(), "r") as fh:
    snp_marker_data = json.load(fh)


def kintelligence_filtering(input):
    return input


def create_output_table(sample_df, orientation, separate, output_type, software):
    allele_des = {"A": "1", "C": "2", "G": "3", "T": "4"}
    if orientation == "uas":
        allele_col = "UAS_Allele"
    else:
        allele_col = "Forward_Strand_Allele"
    all_samples_df = pd.DataFrame()
    for sample in sample_df["SampleID"].unique():
        indiv_df = sample_df[sample_df["SampleID"] == sample]
        compiled_table = create_sample_df(indiv_df, output_type, allele_col)
        if software != "uas":
            compiled_table = check_allele_calls(compiled_table, output_type)
        compiled_table = compiled_table.replace(allele_des)
        compiled_table.insert(0, "Sample.Name", sample)
        all_samples_df = pd.concat([all_samples_df, compiled_table])
        if separate:
            Path(f"{output_type}_samples").mkdir(parents=True, exist_ok=True)
            if output_type == "evidence":
                separated_table = bin_snps(compiled_table, output_type, sample)
                separated_table.to_csv(
                    f"evidence_samples/{sample}_snpsetscombined_evidence.tsv",
                    index=False,
                    sep="\t",
                )
            compiled_table.to_csv(
                f"{output_type}_samples/{sample}_snp_{output_type}.tsv", index=False, sep="\t"
            )
    return all_samples_df


def bin_snps(sample_file, output_type, sample):
    height_cols = [col for col in sample_file.columns if "Height" in col]
    sample_file["Total_Reads"] = sample_file[height_cols].sum(axis=1)
    sorted_file = sample_file.sort_values(by=["Total_Reads", "Marker"])
    compiled_table = pd.DataFrame()
    for snp_num in range(0, 10):
        start = snp_num * 1000
        if snp_num != 9:
            end = start + 1000
            bin_df = sorted_file.iloc[
                start:end,
            ].reset_index(drop=True)
        else:
            bin_df = sorted_file.iloc[
                start : len(sorted_file),
            ].reset_index(drop=True)
        bin_df["Sample.Name"] = bin_df["Sample.Name"] + "_set" + str((snp_num + 1))
        compiled_table = pd.concat([compiled_table, bin_df])
        bin_df.to_csv(
            f"{output_type}_samples/{sample}_set{snp_num+1}.tsv",
            index=False,
            sep="\t",
        )
    return compiled_table


def create_sample_df(indiv_df, output_type, all_col):
    compiled_table = (
        indiv_df.groupby([indiv_df.groupby(["SampleID", "SNP"]).cumcount() + 1, "SNP"])
        .first()[[all_col, "Reads"]]
        .unstack(0)
        .reset_index()
    )
    try:
        compiled_table.columns = ["Marker", "Allele.1", "Allele.2", "Height.1", "Height.2"]
    except ValueError:
        try:
            compiled_table.columns = [
                "Marker",
                "Allele.1",
                "Allele.2",
                "Allele.3",
                "Height.1",
                "Height.2",
                "Height.3",
            ]
        except ValueError:
            compiled_table.columns = [
                "Marker",
                "Allele.1",
                "Allele.2",
                "Allele.3",
                "Allele.4",
                "Height.1",
                "Height.2",
                "Height.3",
                "Height.4",
            ]
            if len(compiled_table[compiled_table["Allele.4"].notna()]) > 0:
                compiled_table = compiled_table.drop(compiled_table.columns[[4, 8]], axis=1)
        if len(compiled_table[compiled_table["Allele.3"].notna()]) > 0:
            compiled_table = compiled_table.drop(compiled_table.columns[[3, 6]], axis=1)
    if output_type == "reference":
        for i, row in compiled_table.iterrows():
            if pd.isnull(compiled_table.loc[i, "Height.2"]):
                compiled_table.loc[i, "Allele.2"] = compiled_table.loc[i, "Allele.1"]
        compiled_table = compiled_table[["Marker", "Allele.1", "Allele.2"]]
    return compiled_table


def check_allele_calls(df, output_type):
    for i, row in df.iterrows():
        snpid = df.loc[i, "Marker"]
        marker_info = snp_marker_data[snpid]
        real_alleles = marker_info["Alleles"]
        if pd.isnull(df.loc[i, "Allele.2"]):
            if output_type == "evidence":
                if df.loc[i, "Allele.1"] == real_alleles[0]:
                    df.loc[i, "Allele.2"] = real_alleles[1]
                elif df.loc[i, "Allele.1"] == real_alleles[1]:
                    df.loc[i, "Allele.2"] = real_alleles[0]
                else:
                    print(f"{snpid} does not contain the correct alleles!")
                df.loc[i, "Height.2"] = 0
            else:
                df.loc[i, "Allele.2"] = df.loc[i, "Allele.1"]
    return df


def straitrazor_filtering(sr_df, thresh):
    thresh = float(thresh)
    for i, row in sr_df.iterrows():
        snpid = sr_df.loc[i, "SNP"]
        sampleid = sr_df.loc[i, "SampleID"]
        allele_reads = sr_df.loc[i, "Reads"]
        total_snp_reads = sr_df[(sr_df["SampleID"] == sampleid) & (sr_df["SNP"] == snpid)][
            "Reads"
        ].sum()
        if sr_df.loc[i, "Issues"] == "Allele call does not match expected allele!":
            sr_df = sr_df.drop(i)
        else:
            if allele_reads < (total_snp_reads * thresh):
                sr_df.loc[i, "Reads"] = 0
    return sr_df


def main(input, output, kit, strand, separate, refs, software, thresh):
    input = str(input)
    output_name = os.path.splitext(output)[0]
    input_file = pd.read_csv(input, sep="\t")
    if software == "uas":
        results = input_file
    else:
        results = straitrazor_filtering(input_file, thresh)
    if refs is None:
        refs = []
    if "," in refs:
        ref_ids = []
        ref_samples = pd.DataFrame()
        for ref in refs.split(","):
            ref = ref.strip()
            ref_samples = pd.concat([ref_samples, results[results["SampleID"].isin([ref])]])
            ref_ids.append(ref)
    else:
        ref_ids = [refs]
        ref_samples = results[results["SampleID"].isin([refs])]
    if len(ref_samples) > 0:
        ref_table = create_output_table(ref_samples, strand, separate, "reference", software)
        ref_table.to_csv(f"{output_name}_snp_reference.csv", index=False, sep="\t")
    evid_samples = results[~results.SampleID.isin(ref_ids)]
    if len(evid_samples) > 0:
        evid_table = create_output_table(evid_samples, strand, separate, "evidence", software)
        evid_table.to_csv(f"{output}_snp_evidence.csv", index=False, sep="\t")


if __name__ == "__main__":
    main(
        snakemake.input,
        output=snakemake.params.outputid,
        kit=snakemake.params.kit,
        strand=snakemake.params.strand,
        separate=snakemake.params.separate,
        refs=snakemake.params.refs,
        software=snakemake.params.software,
        thresh=snakemake.params.thresh,
    )
