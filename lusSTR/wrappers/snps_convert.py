import json
import os
import pandas as pd
from pathlib import Path
from pkg_resources import resource_filename


def get_snp_metadata_file():
    return resource_filename("lusSTR", "data/snp_data.json")


with open(get_snp_metadata_file(), "r") as fh:
    snp_marker_data = json.load(fh)


def kintelligence_filtering(input):
    return input


def create_output_table(sample_df, orientation, separate, output_type, uas):
    allele_des = {"A": "1", "C": "2", "G": "3", "T": "4"}
    if orientation == "uas":
        allele_col = "UAS_Allele"
    else:
        allele_col = "Forward_Strand_Allele"
    all_samples_df = pd.DataFrame()
    for sample in sample_df["SampleID"].unique():
        indiv_df = sample_df[sample_df["SampleID"] == sample]
        compiled_table = create_sample_df(indiv_df, output_type, allele_col)
        if not uas:
            compiled_table = check_allele_calls(compiled_table, output_type)
        compiled_table = compiled_table.replace(allele_des)
        compiled_table.insert(0, "Sample Name", sample)
        all_samples_df = all_samples_df.append(compiled_table)
        if separate:
            Path(f"{output_type}_samples").mkdir(parents=True, exist_ok=True)
            compiled_table.to_csv(
                f"{output_type}_samples/{sample}_snp_{output_type}.csv", index=False, sep="\t"
            )
    return all_samples_df


def create_sample_df(indiv_df, output_type, all_col):
    compiled_table = (
        indiv_df.groupby([indiv_df.groupby(["SampleID", "SNP"]).cumcount() + 1, "SNP"])
        .first()[[all_col, "Reads"]]
        .unstack(0)
        .reset_index()
    )
    print(compiled_table)
    try:
        compiled_table.columns = ["Marker", "Allele 1", "Allele 2", "Height 1", "Height 2"]
    except ValueError:
        print("Too many alleles!")
    if output_type == "reference":
        for i, row in compiled_table.iterrows():
            if compiled_table.loc[i, "Height 2"] == 0:
                compiled_table.loc[i, "Allele 2"] = compiled_table.loc[i, "Allele 1"]
        compiled_table = compiled_table[["Marker", "Allele 1", "Allele 2"]]
    return compiled_table


def check_allele_calls(df, output_type):
    for i, row in df.iterrows():
        snpid = df.loc[i, "Marker"]
        marker_info = snp_marker_data[snpid]
        real_alleles = marker_info["Alleles"]
        if pd.isnull(df.loc[i, "Allele 2"]):
            if output_type == "evidence":
                if df.loc[i, "Allele 1"] == real_alleles[0]:
                    df.loc[i, "Allele 2"] = real_alleles[1]
                elif df.loc[i, "Allele 1"] == real_alleles[1]:
                    df.loc[i, "Allele 2"] = real_alleles[0]
                else:
                    print(f"{snpid} does not contain the correct alleles!")
                df.loc[i, "Height 2"] = 0
            else:
                df.loc[i, "Allele 2"] = df.loc[i, "Allele 1"]
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


def main(input, output, kit, strand, separate, refs, uas, thresh):
    input = str(input)
    output_name = os.path.splitext(output)[0]
    input_file = pd.read_csv(input, sep="\t")
    if uas:
        results = input_file
    else:
        results = straitrazor_filtering(input_file, thresh)
    ref_samples = results[results["SampleID"].isin([refs])]
    if len(ref_samples) > 0:
        ref_table = create_output_table(ref_samples, strand, separate, "reference", uas)
        ref_table.to_csv(f"{output_name}_snp_reference.csv", index=False, sep="\t")
    evid_samples = results[~results["SampleID"].isin([refs])]
    if len(evid_samples) > 0:
        evid_table = create_output_table(evid_samples, strand, separate, "evidence", uas)
        evid_table.to_csv(f"{output}_snp_evidence.csv", index=False, sep="\t")


if __name__ == "__main__":
    main(
        snakemake.input,
        output=snakemake.params.outputid,
        kit=snakemake.params.kit,
        strand=snakemake.params.strand,
        separate=snakemake.params.separate,
        refs=snakemake.params.refs,
        uas=snakemake.params.uas,
        thresh=snakemake.params.thresh,
    )