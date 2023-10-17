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

import glob
import json
import lusSTR
import pandas as pd
import openpyxl
import os
from pkg_resources import resource_filename


def get_snp_metadata_file():
    return resource_filename("lusSTR", "data/snp_data.json")


with open(get_snp_metadata_file(), "r") as fh:
    snp_marker_data = json.load(fh)


snp_type_dict = {"a": "Ancestry", "i": "Identity", "p": "Phenotype", "p/a": "Phenotype;Ancestry"}


snps_within_loci = {
    "mh16-MC1RB": {"SNPs": ["rs1805005", "rs1805006", "rs2228479"]},
    "mh16-MC1RC": {
        "SNPs": [
            "rs11547464",
            "rs1805007",
            "rs201326893_Y152OCH",
            "rs1110400",
            "rs1805008",
            "rs885479",
        ]
    },
}


def complement_base(base):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    comp_base = complement[base]
    return comp_base


def uas_format(infile, snp_type_arg, nofilter):
    """
    This function begins with the compiled data from all files within the specified directory.
    It removes any allele with Reads of 0; identifies whether the allele call needs to be reverse
    complemented to be reported on the forward strand; and checks that the called allele is one of
    two expected alleles for the SNP (and flags any SNP call which is unexpected).
    """
    data_df = uas_load(infile, nofilter, snp_type_arg)
    data_df.columns = [
        "SampleID",
        "Project",
        "Analysis",
        "SNP",
        "Reads",
        "Forward_Strand_Allele",
        "UAS_Allele",
        "Type",
        "Issues",
    ]
    data_dropdups = data_df.drop_duplicates()
    data_final = data_dropdups.sort_values(
        by=["SampleID", "Project", "Analysis", "SNP", "Reads"], ascending=False
    ).reset_index(drop=True)
    return data_final


def uas_load(input, nofilter, type="i"):
    """
    This function lists input .xlsx files within the specified directory and performs a check to
    ensure the correct file is processed (must contain either "Phenotype" or "Sample Details").
    This also compiles the SNP data for each file within the directory.
    """
    if os.path.isfile(input):
        snp_final_output = uas_types(input, type, nofilter)
    else:
        snp_final_output = pd.DataFrame()
        files = glob.glob(os.path.join(input, "[!~]*.xlsx"))
        for filename in sorted(files):
            print(filename)
            if "Phenotype" in filename or "Sample" in filename:
                snps = uas_types(filename, type, nofilter)
            if snps is not None:
                snp_final_output = snp_final_output.append(snps)
            else:
                continue
    return snp_final_output


def uas_types(infile, snp_type, nofilter):
    """
    This function determines which tab within the specified file is required to extract the SNP
    data from based on the name of the file.
    """
    if "Sample Details" in infile and (snp_type == "all" or snp_type == "i"):
        snp_data = parse_snp_table_from_sheet(infile, "iSNPs", snp_type, nofilter)
    elif "Phenotype" in infile and (snp_type == "all" or "a" in snp_type or "p" in snp_type):
        snp_data = parse_snp_table_from_sheet(infile, "SNP Data", snp_type, nofilter)
    elif "Sample Report" in infile:
        snp_data = process_kin(infile, nofilter)
    else:
        snp_data = None
    return snp_data


def parse_snp_table_from_sheet(infile, sheet, snp_type_arg, nofilter):
    """
    This function formats the SNP data from the original file and filters the SNPs based on the
    indicated SNP type.
    """
    file = openpyxl.load_workbook(infile)
    file_sheet = file[sheet]
    table = pd.DataFrame(file_sheet.values)
    offset = table[table.iloc[:, 0] == "Coverage Information"].index.tolist()[0]
    data = table.iloc[offset + 2 :]
    data.columns = table.iloc[offset + 1]
    data = data[["Locus", "Reads", "Allele Name", "Typed Allele?"]]
    if nofilter:
        data_typed = data
    else:
        data_typed = data[data["Typed Allele?"] == "Yes"]
    concat_df = pd.DataFrame()
    if snp_type_arg == "all":
        concat_df = data_typed
    else:
        for snp_type in snp_type_arg:
            filtered_dict = {k: v for k, v in snp_marker_data.items() if snp_type in v["Type"]}
            filtered_data = data_typed[data_typed["Locus"].isin(filtered_dict)].reset_index(
                drop=True
            )
            concat_df = concat_df.append(filtered_data)
    sampleid = table.iloc[2, 1]
    projectid = table.iloc[3, 1]
    analysisid = table.iloc[4, 1]
    final_df = process_sigprep_snps(
        concat_df.reset_index(drop=True), sampleid, projectid, analysisid
    )
    return final_df


def process_sigprep_snps(foren_df, sampid, projid, analyid):
    data_list = []
    for j, row in foren_df.iterrows():
        snpid = foren_df.loc[j, "Locus"]
        metadata = snp_marker_data[snpid]
        type = metadata["Type"]
        uas_allele = foren_df.loc[j, "Allele Name"]
        forward_strand_allele = check_rev_comp(uas_allele, snpid, metadata)
        flag = flags(foren_df, forward_strand_allele, j, metadata)
        row_tmp = [
            sampid,
            projid,
            analyid,
            snpid,
            foren_df.loc[j, "Reads"],
            forward_strand_allele,
            uas_allele,
            snp_type_dict[type],
            flag,
        ]
        data_list.append(row_tmp)
    return data_list


def process_kin(input, nofilter):
    """
    This function processes the Kintelligence Sample Report.
    """
    file = openpyxl.load_workbook(input)
    sheet_names = ["Ancestry SNPs", "Phenotype SNPs", "Identity SNPs", "Kinship SNPs"]
    data_filt = pd.DataFrame()
    uas_version = determine_version(file)
    for sheet in sheet_names:
        file_sheet = file[sheet]
        table = pd.DataFrame(file_sheet.values)
        if uas_version == "2.5.0":
            data = process_v5(table)
        else:
            data = process_v0(table)
        print(data)
        data = data[["Locus", "Reads", "Allele Name", "Typed Allele?"]]
        if nofilter:
            data_typed = data
        else:
            data_typed = data[data["Typed Allele?"] == "Yes"]
        print(data_typed)
        data_filt = data_filt.append(data_typed).reset_index(drop=True)
    sampid = table.iloc[2, 1]
    projid = table.iloc[3, 1]
    analyid = table.iloc[4, 1]
    data_df = []
    print(data_filt)
    for j, row in data_filt.iterrows():
        tmp_row = create_row(data_filt, j, sampid, projid, analyid)
        data_df.append(tmp_row)
    data_final = pd.DataFrame(
        data_df,
        columns=[
            "SampleID",
            "Project",
            "Analysis",
            "SNP",
            "Reads",
            "Forward_Strand_Allele",
            "UAS_Allele",
            "Type",
            "Issues",
        ],
    )
    data_final_sort = data_final.sort_values(
        by=["Project", "SampleID", "SNP", "Reads"], ascending=False
    ).reset_index(drop=True)
    return data_final_sort


def determine_version(file):
    file_sheet = file["Settings"]
    table = pd.DataFrame(file_sheet.values)
    try:
        version = table.loc[table[0] == "Software Version", 1].iloc[0]
    except IndexError:
        version = 2.0
    return version


def process_v0(table):
    offset = table[table.iloc[:, 0] == "Locus"].index.tolist()[0]
    data = table.iloc[offset + 1 :]
    data.columns = table.iloc[offset]
    return data


def process_v5(table):
    offset = table[table.iloc[:, 4] == "Locus"].index.tolist()[0]
    data = table.iloc[offset + 1 :, 4:8]
    data.columns = ["Locus", "Allele Name", "Typed Allele?", "Reads"]
    return data


def create_row(df, j, sampleid, projectid, analysisid):
    """
    This function first it identifies the Sig Prep SNPs (reverse complements those SNPs if
    neccesary and checks SNP allele calls for incorrect allele calls), and reports all SNP calls
    in the same general format.
    """
    snpid = df.loc[j, "Locus"]
    uas_allele = df.loc[j, "Allele Name"]
    print(snpid)
    try:
        metadata = snp_marker_data[snpid]
        forward_strand_allele = check_rev_comp(uas_allele, snpid, metadata)
        flag = flags(df, forward_strand_allele, j, metadata)
        type = snp_type_dict[metadata["Type"]]
    except KeyError:
        forward_strand_allele = uas_allele
        if df.loc[j, "Typed Allele?"] == "No":
            flag = "Contains untyped allele"
        else:
            flag = None
        type = "Kintelligence"
    new_row = [
        sampleid,
        projectid,
        analysisid,
        snpid,
        df.loc[j, "Reads"],
        forward_strand_allele,
        uas_allele,
        type,
        flag,
    ]
    return new_row


def check_rev_comp(allele, snpid, metadata):
    if metadata["ReverseCompNeeded"] == "Yes":
        f_strand_allele = complement_base(allele)
    else:
        f_strand_allele = allele
    return f_strand_allele


def flags(df, allele, j, metadata):
    """
    Checks that allele call is one of two known alleles for that SNP.
    """
    if df.loc[j, "Typed Allele?"] == "No":
        flag = "Contains untyped allele"
    elif allele in metadata["Alleles"]:
        flag = None
    else:
        flag = "Allele call does not match expected allele!"
    return flag


def strait_razor_format(infile, snp_type_arg):
    """
    This function formats STRait Razor input data for two separate reports. The full output
    includes all reads, the SNP allele calls and any results flags. In the main report, the reads
    are summed for identical allele calls per SNP. This function also checks that the allele call
    is one of two expected alleles for the SNP (and flags the allele if not).
    """
    results = strait_razor_concat(infile, snp_type_arg)
    results_sort = results.sort_values(
        by=["SampleID", "Project", "Analysis", "SNP", "Reads"], ascending=False
    )
    results_combine = results_sort.groupby(
        ["SNP", "Forward_Strand_Allele", "UAS_Allele", "Type", "SampleID", "Project", "Analysis"],
        as_index=False,
    )["Reads"].sum()
    results_combine = results_combine[
        [
            "SampleID",
            "Project",
            "Analysis",
            "SNP",
            "Reads",
            "Forward_Strand_Allele",
            "UAS_Allele",
            "Type",
        ]
    ]
    results_combine.loc[:, "Issues"] = None
    for j, row in results_combine.iterrows():
        snpid = results_combine.iloc[j, 3]
        metadata = snp_marker_data[snpid]
        if results_combine.iloc[j, 5] not in metadata["Alleles"]:
            results_combine.iloc[j, 8] = "Allele call does not match expected allele!"
    results_combine_sort = results_combine.sort_values(
        by=["SampleID", "Project", "Analysis", "SNP", "Reads"], ascending=False
    )
    return results_sort, results_combine_sort


def strait_razor_concat(input, snp_type_arg):
    """
    This function reads in all .txt files within the specified directory. For each file, the
    forward and reverse reads are summed and each sequence is processed and compiled into one
    final dataframe.
    """
    if os.path.isdir(input):
        analysisID = os.path.basename(input.rstrip(os.sep))
        files = glob.glob(os.path.join(input, "[!~]*.txt"))
        all_snps = pd.DataFrame()
        for filename in sorted(files):
            snps = read_straitrazor_data(filename, snp_type_arg, analysisID)
            all_snps = all_snps.append(snps)
    else:
        all_snps = read_straitrazor_data(input, snp_type_arg, None)
    all_snps.columns = [
        "SampleID",
        "Project",
        "Analysis",
        "SNP",
        "Sequence",
        "Reads",
        "Forward_Strand_Allele",
        "UAS_Allele",
        "Type",
        "Issues",
    ]
    return all_snps


def read_straitrazor_data(filename, snp_type_arg, analysisid):
    name = filename.replace(".txt", "").split(os.sep)[-1]
    if analysisid is None:
        analysisid = name
    table = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=["Locus_allele", "Length", "Sequence", "Forward_Reads", "Reverse_Reads"],
    )
    row = process_straitrazor_data(table, snp_type_arg, name, analysisid)
    snps = pd.DataFrame()
    if row is not None:
        snps = snps.append(row)
    return snps


def process_straitrazor_data(table, snp_type_arg, name, analysisid):
    try:
        snp_df = pd.DataFrame()
        table[["SNP", "Bases_off"]] = table.Locus_allele.str.split(":", expand=True)
        table["Total_Reads"] = table["Forward_Reads"] + table["Reverse_Reads"]
        snps_only = pd.DataFrame(table[table["SNP"].str.contains("rs|mh16|insA")]).reset_index(
            drop=True
        )
        for j, row in snps_only.iterrows():
            snpid = snps_only.loc[j, "SNP"]
            try:
                row = compile_row_of_snp_data(snps_only, snpid, j, snp_type_arg, name, analysisid)
            except KeyError:
                row = None
            snp_df = snp_df.append(row)
    except:
        print(
            f"Error found with {name}. Will bypass and continue. Please check file"
            f" and rerun the command, if necessary."
        )
        snp_df = None
    return snp_df


def compile_row_of_snp_data(infile, snp, table_loc, type, name, analysis):
    """
    This function is necessary to account for the two sets of SNPs reported from the same
    sequence amplicon. Sequences labeled as mh16-MC1RB and mh16-MC1RC contain 3 and 6 SNPs,
    respectively. This function reports out each SNP from the sequence amplicon as individual
    rows and calls another function to compile data on each SNP.
    """
    snp_df = pd.DataFrame()
    if "mh16" in snp:
        locus_data = snps_within_loci[snp]
        for k in range(0, len(locus_data["SNPs"])):
            snp_id = locus_data["SNPs"][k]
            row_tmp = collect_snp_info(infile, snp_id, table_loc, type, name, analysis)
            if row_tmp is not None:
                snp_df = snp_df.append(row_tmp)
    else:
        row_tmp = collect_snp_info(infile, snp, table_loc, type, name, analysis)
        if row_tmp is not None:
            snp_df = snp_df.append(row_tmp)
    return snp_df


def collect_snp_info(infile, snpid, j, allowed_snptype, name, analysis):
    """
    This function compiles allele calls, reads, reverse complements allele call if necessary to
    match how the UAS reports the allele, and any flags associated with the allele call. The flags
    indicate potential issues, including an unexpected allele call (not one of two expected
    alleles for the SNP) or unexpected length of the sequence amplicon which could result in an
    incorrect allele call. This function also determines if the SNP should be included in the
    final table based on the specified SNP type from the CLI.
    """
    if snpid == "N29insA":
        snpid = "rs312262906_N29insA"
    metadata = snp_marker_data[snpid]
    current_snp_type = metadata["Type"]
    seq = infile.loc[j, "Sequence"]
    expected_alleles = metadata["Alleles"]
    snp_loc = metadata["Coord"]
    all_rows = []
    if len(seq) > snp_loc:
        snp_call = seq[snp_loc]
        if snpid == "rs312262906_N29insA" and snp_call == "A":
            snp_call = "insA"
        if metadata["ReverseCompNeeded"] == "Yes":
            snp_call_uas = complement_base(snp_call)
        else:
            snp_call_uas = snp_call
        if snpid == "rs2402130":
            differ_length = len(seq) - 73
        else:
            differ_length = int(infile.iloc[j, 6])
        if snp_call not in expected_alleles and differ_length != 0:
            if snpid == "rs1821380":
                snp_call, allele_flag = snp_call_exception(seq, differ_length, metadata, snp_call)
                snp_call_uas = complement_base(snp_call)
            else:
                allele_flag = (
                    "Allele call does not match expected allele! Check for indels "
                    "(does not match expected sequence length)"
                )
        elif snp_call not in expected_alleles:
            allele_flag = "Allele call does not match expected allele!"
        elif differ_length != 0:
            allele_flag = "Check for indels (does not match expected sequence length)"
        else:
            allele_flag = None
        if allowed_snptype == "all":
            row_tmp = [
                name,
                analysis,
                analysis,
                snpid,
                seq,
                infile.loc[j, "Total_Reads"],
                snp_call,
                snp_call_uas,
                snp_type_dict[current_snp_type],
                allele_flag,
            ]
            all_rows.append(row_tmp)
        else:
            for snp in allowed_snptype:
                if snp in current_snp_type:
                    row_tmp = [
                        name,
                        analysis,
                        analysis,
                        snpid,
                        seq,
                        infile.loc[j, "Total_Reads"],
                        snp_call,
                        snp_call_uas,
                        snp_type_dict[current_snp_type],
                        allele_flag,
                    ]
                    all_rows.append(row_tmp)
    else:
        all_rows = None
    return all_rows


def snp_call_exception(seq, expected_size, metadata, base):
    """
    This function accounts for insertions and deletions in sequences to identify the correct base
    coordinate for the SNP. If the identified allele is still not one of the expected alleles, the
    sequence will be flagged appropriately.
    """
    new_size = len(seq) + expected_size
    new_base_call = seq[new_size]
    if new_base_call in metadata["Alleles"]:
        flag = (
            "Sequence length different than expected (check for indels); allele position adjusted"
        )
        return new_base_call, flag
    else:
        flag = (
            "Allele call does not match expected allele! Check for indels "
            "(does not match expected sequence length)"
        )
        return base, flag


def main(input, output, kit, uas, snptypes, nofilter):
    output = str(output)
    input = str(input)
    output_name = os.path.splitext(output)[0]
    if uas:
        results = uas_format(input, snptypes, nofilter)
        results.to_csv(output, index=False, sep="\t")
    else:
        results, results_combined = strait_razor_format(input, snptypes)
        results_combined.to_csv(output, index=False, sep="\t")
        results.to_csv(f"{output_name}_full_output.txt", index=False, sep="\t")


if __name__ == "__main__":
    main(
        snakemake.input,
        snakemake.output,
        kit=snakemake.params.kit,
        uas=snakemake.params.uas,
        snptypes=snakemake.params.types,
        nofilter=snakemake.params.nofilter,
    )
