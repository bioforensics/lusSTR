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

import lusSTR
import argparse
import glob
import openpyxl
import os
import pandas as pd
import sys


def uas_load(inpath, sexloci=False):
    """Format a UAS Sample Details Report (.xlsx).

    The `inpath` argument can refer to a report file or a directory of report files. Any files
    without the `.xlsx` file extension are ignored. The `sexloci` argument determines whether X and
    Y chromosome STRs are included along with autosomal STRs.
    """
    if os.path.isdir(inpath):
        auto_strs = pd.DataFrame()
        sex_strs = pd.DataFrame() if sexloci is True else None
        files = glob.glob(os.path.join(inpath, "[!~]*.xlsx"))
        for filename in sorted(files):
            if "Sample Details" not in filename:
                continue
            autodata, sexdata = uas_format(filename, sexloci)
            auto_strs = pd.concat([auto_strs, autodata])
            if sexloci is True:
                sex_strs = pd.concat([sex_strs, sexdata])
    else:
        auto_strs, sex_strs = uas_format(inpath, sexloci)
    return auto_strs, sex_strs


def parse_str_table_from_sheet(infile, sheet, exclude=None):
    file = openpyxl.load_workbook(infile)
    file_sheet = file[sheet]
    table = pd.DataFrame(file_sheet.values)
    offset = table[table.iloc[:, 0] == "Coverage Information"].index.tolist()[0]
    data = table.iloc[offset + 2 :]
    data.columns = table.iloc[offset + 1]
    if exclude is not None:
        data = data[~data.Locus.isin(exclude)]
    data = data[["Locus", "Reads", "Repeat Sequence"]]
    data["SampleID"] = table.iloc[2, 1]
    data["Project"] = table.iloc[3, 1]
    data["Analysis"] = table.iloc[4, 1]
    return data


def uas_format(infile, sexloci=False):
    auto_strs = parse_str_table_from_sheet(infile, sheet="Autosomal STRs")
    sex_strs = None
    if sexloci is True:
        y_strs = parse_str_table_from_sheet(infile, "Y STRs")
        x_strs = parse_str_table_from_sheet(infile, "X STRs")
        sex_strs = pd.concat([y_strs, x_strs], ignore_index=True)
    return auto_strs, sex_strs


def nonuas_load(inpath, software, sexloci=False):
    """Format a directory of STRait Razor/GeneMarker output files."""
    locus_list = [
        "Amelogenin",
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
        "PentaD",
        "PENTAD",
        "Penta D",
        "PentaE",
        "PENTAE",
        "Penta E",
        "TH01",
        "TPOX",
        "vWA",
        "VWA",
    ]
    sex_locus_list = [
        "Y-GATA-H4",
        "DYS643",
        "DYS635",
        "DYS612",
        "DYS576",
        "DYS570",
        "DYS549",
        "DYS533",
        "DYS522",
        "DYS505",
        "DYS481",
        "DYS460",
        "DYS448",
        "DYS439",
        "DYS438",
        "DYS437",
        "DYS392",
        "DYS391",
        "DYS390",
        "DYS389I",
        "DYS389II",
        "DYS385a-b",
        "DYS19",
        "DYF387S1",
        "DYS393",
        "DYS458",
        "DYS456",
        "HPRTB",
        "DXS8378",
        "DXS7423",
        "DXS7132",
        "DXS10135",
        "DXS10074",
        "DXS10103",
        "DYS385",
    ]
    auto_strs = pd.DataFrame()
    sex_strs = pd.DataFrame() if sexloci is True else None
    if os.path.isdir(inpath):
        analysisID = os.path.basename(inpath.rstrip(os.sep))
        suf = ".txt" if software == "straitrazor" else ".csv"
        files = glob.glob(os.path.join(inpath, f"[!~]*{suf}"))
        for filename in sorted(files):
            try:
                if software == "straitrazor":
                    table = strait_razor_table(filename, analysisID, sexloci)
                elif software == "genemarker":
                    table = genemarker_table(filename, analysisID, sexloci)
                else:
                    print("Incorrect analysis software specified!")
                    break
            except ValueError:
                print(
                    f"Error found with {filename}. Will bypass and continue. Please check file"
                    f" and rerun the command, if necessary."
                )
                continue
            auto_strs = pd.concat([auto_strs, table[table.Locus.isin(locus_list)]])
            if sexloci is True:
                sex_strs = pd.concat([sex_strs, table[table.Locus.isin(sex_locus_list)]])
    else:
        if software == "straitrazor":
            table = strait_razor_table(inpath, "NA", sexloci)
        elif software == "genemarker":
            table = genemarker_table(inpath, "NA", sexloci)
        else:
            print("Incorrect analysis software specified!")
        auto_strs = pd.concat([auto_strs, table[table.Locus.isin(locus_list)]])
        if sexloci is True:
            sex_strs = pd.concat([sex_strs, table[table.Locus.isin(sex_locus_list)]])
    return auto_strs, sex_strs


def strait_razor_table(filename, analysisID, sexloci=False):
    name = filename.replace(".txt", "").split(os.sep)[-1]
    table = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=["Locus_allele", "Length", "Sequence", "Forward_Reads", "Reverse_Reads"],
    )
    try:
        table[["Locus", "Allele"]] = table.Locus_allele.str.split(":", expand=True)
    except ValueError:
        table = []
    table["Total_Reads"] = table["Forward_Reads"] + table["Reverse_Reads"]
    table["SampleID"] = name
    table["Project"] = analysisID
    table["Analysis"] = analysisID
    table = table[["Locus", "Total_Reads", "Sequence", "SampleID", "Project", "Analysis"]]
    return table


def genemarker_table(filename, analysisID, sexloci=False):
    name = filename.replace(".csv", "").split(os.sep)[-1]
    table = pd.read_csv(
        filename, comment="#", usecols=["Marker", "Sequence Total Count", "Sequence"]
    )
    table["SampleID"] = name
    table["Project"] = analysisID
    table["Analysis"] = analysisID
    table["Sequence"] = table["Sequence"].str.upper()
    table = table.rename(columns={"Marker": "Locus", "Sequence Total Count": "Total_Reads"})
    return table[["Locus", "Total_Reads", "Sequence", "SampleID", "Project", "Analysis"]]


def main(input, outfile, a_software, sex=False):
    if a_software == "uas":
        results, sex_results = uas_load(str(input), sex)
    else:
        results, sex_results = nonuas_load(str(input), a_software, sex)
    if outfile is None:
        outfile = sys.stdout
    results.to_csv(str(outfile), index=False)
    if sex:
        name = os.path.splitext(str(outfile))[0]
        sex_results.to_csv(f"{name}_sexloci.csv", index=False)


if __name__ == "__main__":
    main(
        snakemake.input,
        snakemake.output,
        a_software=snakemake.params.a_software,
        sex=snakemake.params.sex,
    )
