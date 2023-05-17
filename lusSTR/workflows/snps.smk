import glob
import lusSTR
import openpyxl
import os
import pandas as pd
from pathlib import Path
import re


## placeholder until I update for snps

configfile: "snp_config.yaml"
output_name = config["output"]
input_name = config["samp_input"]
separate = config["separate"]
refs = config["references"]


def get_sample_IDs(input, uas, output, separate, refs):
    file_ext = ".xlsx" if uas is True else ".txt"
    if separate is False:
        return os.path.basename(output)
    else:
        if uas is True:
            ID_list = get_uas_ids(files, refs)
        else:
            ID_list = get_straitrazor_ids(input, refs)
        print(ID_list)
        return ID_list


def get_uas_ids(input, refs):
    samplelist = []
    if os.path.isdir(input):
        files = glob.glob(os.path.join(input, f"[!~]*{file_ext}"))
        for filename in sorted(files):
            if "Sample Details" in filename:
                sampleID = parse_uas(filename, "Autosomal STRs", refs)
            elif "Phenotype" in filename:
                sampleID = parse_uas(filename, "SNP Data", refs)
            else:
                continue
            sampleID = parse_sample_details(filename)
            samplelist.append(sampleID)
    else:
        samplelist = parse_sample_details(input)
    return samplelist


def parse_sample_details(filename, sheet, refs):
    file = openpyxl.load_workbook(filename)
    file_sheet = file[sheet]
    table = pd.DataFrame(file_sheet.values)
    sampleID = re.sub(" ", "_", table.iloc[2, 1])
    if sampleID in refs:
        total_ID = f"{sampleID}_reference"
    else:
        total_ID = f"{sampleID}_evidence"
    return total_ID


def get_straitrazor_ids(input, refs):
    if os.path.isdir(input):
        samplelist = []
        files = glob.glob(os.path.join(input, f"[!~]*{file_ext}"))
        for filename in sorted(files):
            total_ID = extract_ID(filename, refs)
            samplelist.append(total_ID)   
    else:
        samplelist = extract_ID(input, refs)
    return samplelist


def extract_ID(filename, refs):
    drop_path = filename.replace(dir, "")
    drop_ext = drop_path.replace(".txt", "")
    if drop_ext in refs:
        total_ID = f"{drop_ext}_reference"
    else:
        total_ID = f"{drop_ext}_evidence"
    return total_ID

rule all:
    input:
        expand("{name}.txt", name=output_name),
        expand(
            "{outdir}/{samplename}.csv", outdir=output_name,
            samplename=get_sample_IDs(input_name, config["uas"], output_name,
            separate, refs)
        )


rule convert:
    input:
       expand("{samp_input}", samp_input=input_name)
    output:
        expand("{name}.txt", name=output_name)
    params:
       uas=config["uas"],
       kit=config["kit"],
       types=config["types"],
       nofilter=config["nofilter"]
    script:
        lusSTR.wrapper("snps_convert")


rule format:
    input:
        rules.convert.output
    output:
        expand(
                "{outdir}/{samplename}.csv", outdir=output_name,
                samplename=get_sample_IDs(input_name, config["uas"], output_name,
                separate, refs)
            )
    params:
        uas=config["uas"],
        separate=config["separate"],
        kit=config["kit"]
    script:
        lusSTR.wrapper("snps_format")
