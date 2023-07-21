import glob
import lusSTR
import openpyxl
import os
import pandas as pd
from pathlib import Path
import re


configfile: "config.yaml"
output_name = config["output"]
input_name = config["samp_input"]
software = config["output_type"]
prof = config["profile_type"]
data = config["data_type"]
separate = config["separate"]


def get_sample_IDs(input, uas, output, software, separate):
    convert_out = f"{output}.txt"
    format_out = f"{output}.csv"
    if software == "efm" and separate is False:
        ID_list = os.path.basename(output)
    elif os.path.exists(convert_out):
        ID_list = get_existing_IDs(convert_out, "\t")
    elif os.path.exists(format_out):
        ID_list = get_existing_IDs(format_out, ",")
    else:
        file_ext = ".xlsx" if uas is True else ".txt"
        if uas is True:
            if os.path.isdir(input):
                files = glob.glob(os.path.join(input, f"[!~]*{file_ext}"))
            else:
                files = input
            ID_list = get_uas_ids(files)
        else:
            if os.path.isdir(input):
                files = glob.glob(os.path.join(input, f"[!~]*{file_ext}"))
                files = [sub.replace(input, "") for sub in files]
                ID_list = [sub.replace(file_ext, "") for sub in files]
            else:
                files = os.path.basename(input)
                ID_list = files.replace(file_ext, "")
    return ID_list


def get_existing_IDs(infile, separator):
    data = pd.read_csv(infile, sep=separator)
    IDs = data["SampleID"].unique()
    return IDs


def get_uas_ids(files):
    samplelist = []
    if isinstance(files, list):
        for filename in sorted(files):
            if "Sample Details" not in filename:
                continue
            sampleID = parse_sample_details(filename)
            samplelist.append(sampleID)
    else:
        samplelist = parse_sample_details(files)
    return samplelist


def parse_sample_details(filename):
    file = openpyxl.load_workbook(filename)
    file_sheet = file["Autosomal STRs"]
    table = pd.DataFrame(file_sheet.values)
    sampleID = re.sub(" ", "_", table.iloc[2, 1])
    return sampleID


rule all:
    input:
        expand("{name}.csv", name=output_name),
        expand("{name}.txt", name=output_name),
        expand(
            "{outdir}/{samplename}_{prof_t}_{data_t}.csv", outdir=output_name,
            samplename=get_sample_IDs(input_name, config["uas"], output_name, software, 
            separate), prof_t=prof, data_t=data
        )


rule format:
    input:
       expand("{samp_input}", samp_input=input_name)
    output:
        expand("{name}.csv", name=output_name)
    params:
       uas=config["uas"],
       sex=config["sex"]
    script:
        lusSTR.wrapper("format")


rule convert:
    input:
        rules.format.output
    output:
        expand("{name}.txt", name=output_name)
    params:
        uas=config["uas"],
        sex=config["sex"],
        nocombine=config["nocombine"],
        kit=config["kit"]
    script:
        lusSTR.wrapper("convert")


rule filter:
    input:
        rules.convert.output
    output:
        expand(
            "{outdir}/{samplename}_{prof_t}_{data_t}.csv", outdir=output_name,
            samplename=get_sample_IDs(input_name, config["uas"], output_name, software, 
            separate), prof_t=prof, data_t=data
        )
    params:
        output_type=config["output_type"],
        profile_type=config["profile_type"],
        data_type=config["data_type"],
        output_dir=config["output"],
        info=config["info"],
        separate=config["separate"],
        filters=config["nofilters"],
        strand=config["strand"]
    script:
        lusSTR.wrapper("filter")
    