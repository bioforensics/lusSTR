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
filter_sep = config["filter_sep"]


def get_sample_IDs(input, uas, output, software, separate):
    file_ext = ".xlsx" if uas is True else ".txt"
    if software == "efm" and separate is False:
        return os.path.basename(output)
    else:
        if uas is True:
            if os.path.isdir(input):
                files = glob.glob(os.path.join(input, f"[!~]*{file_ext}"))
            else:
                files = input
            ID_list = get_uas_ids(files)
        else:
            if os.path.isdir(input):
                files = glob.glob(os.path.join(input, f"[!~]*{file_ext}"))
            else:
                files = input
            files = [sub.replace(dir, "") for sub in files]
            ID_list = [sub.replace(file_ext, "") for sub in files]
        return ID_list


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
            filter_sep), prof_t=prof, data_t=data
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


rule annotate:
    input:
        rules.format.output
    output:
        expand("{name}.txt", name=output_name)
    params:
        uas=config["uas"],
        sex=config["sex"],
        combine=config["nocombine"],
        separate=config["separate"],
        kit=config["kit"]
    script:
        lusSTR.wrapper("annot")


rule filter:
    input:
        rules.annotate.output
    output:
        expand(
            "{outdir}/{samplename}_{prof_t}_{data_t}.csv", outdir=output_name,
            samplename=get_sample_IDs(input_name, config["uas"], output_name, software, 
            filter_sep), prof_t=prof, data_t=data
        )
    params:
        output_type=config["output_type"],
        profile_type=config["profile_type"],
        data_type=config["data_type"],
        output_dir=config["output"],
        info=config["info"],
        filter_sep=config["filter_sep"],
        filters=config["nofilters"],
    script:
        lusSTR.wrapper("filter")
    