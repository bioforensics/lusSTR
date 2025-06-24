from datetime import datetime
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
custom = config["custom_ranges"]


def get_sample_IDs(input, a_software, output, software, separate):
    convert_out = f"{output}.txt"
    format_out = f"{output}.csv"
    if os.path.exists(convert_out):
        ID_list = get_existing_IDs(convert_out, "\t")
    elif os.path.exists(format_out):
        ID_list = get_existing_IDs(format_out, ",")
    else:
        file_ext = ".xlsx" if a_software == "uas" else ".txt"
        if a_software == "uas":
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


def create_log(log):
    now = datetime.now()
    dt = now.strftime("%m%d%Y_%H_%M_%S")
    shell("mkdir -p logs/{dt}/input/")
    shell("cp '{log}' logs/{dt}/")
    if os.path.isdir(input_name):
        shell("cp '{input_name}'/*.* logs/{dt}/input/")
    else:
        shell("cp '{input_name}' logs/{dt}/input/")
    shell("cp config.yaml logs/{dt}/")


def get_output():
    if custom:
        outname = expand("{name}_custom_range.txt", name=output_name)
    else:
        outname = expand("{name}.txt", name=output_name)
    return outname


def get_markerplot_name(output, custom):
    if custom:
        return f"{output}_custom_range"
    else:
        return output


rule all:
    input:
        expand("{name}.csv", name=output_name),
        expand("{name}.txt", name=output_name),
        expand(
            "{outdir}/{samplename}_{prof_t}_{data_t}.csv", outdir=output_name,
            samplename=get_sample_IDs(input_name, config["analysis_software"], output_name, software, 
            separate), prof_t=prof, data_t=data
        )


rule format:
    input:
       expand("{samp_input}", samp_input=input_name)
    output:
        expand("{name}.csv", name=output_name)
    params:
       a_software=config["analysis_software"],
       sex=config["sex"]
    script:
        lusSTR.wrapper("format")


rule convert:
    input:
        rules.format.output
    output:
        get_output()
    params:
        a_software=config["analysis_software"],
        sex=config["sex"],
        nocombine=config["nocombine"],
        kit=config["kit"],
        custom=config["custom_ranges"]
    script:
        lusSTR.wrapper("convert") 


rule filter:
    input:
        rules.convert.output
    output:
        expand(
            "MarkerPlots/{output_name}_{samplename}_marker_plots.pdf", output_name=get_markerplot_name(config["output"], config["custom_ranges"]), 
            samplename=get_sample_IDs(input_name, config["analysis_software"], output_name, software, 
            separate)
        )
    params:
        output_type=config["output_type"],
        profile_type=config["profile_type"],
        data_type=config["data_type"],
        output_dir=config["output"],
        info=config["info"],
        separate=config["separate"],
        filters=config["nofilters"],
        strand=config["strand"],
        custom=config["custom_ranges"],
        sex=config["sex"],
        kit=config["kit"]
    script:
        lusSTR.wrapper("filter")


onsuccess:
    create_log(log)

onerror:
    create_log(log)
    