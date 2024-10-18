# -------------------------------------------------------------------------------------------------
# Copyright (c) 2024, DHS.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

from . import select
from .util import generate_config_file, validate_prefix
from datetime import datetime
from lusSTR.wrappers.filter import get_at, EFM_output, marker_plots, STRmix_output
import math
import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
import streamlit as st
import yaml
import subprocess
import os


def str_workflow_display():
    st.title("STR Workflow")
    st.info(
        "Please Select STR Settings Below for lusSTR. For information regarding the settings, see the How to Use tab."
    )
    str_input()
    str_output()
    str_general_settings()
    str_filter_settings()
    str_footer()


def str_input():
    st.subheader("Specify Worfklow Inputs")
    st.info(
        "Indicate whether you are providing an individual file or a folder containing multiple files"
    )
    if "samp_input_str" not in st.session_state:
        st.session_state.samp_input_str = None
    input_option = st.radio(
        "Select Input Option:", ("Individual file", "Folder with multiple files")
    )
    if input_option == "Folder with Multiple Files":
        clicked = st.button("Select a folder")
        if clicked:
            st.session_state.samp_input_str = select.folder()
    else:
        clicked_file = st.button("Select a file")
        if clicked_file:
            st.session_state.samp_input_str = select.file()
    if st.session_state.samp_input_str:
        st.text_input("Location of your input file(s):", st.session_state.samp_input_str)


def str_output():
    st.subheader("Specify Working Directory for Workflow Outputs")
    columns = st.columns(5)
    if "wd_dirname_str" not in st.session_state:
        st.session_state.wd_dirname_str = None
    clicked_wd = columns[0].button("Select Output Folder")
    if clicked_wd:
        st.session_state.wd_dirname_str = select.folder()
    if st.session_state.wd_dirname_str:
        st.text_input("Working directory:", st.session_state.wd_dirname_str)


def str_general_settings():
    st.subheader("General Settings")
    columns = st.columns(5)
    if "analysis_software" not in st.session_state:
        st.session_state.analysis_software = None
    selected_software = columns[0].selectbox(
        "Analysis Software",
        options=["UAS", "STRait Razor v3", "GeneMarker HTS"],
        help="Indicate the analysis software used prior to lusSTR.",
    )
    software = {
        "UAS": "uas",
        "STRait Razor v3": "straitrazor",
        "GeneMarker HTS": "genemarker",
    }
    st.session_state.analysis_software = software[selected_software]
    if "custom_ranges" not in st.session_state:
        st.session_state.custom_ranges = None
    st.session_state.custom_ranges = st.checkbox(
        "Use Custom Sequence Ranges",
        help="Check the box to use the specified custom sequence ranges as defined in the `str_markers.json` file.",
    )
    if "sex" not in st.session_state:
        st.session_state.sex = None
    st.session_state.sex = st.checkbox(
        "Include X- and Y-STRs",
        help="Check the box to include X- and Y-STRs, otherwise leave unchecked.",
    )
    if "kit" not in st.session_state:
        st.session_state.kit = None
    selected_kit = columns[1].selectbox(
        "Library Preparation Kit",
        options=["ForenSeq Signature Prep", "PowerSeq 46GY"],
        help="Specify the library preparation kit used to generate the sequences.",
    )
    kits = {"ForenSeq Signature Prep": "forenseq", "PowerSeq 46GY": "powerseq"}
    st.session_state.kit = kits[selected_kit]
    if "output" not in st.session_state:
        st.session_state.output = None
    st.session_state.output = columns[2].text_input(
        "Output File Name",
        "lusstr_output",
        help="Please specify a name for the created files. It can only contain alphanumeric characters, underscores and hyphens. No spaces allowed.",
    )
    if "nocombine" not in st.session_state:
        st.session_state.nocombine = None
    st.session_state.nocombine = st.checkbox(
        "Do Not Combine Identical Sequences",
        help="If using STRait Razor data, by default, identical sequences (after removing flanking sequences) are combined and reads are summed. Checking this will not combine identical sequences.",
    )


def str_filter_settings():
    st.subheader("Filter Settings")
    columns = st.columns(5)
    if "output_type" not in st.session_state:
        st.session_state.output_type = None
    st.session_state.output_type = {
        "STRmix": "strmix",
        "EuroForMix": "efm",
        "MPSproto": "mpsproto",
    }[
        columns[0].selectbox(
            "Probabilistic Genotyping Software",
            options=["STRmix", "EuroForMix", "MPSproto"],
            help="Select which probabilistic genotyping software files to create",
        )
    ]
    if "profile_type" not in st.session_state:
        st.session_state.profile_type = None
    st.session_state.profile_type = {"Evidence": "evidence", "Reference": "reference"}[
        columns[1].selectbox(
            "Profile Type",
            options=["Evidence", "Reference"],
            help="Select the file type (format) to create for the probabilistic genotyping software.",
        )
    ]
    if "data_type" not in st.session_state:
        st.session_state.data_type = None
    st.session_state.data_type = {
        "Sequence": "ngs",
        "CE allele": "ce",
        "LUS+ allele": "lusplus",
    }[
        columns[2].selectbox(
            "Data Type",
            options=["Sequence", "CE allele", "LUS+ allele"],
            help="Select the allele type used to determine sequence type (belowAT, stutter or typed) and used in the final output file.",
        )
    ]
    if "info" not in st.session_state:
        st.session_state.info = None
    st.session_state.info = st.checkbox(
        "Create Allele Information File",
        value=True,
        help="Create file containing information about each sequence, including sequence type (belowAT, stutter or typed), stuttering sequence information and metrics involving stutter and noise.",
    )
    if "separate" not in st.session_state:
        st.session_state.separate = None
    st.session_state.separate = st.checkbox(
        "Create Separate Files for Samples",
        help="If checked, will create individual files for samples; If unchecked, will create one file with all samples.",
    )
    if "nofilters" not in st.session_state:
        st.session_state.nofilters = None
    st.session_state.nofilters = st.checkbox(
        "Skip All Filtering Steps",
        help="Filtering will not be performed but will still create EFM/MPSproto/STRmix output files containing all sequences.",
    )
    if "strand" not in st.session_state:
        st.session_state.strand = None
    st.session_state.strand = {"UAS Orientation": "uas", "Forward Strand": "forward"}[
        columns[3].selectbox(
            "Strand Orientation",
            options=["Forward Strand", "UAS Orientation"],
            help="Indicates the strand orientation in which to report the sequence in the final output table as some markers are reported in the UAS on the reverse strand. Selecting the UAS Orientation will report those markers on the reverse strand while the remaining will be reported on the forward strand. Selecting the Forward Strand will report all markers on the forward strand orientation. This applies to STRmix NGS only.",
        )
    ]


def str_launch_workflow():
    if st.button("Run lusSTR"):
        if (
            st.session_state.analysis_software
            and st.session_state.samp_input_str
            and st.session_state.output
            and st.session_state.wd_dirname_str
        ):
            if not validate_prefix(st.session_state.output):
                st.warning(
                    "Please enter a valid output prefix. Only alphanumeric characters, underscore, and hyphen are allowed."
                )
                st.stop()
            with st.spinner("Processing your data..."):
                config_data = {
                    "analysis_software": st.session_state.analysis_software,
                    "custom_ranges": st.session_state.custom_ranges,
                    "sex": st.session_state.sex,
                    "samp_input_str": st.session_state.samp_input_str,
                    "output": st.session_state.output,
                    "kit": st.session_state.kit,
                    "nocombine": st.session_state.nocombine,
                    "output_type": st.session_state.output_type,
                    "profile_type": st.session_state.profile_type,
                    "data_type": st.session_state.data_type,
                    "info": st.session_state.info,
                    "separate": st.session_state.separate,
                    "nofilters": st.session_state.nofilters,
                    "strand": st.session_state.strand,
                }
                generate_config_file(config_data, st.session_state.wd_dirname_str, "STR")
                command = ["lusstr", "strs", "all"]
                if wd_dirname_str:
                    command.extend(["-w", st.session_state.wd_dirname_str + "/"])
                try:
                    subprocess.run(command, check=True)
                    st.success(
                        "Config file generated and lusSTR executed successfully! Output files have been saved to your designated directory and labeled with your specified prefix."
                    )
                except subprocess.CalledProcessError as e:
                    st.error(f"Error: {e}")
                    st.info(
                        "Please make sure to check the 'How to Use' tab for common error resolutions."
                    )
        else:
            st.warning(
                "Please make sure to fill out all required fields (Analysis Software, Input Directory or File, Prefix for Output, and Specification of Working Directory) before submitting."
            )


def str_footer():
    st.write("---")
    st.write(
        "After running lusSTR, or if lusSTR has been run previously, the user may view and edit the individual STR marker plots and data."
    )
    st.write(
        "If lusSTR has been previously run, only the above ```Output Folder``` containing the run files needs to be specified. Other settings will be automatically loaded from the config.yaml file within the specified folder."
    )
    if "interactive" not in st.session_state:
        st.session_state.interactive = None
    if st.button("See Individual Marker Plots & Data") or st.session_state.interactive:
        st.session_state.interactive = True
        create_settings()
        file = f"{st.session_state.wd_dirname_str}/{st.session_state.output}/{st.session_state.output}"
        if st.session_state.custom_ranges:
            file += "_custom_range"
        try:
            sequence_info = pd.read_csv(f"{file}_sequence_info.csv")
            interactive_setup(sequence_info, file)
        except FileNotFoundError:
            print(f"{file}_sequence_info.csv not found. Please check output folder specification.")


def df_on_change(locus):
    state = st.session_state[f"{locus}_edited"]
    for index, updates in state["edited_rows"].items():
        st.session_state[locus].loc[st.session_state[locus].index == index, "edited"] = True
        for key, value in updates.items():
            st.session_state[locus].loc[st.session_state[locus].index == index, key] = value


def interactive_plots_allmarkers(sample_df, flagged_df):
    cols = st.columns(4)
    max_reads = max(sample_df["Reads"])
    n = 100 if max_reads > 1000 else 10
    max_yvalue = int(math.ceil(max_reads / n)) * n
    increase_value = int(math.ceil((max_yvalue / 5)) / n) * n
    n = 0
    for marker in sample_df["Locus"].unique():
        col = cols[n]
        container = col.container(border=True)
        sample_locus = sample_df["SampleID"].unique() + "_" + marker
        marker_df = sample_df[sample_df["Locus"] == marker].sort_values(by="CE_Allele")
        if sample_locus in flagged_df["key"].values:
            marker = f"⚠️{marker}⚠️"
        plot = interactive_plots(marker_df, marker, max_yvalue, increase_value, all=True)
        container.plotly_chart(plot, use_container_width=True)
        if n == 3:
            n = 0
        else:
            n += 1


def interactive_plots(df, locus, ymax, increase, all=False):
    if "⚠️" in locus:
        locus_at = locus.replace("⚠️", "")
    else:
        locus_at = locus
    at = get_at(df, locus_at)
    for i, row in df.iterrows():
        if "stutter" in df.loc[i, "allele_type"]:
            df.loc[i, "Label"] = "Stutter"
        else:
            df.loc[i, "Label"] = df.loc[i, "allele_type"]
    min_x = round(min(df["CE_Allele"]) - 1)
    max_x = round(max(df["CE_Allele"]) + 1)
    plot = px.bar(
        df,
        x="CE_Allele",
        y="Reads",
        color="Label",
        color_discrete_map={
            "Typed": "green",
            "BelowAT": "red",
            "Stutter": "blue",
            "Deleted": "purple",
        },
        title=locus,
    )
    plot.add_hline(y=at, line_width=3, line_dash="dot", line_color="gray")
    plot.add_annotation(text="AT", x=min_x + 0.1, y=at, showarrow=False, yshift=10)
    plot.update_layout(
        xaxis=dict(range=[min_x, max_x], tickmode="array", tickvals=np.arange(min_x, max_x, 1))
    )
    if all:
        plot.update_layout(
            yaxis=dict(range=[0, ymax], tickmode="array", tickvals=np.arange(0, ymax, increase))
        )
    return plot


def remake_final_files(full_df, outpath):
    if st.session_state.custom_ranges:
        seq_col = "Custom_Range_Sequence"
        brack_col = "Custom_Bracketed_Notation"
    else:
        seq_col = (
            "UAS_Output_Sequence"
            if st.session_state.strand == "uas"
            else "Forward_Strand_Sequence"
        )
        brack_col = (
            "UAS_Output_Bracketed_Notation"
            if st.session_state.strand == "uas"
            else "Forward_Strand_Bracketed_Notation"
        )
    if st.session_state.nofilters:
        full_df["allele_type"] = "Typed"
    if st.session_state.output_type == "efm" or st.session_state.output_type == "mpsproto":
        EFM_output(
            full_df,
            outpath,
            st.session_state.profile_type,
            st.session_state.data_type,
            brack_col,
            st.session_state.sex,
            st.session_state.separate,
        )
    else:
        STRmix_output(
            full_df, outpath, st.session_state.profile_type, st.session_state.data_type, seq_col
        )


def interactive_setup(df1, file):
    col1, col2, col3, col4, col5 = st.columns(5)
    sample = col1.selectbox("Select Sample:", options=df1["SampleID"].unique())
    sample_df = df1[df1["SampleID"] == sample].reset_index(drop=True)
    locus_list = pd.concat([pd.Series("All Markers"), sample_df["Locus"].drop_duplicates()])
    if os.path.isfile(f"{file}_Flagged_Loci.csv"):
        flags = pd.read_csv(f"{file}_Flagged_Loci.csv")
    else:
        flags = pd.DataFrame(columns=["key", "SampleID", "Locus"])
    flags["key"] = flags["SampleID"] + "_" + flags["Locus"]
    flags_sample = flags[flags["SampleID"] == sample].reset_index(drop=True)
    for flagged_locus in flags_sample["Locus"].unique():
        locus_list = locus_list.str.replace(flagged_locus, f"⚠️{flagged_locus}⚠️")
    locus = col2.selectbox("Select Marker:", options=locus_list)
    if "⚠️" in locus:
        locus = locus.replace("⚠️", "")
    if locus == "All Markers":
        if not flags_sample.empty:
            st.write(
                "⚠️ indicates potential problems with the marker. Examine the individual marker plots for more information."
            )
        interactive_plots_allmarkers(sample_df, flags)
    else:
        locus_key = f"{sample}_{locus}"
        if locus_key not in st.session_state:
            st.session_state[locus_key] = sample_df[sample_df["Locus"] == locus].reset_index(
                drop=True
            )
        Type = [
            "Deleted",
            "Typed",
            "-1_stutter",
            "-2_stutter",
            "BelowAT",
            "-1_stutter/+1_stutter",
            "+1_stutter",
        ]
        plot = interactive_plots(st.session_state[locus_key], locus, None, None)
        st.plotly_chart(plot, use_container_width=True)
        col1, col2, col3 = st.columns(3)
        if locus_key in flags["key"].values:
            locus_flags = flags[flags["key"] == locus_key]
            for flag in locus_flags["Flags"].unique():
                col2.write(f"⚠️ Potential issue: {flag} identified!")
        st.data_editor(
            data=st.session_state[locus_key],
            disabled=(
                "SampleID",
                "Locus",
                "UAS_Output_Sequence",
                "CE_Allele",
                "UAS_Output_Bracketed_Notation",
                "Custom_Range_Sequence",
                "Custom_Bracketed_Notation",
                "Reads",
                "parent_allele1",
                "parent_allele2",
                "allele1_ref_reads",
                "allele2_ref_reads",
                "perc_noise",
                "perc_stutter",
            ),
            column_config={
                "allele_type": st.column_config.SelectboxColumn("allele_type", options=Type)
            },
            hide_index=True,
            key=f"{locus_key}_edited",
            on_change=df_on_change,
            args=(locus_key,),
        )
    if st.button("Save Edits"):
        ph = st.empty()
        with ph.container():
            st.write("Saving Changes - May take a minute or two.")
        combined_df = pd.DataFrame()
        for sample in df1["SampleID"].unique():
            sample_df = df1[df1["SampleID"] == sample].reset_index(drop=True)
            for locus in sample_df["Locus"].unique():
                locus_key = f"{sample}_{locus}"
                try:
                    combined_df = pd.concat([combined_df, st.session_state[locus_key]])
                except KeyError:
                    combined_df = pd.concat(
                        [
                            combined_df,
                            sample_df[sample_df["Locus"] == locus].reset_index(drop=True),
                        ]
                    )
        now = datetime.now()
        dt = now.strftime("%m%d%Y_%H_%M_%S")
        del combined_df["Label"]
        Path(f"{st.session_state.wd_dirname_str}/{st.session_state.output}/edited_{dt}").mkdir(
            parents=True, exist_ok=True
        )
        outpath = f"{st.session_state.wd_dirname_str}/{st.session_state.output}/edited_{dt}/"
        marker_plots(combined_df, f"{st.session_state.output}_edited_{dt}", sex=False, wd=outpath)
        combined_df.to_csv(
            f"{st.session_state.wd_dirname_str}/{st.session_state.output}/edited_{dt}/"
            f"{st.session_state.output}_sequence_info_edited_{dt}.csv",
            index=False,
        )
        new_text = (
            f"Changes saved to {st.session_state.wd_dirname_str}/{st.session_state.output}"
            f"/edited_{dt}/{st.session_state.output}_sequence_info_edited_{dt}.csv"
            f"New {st.session_state.output_type} files created in {st.session_state.wd_dirname_str}"
            f"/{st.session_state.output}/edited_{dt}/ folder"
        )
        remake_final_files(combined_df, outpath)
        ph.empty()
        with ph.container():
            st.write(
                f"New files and marker plots with edits saved to {st.session_state.wd_dirname_str}/"
                f"{st.session_state.output}/edited_{dt}/"
            )


def create_settings():
    if os.path.isfile(f"{st.session_state.wd_dirname_str}/config.yaml"):
        st.write(f"Loading settings from {st.session_state.wd_dirname_str}/config.yaml")
        with open(f"{st.session_state.wd_dirname_str}/config.yaml", "r") as file:
            config_settings = yaml.safe_load(file)
        st.session_state.output = config_settings["output"]
        st.session_state.custom_ranges = config_settings["custom_ranges"]
        st.session_state.profile_type = config_settings["profile_type"]
        st.session_state.data_type = config_settings["data_type"]
        st.session_state.sex = config_settings["sex"]
        st.session_state.separate = config_settings["separate"]
        st.session_state.strand = config_settings["strand"]
        st.session_state.output_type = config_settings["output_type"]
