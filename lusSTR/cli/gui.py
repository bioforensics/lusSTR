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
#################################################################
#              Importing Necessary Packages                     #
#################################################################

from datetime import datetime
import json
import importlib.resources
from lusSTR.wrappers.filter import get_at, EFM_output, marker_plots, make_plot, STRmix_output
import math
import numpy as np
import pandas as pd
from pathlib import Path
import plotly.express as px
import streamlit as st
from streamlit_option_menu import option_menu
import yaml
import subprocess
import os
import re


#################################################################
#                        Functions                              #
#################################################################


def get_filter_metadata_file():
    return importlib.resources.files("lusSTR") / "data/filters.json"


with open(get_filter_metadata_file(), "r") as fh:
    filter_marker_data = json.load(fh)


# ------------ Function to Generate config.yaml File ---------- #


def generate_config_file(config_data, working_directory, workflow_type):
    if workflow_type == "STR":
        config_filename = "config.yaml"
    elif workflow_type == "SNP":
        config_filename = "snp_config.yaml"
    else:
        raise ValueError("Invalid workflow type. Please specify either 'STR' or 'SNP'.")

    config_path = os.path.join(working_directory, config_filename)
    with open(config_path, "w") as file:
        yaml.dump(config_data, file)


# ------------ Function for folder selection ------------------ #


def folder_picker_dialog():
    script_path = importlib.resources.files("lusSTR") / "scripts" / "folder_selector.py"
    result = subprocess.run(["python", script_path], capture_output=True, text=True)
    if result.returncode == 0:
        folder_data = json.loads(result.stdout)
        folder_path = folder_data.get("folder_path")
        if folder_path:
            return folder_path
        else:
            st.error("No folder selected")
    else:
        st.error("Error selecting folder")


# ------- Function for individual file selection -------------- #


def file_picker_dialog():
    script_path = importlib.resources.files("lusSTR") / "scripts" / "file_selector.py"
    result = subprocess.run(["python", script_path], capture_output=True, text=True)
    if result.returncode == 0:
        file_data = json.loads(result.stdout)
        file_path = file_data.get("file_path")
        if file_path:
            return file_path
        else:
            st.error("No folder selected")
    else:
        st.error("Error selecting folder")


# ---- Function to validate prefix for output folder ---------- #


def validate_prefix(prefix):
    if re.match(
        r"^[A-Za-z0-9_-]+$", prefix
    ):  # Allow alphanumeric characters, underscore, and hyphen
        return True
    else:
        return False


#################################################################
#              Front-End Logic For Navigation Bar               #
#################################################################


def main():

    # Page Layout (Theme and Fonts have been established in .streamlit/config.toml)
    st.set_page_config(layout="wide", initial_sidebar_state="collapsed")

    # Creating Navigation Bar

    selected = option_menu(
        menu_title=None,
        options=["Home", "STRs", "SNPs", "How to Use", "Contact"],
        icons=["house", "gear", "gear-fill", "book", "envelope"],
        menu_icon="cast",
        default_index=0,
        orientation="horizontal",
    )

    if selected == "Home":
        show_home_page()

    elif selected == "STRs":
        show_STR_page()

    elif selected == "SNPs":
        show_SNP_page()

    elif selected == "How to Use":
        show_how_to_use_page()

    elif selected == "Contact":
        show_contact_page()


#####################################################################
#                     lusSTR Home Page                              #
#####################################################################


def show_home_page():

    image_path = importlib.resources.files("lusSTR") / "cli" / "logo.png"

    # CSS to hide full-screen button
    hide_img_fs = """
    <style>
    button[title="View fullscreen"]{
    visibility: hidden;}
    </style>
    """

    # Define column layout for centering image
    left_co, cent_co, last_co = st.columns([2.5, 8, 2.5])
    with cent_co:
        st.image(str(image_path), use_column_width="auto")

    # Apply CSS to hide full-screen button
    st.markdown(hide_img_fs, unsafe_allow_html=True)

    # -- Welcome Message Stuff

    st.markdown(
        """
        lusSTR is an end-to-end workflow for processing human forensic data (STRs and SNPs) 
        derived from Next Generation Sequencing (NGS) data for use in probabilistic genotyping 
        software. For more information on lusSTR, visit our 
        [GitHub page](https://github.com/bioforensics/lusSTR).
        """,
        unsafe_allow_html=True,
    )

    st.info("Please Select One of the Tabs Above to Get Started on Processing Your Data!")


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
        if df.loc[i, "allele_type"] == "Typed":
            df.loc[i, "Label"] = "Typed"
        elif df.loc[i, "allele_type"] == "BelowAT" or df.loc[i, "allele_type"] == "NotTyped":
            df.loc[i, "Label"] = "BelowAT"
        else:
            df.loc[i, "Label"] = "Stutter"
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
        },
        title=locus,
    )
    plot.add_hline(y=at, line_width=3, line_dash="dot", line_color="gray")
    plot.add_annotation(text=f"AT", x=min_x + 0.1, y=at, showarrow=False, yshift=10)
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
                f"⚠️ indicates potential problems with the marker. Examine the individual marker "
                f"plots for more information."
            )
        interactive_plots_allmarkers(sample_df, flags)
    else:
        locus_key = f"{sample}_{locus}"
        if locus_key not in st.session_state:
            st.session_state[locus_key] = sample_df[sample_df["Locus"] == locus].reset_index(drop=True)
        Type = [
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
        Path(f"{st.session_state.wd_dirname}/{st.session_state.output}/edited_{dt}").mkdir(
            parents=True, exist_ok=True
        )
        outpath = f"{st.session_state.wd_dirname}/{st.session_state.output}/edited_{dt}/"
        marker_plots(combined_df, f"{st.session_state.output}_edited_{dt}", sex=False, wd=outpath)
        combined_df.to_csv(
            f"{st.session_state.wd_dirname}/{st.session_state.output}/edited_{dt}/"
            f"{st.session_state.output}_sequence_info_edited_{dt}.csv",
            index=False,
        )
        new_text = (
            f"Changes saved to {st.session_state.wd_dirname}/{st.session_state.output}"
            f"/edited_{dt}/{st.session_state.output}_sequence_info_edited_{dt}.csv"
            f"New {st.session_state.output_type} files created in {st.session_state.wd_dirname}"
            f"/{st.session_state.output}/edited_{dt}/ folder"
        )
        remake_final_files(combined_df, outpath)
        ph.empty()
        with ph.container():
            st.write(
                f"New files and marker plots with edits saved to {st.session_state.wd_dirname}/"
                f"{st.session_state.output}/edited_{dt}/"
            )


def create_settings():
    if os.path.isfile(f"{st.session_state.wd_dirname}/config.yaml"):
        st.write(f"Loading settings from {st.session_state.wd_dirname}/config.yaml")
        with open(f"{st.session_state.wd_dirname}/config.yaml", "r") as file:
            config_settings = yaml.safe_load(file)
        st.session_state.output = config_settings["output"]
        st.session_state.custom_ranges = config_settings["custom_ranges"]
        st.session_state.profile_type = config_settings["profile_type"]
        st.session_state.data_type = config_settings["data_type"]
        st.session_state.sex = config_settings["sex"]
        st.session_state.separate = config_settings["separate"]
        st.session_state.strand = config_settings["strand"]
        st.session_state.output_type = config_settings["output_type"]


#####################################################################
#                        STR WORKFLOW                               #
#####################################################################

#####################################################################
# Specify STR Settings Which Will Be Used to Generate Config File   #
#####################################################################


def show_STR_page():

    st.title("STR Workflow")
    st.info(
        "Please Select STR Settings Below for lusSTR! For Information Regarding the "
        "Settings, See the How to Use Tab."
    )

    # Input File Specification
    st.subheader("Input Files Selection")

    # Ask user if submitting a directory or individual file
    st.info(
        "Please Indicate If You Are Providing An Individual Input File or a Folder Containing "
        "Multiple Input Files"
    )
    input_option = st.radio(
        "Select Input Option:", ("Individual File", "Folder with Multiple Files")
    )

    # Initialize session state if not already initialized
    if "samp_input" not in st.session_state:
        st.session_state.samp_input = None

    # Logic for Path Picker based on user's input option

    if input_option == "Folder with Multiple Files":
        clicked = st.button("Select a Folder")
        if clicked:
            dirname = folder_picker_dialog()
            st.session_state.samp_input = dirname

    else:
        clicked_file = st.button("Select a File")
        if clicked_file:
            filename = file_picker_dialog()
            st.session_state.samp_input = filename

    # Display The Selected Path
    if st.session_state.samp_input:
        st.text_input("Location Of Your Input File(s):", st.session_state.samp_input)

    # Store the Selected Path to Reference in Config
    samp_input = st.session_state.samp_input

    #####################################################################
    #     STR: Specify Working Directory                                #
    #####################################################################

    st.subheader("Output Folder Selection")

    col1, col2, col3, col4, col5 = st.columns(5)

    # Initialize session state if not already initialized
    if "wd_dirname" not in st.session_state:
        st.session_state.wd_dirname = None

    clicked_wd = col1.button("Select An Output Folder")
    if clicked_wd:
        wd = folder_picker_dialog()
        st.session_state.wd_dirname = wd

    # Display selected path
    if st.session_state.wd_dirname:
        st.text_input("Your Specified Output Folder:", st.session_state.wd_dirname)

    # Store Selected Path to Reference in Config
    wd_dirname = st.session_state.wd_dirname

    #####################################################################
    #      STR: General Software Settings to Generate Config File       #
    #####################################################################

    st.subheader("General Settings")

    col1, col2, col3, col4, col5 = st.columns(5)

    if "analysis_software" not in st.session_state:
        st.session_state.analysis_software = None

    st.session_state.analysis_software = {
        "UAS": "uas",
        "STRait Razor v3": "straitrazor",
        "GeneMarker HTS": "genemarker",
    }[
        col1.selectbox(
            "Analysis Software",
            options=["UAS", "STRait Razor v3", "GeneMarker HTS"],
            help="Indicate the analysis software used prior to lusSTR.",
        )
    ]

    if "custom_ranges" not in st.session_state:
        st.session_state.custom_ranges = None

    st.session_state.custom_ranges = st.checkbox(
        "Use Custom Sequence Ranges",
        help="Check the box to use the specified custom sequence ranges as defined in the "
        "str_markers.json file.",
    )

    if "sex" not in st.session_state:
        st.session_state.sex = None

    st.session_state.sex = st.checkbox(
        "Include X- and Y-STRs",
        help="Check the box to include X- and Y-STRs, otherwise leave unchecked.",
    )

    if "kit" not in st.session_state:
        st.session_state.kit = None

    st.session_state.kit = {"ForenSeq Signature Prep": "forenseq", "PowerSeq 46GY": "powerseq"}[
        col2.selectbox(
            "Library Preparation Kit",
            options=["ForenSeq Signature Prep", "PowerSeq 46GY"],
            help="Specify the library preparation kit used to generate the sequences.",
        )
    ]

    if "output" not in st.session_state:
        st.session_state.output = None

    st.session_state.output = col3.text_input(
        "Output File Name",
        "lusstr_output",
        help="Please specify a name for the created files. It can only contain alphanumeric "
        "characters, underscores and hyphens. No spaces allowed.",
    )

    if "nocombine" not in st.session_state:
        st.session_state.nocombine = None

    st.session_state.nocombine = st.checkbox(
        "Do Not Combine Identical Sequences",
        help="If using STRait Razor data, by default, identical sequences (after removing "
        "flanking sequences) are combined and reads are summed. Checking this will not combine"
        " identical sequences.",
    )

    #####################################################################
    #     STR: Filter Settings to Generate Config File                  #
    #####################################################################

    st.subheader("Filter Settings")

    col1, col2, col3, col4, col5 = st.columns(5)

    if "output_type" not in st.session_state:
        st.session_state.output_type = None

    st.session_state.output_type = {
        "STRmix": "strmix",
        "EuroForMix": "efm",
        "MPSproto": "mpsproto",
    }[
        col1.selectbox(
            "Probabilistic Genotyping Software",
            options=["STRmix", "EuroForMix", "MPSproto"],
            help="Select which probabilistic genotyping software files to create",
        )
    ]

    if "profile_type" not in st.session_state:
        st.session_state.profile_type = None

    st.session_state.profile_type = {"Evidence": "evidence", "Reference": "reference"}[
        col2.selectbox(
            "Profile Type",
            options=["Evidence", "Reference"],
            help="Select the file type (format) to create for the probabilistic genotyping "
            "software.",
        )
    ]

    if "data_type" not in st.session_state:
        st.session_state.data_type = None

    st.session_state.data_type = {"Sequence": "ngs", "CE allele": "ce", "LUS+ allele": "lusplus"}[
        col3.selectbox(
            "Data Type",
            options=["Sequence", "CE allele", "LUS+ allele"],
            help="Select the allele type used to determine sequence type (belowAT, stutter or "
            "typed) and used in the final output file.",
        )
    ]

    if "info" not in st.session_state:
        st.session_state.info = None

    st.session_state.info = st.checkbox(
        "Create Allele Information File",
        value=True,
        help="Create file containing information about each sequence, including sequence type "
        "(belowAT, stutter or typed), stuttering sequence information and metrics involving "
        "stutter and noise.",
    )

    if "separate" not in st.session_state:
        st.session_state.separate = None

    st.session_state.separate = st.checkbox(
        "Create Separate Files for Samples",
        help="If checked, will create individual files for samples; If unchecked, will create "
        "one file with all samples.",
    )

    if "nofilters" not in st.session_state:
        st.session_state.nofilters = None

    st.session_state.nofilters = st.checkbox(
        "Skip All Filtering Steps",
        help="Filtering will not be performed but will still create EFM/MPSproto/STRmix output "
        "files containing all sequences.",
    )

    if "strand" not in st.session_state:
        st.session_state.strand = None

    st.session_state.strand = {"UAS Orientation": "uas", "Forward Strand": "forward"}[
        col4.selectbox(
            "Strand Orientation",
            options=["Forward Strand", "UAS Orientation"],
            help="Indicates the strand orientation in which to report the sequence in the final "
            "output table as some markers are reported in the UAS on the reverse strand. "
            "Selecting the UAS Orientation will report those markers on the reverse strand while"
            " the remaining will be reported on the forward strand. Selecting the Forward Strand "
            "will report all markers on the forward strand orientation. This applies to STRmix "
            "NGS only.",
        )
    ]

    #####################################################################
    #     STR: Generate Config File Based on Settings                   #
    #####################################################################

    # Submit Button Instance
    if st.button("Run lusSTR"):

        # Check if all required fields are filled
        if (
            st.session_state.analysis_software
            and st.session_state.samp_input
            and st.session_state.output
            and st.session_state.wd_dirname
        ):

            # Validate output prefix
            if not validate_prefix(st.session_state.output):
                st.warning(
                    "Please enter a valid output prefix. Only alphanumeric characters, "
                    "underscore, and hyphen are allowed."
                )
                st.stop()  # Stop execution if prefix is invalid

            # Display loading spinner (Continuing Process Checks Above Were Passed)
            with st.spinner("Processing Your Data..."):

                # Construct config data

                config_data = {
                    "analysis_software": st.session_state.analysis_software,
                    "custom_ranges": st.session_state.custom_ranges,
                    "sex": st.session_state.sex,
                    "samp_input": st.session_state.samp_input,
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

                # Generate YAML config file
                generate_config_file(config_data, st.session_state.wd_dirname, "STR")

                # Subprocess lusSTR commands
                command = ["lusstr", "strs", "all"]

                # Specify WD to lusSTR
                if wd_dirname:
                    command.extend(["-w", st.session_state.wd_dirname + "/"])

                # Run lusSTR command in terminal
                try:
                    subprocess.run(command, check=True)
                    st.success(
                        "Config File Generated and lusSTR Executed Successfully! Output Files"
                        "Have Been Saved to Your Designated Directory and Labeled with your "
                        "Specified Prefix"
                    )
                except subprocess.CalledProcessError as e:
                    st.error(f"Error: {e}")
                    st.info(
                        "Please make sure to check the 'How to Use' tab for common error "
                        "resolutions."
                    )

        else:
            st.warning(
                "Please make sure to fill out all required fields (Analysis Software, Input "
                "Directory or File, Prefix for Output, and Specification of Working Directory) "
                "before submitting."
            )
    st.write("---")
    st.write(
        "After running lusSTR, or if lusSTR has been run previously, the user may view and edit "
        "the individual STR marker plots and data."
    )
    st.write(
        "If lusSTR has been previously run, only the above ```Output Folder``` containing the run"
        " files needs to be specified. Other settings will be automatically loaded from the "
        "config.yaml file within the specified folder."
    )
    if "interactive" not in st.session_state:
        st.session_state.interactive = None
    if st.button("See Individual Marker Plots & Data") or st.session_state.interactive:
        st.session_state.interactive = True
        create_settings()
        if st.session_state.custom_ranges:
            file = (
                f"{st.session_state.wd_dirname}/{st.session_state.output}/"
                f"{st.session_state.output}_custom_range"
            )
        else:
            file = f"{wd_dirname}/{st.session_state.output}/{st.session_state.output}"
        try:
            sequence_info = pd.read_csv(f"{file}_sequence_info.csv")
            interactive_setup(sequence_info, file)
        except FileNotFoundError:
            print(f"{file}_sequence_info.csv not found. Please check output folder specification.")


#####################################################################
#                        SNP WORKFLOW                               #
#####################################################################

#####################################################################
# Specify SNP Settings Which Will Be Used to Generate Config File   #
#####################################################################


def show_SNP_page():

    st.title("SNP Workflow")
    st.info(
        "Please Select SNP Settings Below for lusSTR! For Information Regarding the Settings,"
        " See the How to Use Tab."
    )

    # Input File Specification
    st.subheader("Input Files Selection")

    # Ask user if submitting a directory or individual file
    st.info(
        "Please Indicate If You Are Providing An Individual Input File or a Folder Containing "
        "Multiple Input Files"
    )
    input_option = st.radio(
        "Select Input Option:", ("Individual File", "Folder with Multiple Files")
    )

    # Initialize session state if not already initialized
    if "samp_input" not in st.session_state:
        st.session_state.samp_input = None

    # Logic for Path Picker based on user's input option

    if input_option == "Folder with Multiple Files":
        clicked = st.button("Please Select a Folder")
        if clicked:
            dirname = folder_picker_dialog()
            # st.text_input('You Selected The Following folder:', dirname)
            st.session_state.samp_input = dirname

    else:
        clicked_file = st.button("Please Select a File")
        if clicked_file:
            filename = file_picker_dialog()
            # st.text_input('You Selected The Following file:', filename)
            st.session_state.samp_input = filename

    # Display The Selected Path
    if st.session_state.samp_input:
        st.text_input("Location Of Your Input File(s):", st.session_state.samp_input)

    # Store Selected Path to Reference in Config
    samp_input = st.session_state.samp_input

    #####################################################################
    #      SNP: General Software Settings to Generate Config File       #
    #####################################################################

    st.subheader("General Settings")

    col1, col2, col3, col4, col5 = st.columns(5)

    analysis_software = {"UAS": "uas", "STRait Razor v3": "straitrazor"}[
        col1.selectbox(
            "Analysis Software",
            options=["UAS", "STRait Razor v3"],
            help="Indicate the analysis software used prior to lusSTR sex.",
        )
    ]

    output = col2.text_input(
        "Output File Name", "lusstr_output", help="Please specify a name for the created files."
    )

    kit = {"Signature Prep": "sigprep", "Kintelligence": "kintelligence"}[
        col3.selectbox("Library Preparation Kit", options=["Signature Prep", "Kintelligence"])
    ]

    #####################################################################
    #     SNP: Format Settings to Generate Config File                  #
    #####################################################################

    st.subheader("Convert Settings")

    col1, col2, col3, col4, col5 = st.columns(5)

    # -- Select Type (Unique to SNP Workflow)
    types_mapping = {
        "Identify SNPs": "i",
        "Phenotype SNPs": "p",
        "Ancestry SNPs": "a",
        "All SNPs": "all",
    }
    selected_types = col1.multiselect(
        "Select SNP Types:",
        options=types_mapping.keys(),
        help="Select the SNP types to process; can select one or more options",
    )
    types_string = (
        "all"
        if "All" in selected_types
        else ", ".join(types_mapping.get(t, t) for t in selected_types)
    )

    # -- Filter
    nofilters = st.checkbox(
        "Skip all filtering steps",
        help="Specify for no filtering",
    )

    #####################################################################
    #     SNP: Convert Settings to Generate Config File                 #
    #####################################################################

    separate = st.checkbox(
        "Create Separate Files for Samples",
        help="If want to separate samples into individual files for use in EFM",
    )

    strand = {"UAS Orientation": "uas", "Forward Strand": "forward"}[
        col2.selectbox(
            "Strand Orientation",
            options=["UAS Orientation", "Forward Strand"],
            help="Indicate which orientation to report the alleles for the SigPrep SNPs.",
        )
    ]

    # Analytical threshold value
    thresh = col3.number_input("Analytical threshold value:", value=0.03, step=0.01, min_value=0.0)

    #####################################################################
    #     SNP: Specify a Reference File if User Has One                 #
    #####################################################################

    col1, col2, col3 = st.columns(3)

    if "reference" not in st.session_state:
        st.session_state.reference = None

    reference = col1.text_input(
        "Please Specify Your Reference Sample IDs",
        help="List IDs of the samples to be run as references in EFM; default is no "
        "reference samples",
    )

    #####################################################################
    #     SNP: Specify Working Directory                                #
    #####################################################################

    st.subheader("Set Output Folder")

    col1, col2, col3, col4, col5 = st.columns(5)

    # Initialize session state if not already initialized
    if "wd_dirname" not in st.session_state:
        st.session_state.wd_dirname = None

    clicked_wd = col1.button("Please Select An Output Folder")
    if clicked_wd:
        wd = folder_picker_dialog()
        st.session_state.wd_dirname = wd

    # Display selected path
    if st.session_state.wd_dirname:
        st.text_input("Your Specified Output Folder:", st.session_state.wd_dirname)

    #####################################################################
    #     SNP: Generate Config File Based on Settings                   #
    #####################################################################

    # Submit Button Instance
    if st.button("Submit"):

        # Check if all required fields are filled
        if analysis_software and samp_input and output and wd_dirname:

            # Validate output prefix
            if not validate_prefix(output):
                st.warning(
                    "Please enter a valid output prefix. Only alphanumeric characters, "
                    "underscore, and hyphen are allowed."
                )
                st.stop()  # Stop execution if prefix is invalid

            # Display loading spinner (Continuing Process Checks Above Were Passed)
            with st.spinner("Processing Your Data..."):

                # Construct config data

                config_data = {
                    "analysis_software": analysis_software,
                    "samp_input": samp_input,
                    "output": output,
                    "kit": kit,
                    "types": types_string,
                    "thresh": thresh,
                    "separate": separate,
                    "nofilter": nofilters,
                    "strand": strand,
                    "references": None,  # Default value is None
                }

                # If a reference file was specified, add to config
                if reference:
                    config_data["references"] = st.session_state.reference

                # Generate YAML config file
                generate_config_file(config_data, st.session_state.wd_dirname, "SNP")

                # Subprocess lusSTR commands
                command = ["lusstr", "snps", "all"]

                # Specify WD to lusSTR
                if wd_dirname:
                    command.extend(["-w", st.session_state.wd_dirname + "/"])

                # Run lusSTR command in terminal
                try:
                    subprocess.run(command, check=True)
                    st.success(
                        "Config File Generated and lusSTR Executed Successfully! Output Files "
                        "Have Been Saved to Your Designated Directory and Labeled with your "
                        "Specified Prefix"
                    )
                except subprocess.CalledProcessError as e:
                    st.error(f"Error: {e}")
                    st.info(
                        "Please make sure to check the 'How to Use' tab for common error "
                        "resolutions."
                    )

        else:
            st.warning(
                "Please make sure to fill out all required fields (Analysis Software, Input "
                "Directory or File, Prefix for Output, and Specification of Working Directory) "
                "before submitting."
            )


#####################################################################
#                        How To Use Page                            #
#####################################################################


def show_how_to_use_page():

    st.title("Common Errors and Best Practices for Using lusSTR")

    st.header("1. File/Folder Path Formatting")

    st.write(
        "Please ensure that the displayed path accurately reflects your selection. When using the"
        " file or folder picker, navigate to the desired location and click 'OK' to confirm your "
        "selection."
    )

    st.header("2. Specifying Output Prefix")

    st.write(
        "The purpose of specifying the output prefix is for lusSTR to create result files and "
        "folders with that prefix in your working directory. Please ensure that you are following"
        " proper file naming formatting and rules when specifying this prefix. Avoid using "
        "characters such as '/', '', '.', and others. Note: To avoid potential errors, you can "
        "simply use the default placeholder for output."
    )

    st.code("Incorrect: 'working_directory/subfolder/subfolder'\nCorrect: output")

    st.write(
        "Note that some result files may be saved directly in the working directory with the "
        "specified prefix, while others will be populated in a folder labeled with the prefix "
        "in your working directory."
    )
    st.write("Be aware of this behavior when checking for output files.")

    st.header("3. Specifying Output Folder")
    st.write(
        "Please Ensure That You Properly Specify an Output Folder. This is where all lusSTR "
        "output files will be saved. To avoid potential errors, specifying a working directory "
        "is required."
    )

    st.title("About lusSTR")

    st.markdown(
        """

    **_lusSTR Accommodates Four Different Input Formats:_**

    (1) UAS Sample Details Report, UAS Sample Report, and UAS Phenotype Report (for SNP "
    "processing) in .xlsx format (a single file or directory containing multiple files)

    (2) STRait Razor v3 output with one sample per file (a single file or directory containing"
    " multiple files)

    (3) GeneMarker v2.6 output (a single file or directory containing multiple files)

    (4) Sample(s) sequences in CSV format; first four columns must be Locus, NumReads, Sequence, "
    "SampleID; Optional last two columns can be Project and Analysis IDs.


    """,
        unsafe_allow_html=True,
    )


#####################################################################
#                        Contact Page                            #
#####################################################################


def show_contact_page():
    st.title("Contact Us")
    st.write(
        "For any questions or issues, please contact rebecca.mitchell@st.dhs.gov, "
        "daniel.standage@st.dhs.gov, or s.h.syed@email.msmary.edu"
    )


#####################################################################
#     lusSTR RUN                                                    #
#####################################################################

if __name__ == "__main__":
    main()


def subparser(subparsers):
    parser = subparsers.add_parser("gui", help="Launch the Streamlit GUI")
    parser.set_defaults(func=main)
