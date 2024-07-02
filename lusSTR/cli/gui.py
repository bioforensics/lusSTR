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

import json
import importlib.resources
import streamlit as st
from streamlit_option_menu import option_menu
import yaml
import subprocess
import os
import re

# ------ Packages For File/Folder Directory Selection --------- #

import tkinter as tk
from tkinter import filedialog

# Create a global Tkinter root window
root = tk.Tk()
root.withdraw()  # Hide the root window

#################################################################
#                        Functions                              #
#################################################################

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
        lusSTR is an end-to-end workflow for processing human forensic data (STRs and SNPs) derived from Next Generation Sequencing (NGS) data for use in probabilistic genotyping software.
        For more information on lusSTR, visit our [GitHub page](https://github.com/bioforensics/lusSTR).
        """,
        unsafe_allow_html=True,
    )

    st.info("Please Select One of the Tabs Above to Get Started on Processing Your Data!")


#####################################################################
#                        STR WORKFLOW                               #
#####################################################################

#####################################################################
# Specify STR Settings Which Will Be Used to Generate Config File   #
#####################################################################


def show_STR_page():

    st.title("STR Workflow")
    st.info(
        "Please Select STR Settings Below for lusSTR! For Information Regarding the Settings, See the How to Use Tab."
    )

    # Input File Specification
    st.subheader("Input Files Selection")

    # Ask user if submitting a directory or individual file
    st.info(
        "Please Indicate If You Are Providing An Individual Input File or a Folder Containing Multiple Input Files"
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
            st.session_state.samp_input = dirname

    else:
        clicked_file = st.button("Please Select a File")
        if clicked_file:
            filename = file_picker_dialog()
            st.session_state.samp_input = filename

    # Display The Selected Path
    if st.session_state.samp_input:
        st.text_input("Location Of Your Input File(s):", st.session_state.samp_input)

    # Store the Selected Path to Reference in Config
    samp_input = st.session_state.samp_input

    #####################################################################
    #      STR: General Software Settings to Generate Config File       #
    #####################################################################

    st.subheader("General Settings")

    col1, col2, col3, col4, col5 = st.columns(5)

    analysis_software = {
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

    sex = st.checkbox(
        "Include X- and Y-STRs",
        help="Check the box to include X- and Y-STRs, otherwise leave unchecked.",
    )

    kit = {"ForenSeq Signature Prep": "forenseq", "PowerSeq 46GY": "powerseq"}[
        col2.selectbox(
            "Library Preparation Kit", options=["ForenSeq Signature Prep", "PowerSeq 46GY"], help="Specify the library preparation kit used to generate the sequences."
        )
    ]

    output = col3.text_input(
        "Output File Name", "lusstr_output", help="Please specify a name for the created files. It can only contain alphanumeric characters, underscores and hyphens. No spaces allowed."
    )


    nocombine = st.checkbox(
        "Do Not Combine Identical Sequences",
        help="If using STRait Razor data, by default, identical sequences (after removing flanking sequences) are combined and reads are summed. Checking this will not combine identical sequences.",
    )

    #####################################################################
    #     STR: Filter Settings to Generate Config File                  #
    #####################################################################

    st.subheader("Filter Settings")

    col1, col2, col3, col4, col5 = st.columns(5)

    output_type = {"STRmix": "strmix", "EuroForMix": "efm", "MPSproto": "mpsproto"}[
        col1.selectbox(
            "Probabilistic Genotyping Software",
            options=["STRmix", "EuroForMix", "MPSproto"],
            help="Select which probabilistic genotyping software files to create",
        )
    ]

    profile_type = {"Evidence": "evidence", "Reference": "reference"}[
        col2.selectbox(
            "Profile Type",
            options=["Evidence", "Reference"],
            help="Select the file type (format) to create for the probabilistic genotyping software.",
        )
    ]

    data_type = {"Sequence": "ngs", "CE allele": "ce", "LUS+ allele": "lusplus"}[
        col3.selectbox(
            "Data Type",
            options=["Sequence", "CE allele", "LUS+ allele"],
            help="Select the allele type used to determine sequence type (belowAT, stutter or typed) and used in the final output file.",
        )
    ]

    info = st.checkbox(
        "Create Allele Information File",
        value=True,
        help="Create file containing information about each sequence, including sequence type (belowAT, stutter or typed), stuttering sequence information and metrics involving stutter and noise.",
    )

    separate = st.checkbox(
        "Create Separate Files for Samples",
        help="If checked, will create individual files for samples; If unchecked, will create one file with all samples.",
    )

    nofilters = st.checkbox(
        "Skip all filtering steps",
        help="Will not perform filtering; will still create EFM/MPSproto/STRmix output files",
    )

    strand = {"UAS Orientation": "uas", "Forward Strand": "forward"}[
        col4.selectbox(
            "Strand Orientation",
            options=["UAS Orientation", "Forward Strand"],
            help="Indicates the strand orientation in which to report the sequence in the final output table; for STRmix NGS only.",
        )
    ]

    #####################################################################
    #     STR: Specify Working Directory                                #
    #####################################################################

    st.subheader("Output Folder Selection")

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

    # Store Selected Path to Reference in Config
    wd_dirname = st.session_state.wd_dirname

    #####################################################################
    #     STR: Generate Config File Based on Settings                   #
    #####################################################################

    # Submit Button Instance
    if st.button("Submit"):

        # Check if all required fields are filled
        if analysis_software and samp_input and output and wd_dirname:

            # Validate output prefix
            if not validate_prefix(output):
                st.warning(
                    "Please enter a valid output prefix. Only alphanumeric characters, underscore, and hyphen are allowed."
                )
                st.stop()  # Stop execution if prefix is invalid

            # Display loading spinner (Continuing Process Checks Above Were Passed)
            with st.spinner("Processing Your Data..."):

                # Construct config data

                config_data = {
                    "analysis_software": analysis_software,
                    "sex": sex,
                    "samp_input": samp_input,
                    "output": output,
                    "kit": kit,
                    "nocombine": nocombine,
                    "output_type": output_type,
                    "profile_type": profile_type,
                    "data_type": data_type,
                    "info": info,
                    "separate": separate,
                    "nofilters": nofilters,
                    "strand": strand,
                }

                # Generate YAML config file
                generate_config_file(config_data, wd_dirname, "STR")

                # Subprocess lusSTR commands
                command = ["lusstr", "strs", "all"]

                # Specify WD to lusSTR
                if wd_dirname:
                    command.extend(["-w", wd_dirname + "/"])

                # Run lusSTR command in terminal
                try:
                    subprocess.run(command, check=True)
                    st.success(
                        "Config File Generated and lusSTR Executed Successfully! Output Files Have Been Saved to Your Designated Directory and Labeled with your Specified Prefix"
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


#####################################################################
#                        SNP WORKFLOW                               #
#####################################################################

#####################################################################
# Specify SNP Settings Which Will Be Used to Generate Config File   #
#####################################################################


def show_SNP_page():

    st.title("SNP Workflow")
    st.info(
        "Please Select SNP Settings Below for lusSTR! For Information Regarding the Settings, See the How to Use Tab."
    )

    # Input File Specification
    st.subheader("Input Files Selection")

    # Ask user if submitting a directory or individual file
    st.info(
        "Please Indicate If You Are Providing An Individual Input File or a Folder Containing Multiple Input Files"
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
        help="List IDs of the samples to be run as references in EFM; default is no reference samples",
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

    # Store Selected Path to Reference in Config
    wd_dirname = st.session_state.wd_dirname

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
                    "Please enter a valid output prefix. Only alphanumeric characters, underscore, and hyphen are allowed."
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
                    config_data["references"] = reference

                # Generate YAML config file
                generate_config_file(config_data, wd_dirname, "SNP")

                # Subprocess lusSTR commands
                command = ["lusstr", "snps", "all"]

                # Specify WD to lusSTR
                if wd_dirname:
                    command.extend(["-w", wd_dirname + "/"])

                # Run lusSTR command in terminal
                try:
                    subprocess.run(command, check=True)
                    st.success(
                        "Config File Generated and lusSTR Executed Successfully! Output Files Have Been Saved to Your Designated Directory and Labeled with your Specified Prefix"
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


#####################################################################
#                        How To Use Page                            #
#####################################################################


def show_how_to_use_page():

    st.title("Common Errors and Best Practices for Using lusSTR")

    st.header("1. File/Folder Path Formatting")

    st.write(
        "Please ensure that the displayed path accurately reflects your selection. When using the file or folder picker, navigate to the desired location and click 'OK' to confirm your selection."
    )

    st.header("2. Specifying Output Prefix")

    st.write(
        "The purpose of specifying the output prefix is for lusSTR to create result files and folders with that prefix in your working directory. Please ensure that you are following proper file naming formatting and rules when specifying this prefix. Avoid using characters such as '/', '', '.', and others. Note: To avoid potential errors, you can simply use the default placeholder for output."
    )

    st.code("Incorrect: 'working_directory/subfolder/subfolder'\nCorrect: output")

    st.write(
        "Note that some result files may be saved directly in the working directory with the specified prefix, while others will be populated in a folder labeled with the prefix in your working directory."
    )
    st.write("Be aware of this behavior when checking for output files.")

    st.header("3. Specifying Output Folder")
    st.write(
        "Please Ensure That You Properly Specify an Output Folder. This is where all lusSTR output files will be saved. To avoid potential errors, specifying a working directory is required."
    )

    st.title("About lusSTR")

    st.markdown(
        """

    **_lusSTR Accommodates Four Different Input Formats:_**

    (1) UAS Sample Details Report, UAS Sample Report, and UAS Phenotype Report (for SNP processing) in .xlsx format (a single file or directory containing multiple files)

    (2) STRait Razor v3 output with one sample per file (a single file or directory containing multiple files)

    (3) GeneMarker v2.6 output (a single file or directory containing multiple files)

    (4) Sample(s) sequences in CSV format; first four columns must be Locus, NumReads, Sequence, SampleID; Optional last two columns can be Project and Analysis IDs.


    """,
        unsafe_allow_html=True,
    )


#####################################################################
#                        Contact Page                            #
#####################################################################


def show_contact_page():
    st.title("Contact Us")
    st.write(
        "For any questions or issues, please contact rebecca.mitchell@st.dhs.gov, daniel.standage@st.dhs.gov, or s.h.syed@email.msmary.edu"
    )


#####################################################################
#     lusSTR RUN                                                    #
#####################################################################

if __name__ == "__main__":
    main()


def subparser(subparsers):
    parser = subparsers.add_parser("gui", help="Launch the Streamlit GUI")
    parser.set_defaults(func=main)
