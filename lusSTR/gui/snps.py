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
import streamlit as st
import subprocess


def snp_workflow_display():
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
    if "samp_input_snp" not in st.session_state:
        st.session_state.samp_input_snp = None

    # Logic for Path Picker based on user's input option

    if input_option == "Folder with Multiple Files":
        clicked = st.button("Please Select a Folder")
        if clicked:
            st.session_state.samp_input_snp = select.folder()

    else:
        clicked_file = st.button("Please Select a File")
        if clicked_file:
            st.session_state.samp_input_snp = select.file()

    # Display The Selected Path
    if st.session_state.samp_input_snp:
        st.text_input("Location Of Your Input File(s):", st.session_state.samp_input_snp)

    # Store Selected Path to Reference in Config
    samp_input_snp = st.session_state.samp_input_snp

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
    if "wd_dirname_snp" not in st.session_state:
        st.session_state.wd_dirname_snp = None

    clicked_wd = col1.button("Please Select An Output Folder")
    if clicked_wd:
        st.session_state.wd_dirname_snp = select.folder()

    # Display selected path
    if st.session_state.wd_dirname_snp:
        st.text_input("Your Specified Output Folder:", st.session_state.wd_dirname_snp)

    #####################################################################
    #     SNP: Generate Config File Based on Settings                   #
    #####################################################################

    # Submit Button Instance
    if st.button("Submit"):

        # Check if all required fields are filled
        if (
            analysis_software
            and st.session_state.samp_input_snp
            and output
            and st.session_state.wd_dirname_snp
        ):

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
                    "samp_input_snp": samp_input_snp,
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
                generate_config_file(config_data, st.session_state.wd_dirname_snp, "SNP")

                # Subprocess lusSTR commands
                command = ["lusstr", "snps", "all"]

                # Specify WD to lusSTR
                if st.session_state.wd_dirname_snp:
                    command.extend(["-w", st.session_state.wd_dirname_snp + "/"])

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
