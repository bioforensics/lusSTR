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

from importlib.resources import files
import streamlit as st


def home_page():
    image_path = files("lusSTR") / "cli" / "logo.png"
    left_column, center_column, right_column = st.columns([2.5, 8, 2.5])
    with center_column:
        st.image(str(image_path), use_column_width="auto")

    # CSS to hide full-screen button
    hide_img_fs = """
    <style>
    button[title="View fullscreen"]{
    visibility: hidden;}
    </style>
    """
    st.markdown(hide_img_fs, unsafe_allow_html=True)

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


def how_to_use_page():
    st.title("Common Errors and Best Practices for Using lusSTR")
    st.header("1. File/Folder Path Formatting")
    st.write(
        "Please ensure that the displayed path accurately reflects your selection. When using "
        "the file or folder picker, navigate to the desired location and click 'OK' to "
        "confirm your selection."
    )
    st.header("2. Specifying Output Prefix")
    st.write(
        "The purpose of specifying the output prefix is for lusSTR to create result files and "
        "folders with that prefix in your working directory. Please ensure that you are "
        "following proper file naming formatting and rules when specifying this prefix. Avoid "
        "using characters such as '/', '', '.', and others. Note: To avoid potential errors, "
        "you can simply use the default placeholder for output."
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
        "output files will be saved. To avoid potential errors, specifying a working "
        "directory is required."
    )
    st.title("About lusSTR")
    st.markdown("""
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


def contact_page():
    st.title("Contact Us")
    st.write("For any questions or issues, please contact rebecca.mitchell@st.dhs.gov or daniel.standage@st.dhs.gov.")
