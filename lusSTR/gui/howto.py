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

import streamlit as st


def howto_page_display():
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
