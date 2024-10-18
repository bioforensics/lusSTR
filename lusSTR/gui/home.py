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


def home_page_display():
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
