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

from .home import HomePage
from .howto import HowToPage
from .contact import ContactPage
from .strs import STRWorkflow
import streamlit as st
from streamlit_option_menu import option_menu


apps = {
    "Home": HomePage,
    "STRs": STRWorkflow,
    "SNPs": None,
    "How to Use": HowToPage,
    "Contact": ContactPage,
}


def initialize():
    st.set_page_config(layout="wide", initial_sidebar_state="collapsed")
    selected = option_menu(
        menu_title=None,
        options=["Home", "STRs", "SNPs", "How to Use", "Contact"],
        icons=["house", "gear", "gear-fill", "book", "envelope"],
        menu_icon="cast",
        default_index=0,
        orientation="horizontal",
    )
    appname = str(selected)
    return apps[appname]
