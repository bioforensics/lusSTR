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

from .contact import contact_page_display
from .home import home_page_display
from .howto import howto_page_display
from .snps import snp_workflow_display
from .strs import str_workflow_display
import streamlit as st
from streamlit_option_menu import option_menu


pages = {
    "Home": home_page_display,
    "STRs": str_workflow_display,
    "SNPs": snp_workflow_display,
    "How to Use": howto_page_display,
    "Contact": contact_page_display,
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
    renderer = pages[appname]
    renderer()
