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
import json
import streamlit as st
from subprocess import run


def folder():
    script = files("lusSTR") / "scripts" / "folder_selector.py"
    result = run(["python", script], capture_output=True, text=True)
    if result.returncode == 0:
        folder_data = json.loads(result.stdout)
        folder_path = folder_data.get("folder_path")
        if folder_path:
            return folder_path
        else:
            st.error("No folder selected")
    else:
        st.error("Error selecting folder")


def file():
    script = files("lusSTR") / "scripts" / "file_selector.py"
    result = run(["python", script], capture_output=True, text=True)
    if result.returncode == 0:
        file_data = json.loads(result.stdout)
        file_path = file_data.get("file_path")
        if file_path:
            return file_path
        else:
            st.error("No folder selected")
    else:
        st.error("Error selecting folder")
