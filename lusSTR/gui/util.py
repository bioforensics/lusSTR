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

from pathlib import Path
import re
import yaml


def generate_config_file(config_data, working_directory, workflow_type):
    if workflow_type == "STR":
        config_filename = "config.yaml"
    elif workflow_type == "SNP":
        config_filename = "snp_config.yaml"
    else:
        raise ValueError("Invalid workflow type. Please specify either 'STR' or 'SNP'.")
    config_path = Path(working_directory) / config_filename
    with open(config_path, "w") as file:
        yaml.dump(config_data, file)


def validate_prefix(prefix):
    return re.match(r"^[A-Za-z0-9_-]+$", prefix) is not None
