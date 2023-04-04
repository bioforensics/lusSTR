# -------------------------------------------------------------------------------------------------
# Copyright (c) 2020, DHS.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR) and is licensed under
# the BSD license: see LICENSE.txt.
#
# This software was prepared for the Department of Homeland Security (DHS) by the Battelle National
# Biodefense Institute, LLC (BNBI) as part of contract HSHQDC-15-C-00064 to manage and operate the
# National Biodefense Analysis and Countermeasures Center (NBACC), a Federally Funded Research and
# Development Center.
# -------------------------------------------------------------------------------------------------

import argparse
import lusSTR
from pathlib import Path
from pkg_resources import resource_filename
import shutil


def main(args):
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    final_dest = f"{args.workdir}/config.yaml"
    config = resource_filename("lusSTR", "data/config.yaml")
    shutil.copyfile(config, final_dest)


def subparser(subparsers):
    p = subparsers.add_parser("config", description="Create config file for running STR pipeline")
    p.add_argument(
        "-w", "--workdir", metavar="W", default=".",
        help="directory to add config file; default is current working directory")
    
