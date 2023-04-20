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
import os
from pathlib import Path
from pkg_resources import resource_filename
import yaml


def main(args):
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    final_dest = f"{args.workdir}/config.yaml"
    config = resource_filename("lusSTR", "data/config.yaml")
    final_config = edit_config(config, args)
    with open(final_dest, "w") as file:
        yaml.dump(final_config, file)

def edit_config(config, args):
    with open(config, "r") as file:
        data = yaml.safe_load(file)
    if args.straitrazor:
        data["uas"] = False
    if args.powerseq:
        data["kit"] = "powerseq"
    if args.input:
        data["samp_input"] = args.input
    else:
        data["samp_input"] = os.getcwd()
    if args.out:
        data["output"] = args.out
    if args.sex:
        data["sex"] = True
    if args.separate:
        data["separate"] = True
    if args.nocombine:
        data["nocombine"] = True
    if args.nofiltering:
        data["nofilters"] = True
    if args.noinfo:
        data["info"] = False
    if args.reference:
        data["profile_type"] = "reference"
    if args.ce:
        data["data_type"] = "ce"
    if args.efm:
        data["output_type"] = "efm"
    return data


def subparser(subparsers):
    p = subparsers.add_parser("config", description="Create config file for running STR pipeline")
    p.add_argument(
        "-w", "--workdir", metavar="W", default=".",
        help="directory to add config file; default is current working directory")
    p.add_argument(
        "--straitrazor", action="store_true",
        help="Use if sequences have been previously run through STRait Razor."
    )
    p.add_argument("--input", help="Input file or directory")
    p.add_argument("--out", "-o", help="Output file/directory name")
    p.add_argument(
        "--powerseq", action="store_true",
        help="Use to indicate sequences were created using the PowerSeq Kit."
    )
    p.add_argument(
        "--sex", action="store_true",
        help="Use if including the X and Y STR markers. Separate reports for these markers "
        "will be created.",
    )
    p.add_argument(
        "--nocombine", action="store_true",
        help="Do not combine read counts for duplicate sequences within the UAS region "
        "during the 'annotate' step. By default, read counts are combined for sequences "
        "not run through the UAS.",
    )
    p.add_argument(
        "--reference", action="store_true", 
        help="Use for creating Reference profiles"
    )
    p.add_argument("--efm", action="store_true",help="Use to create EuroForMix profiles")
    p.add_argument("--ce", action="store_true", help="Use for CE data")
    p.add_argument(
        "--noinfo", action="store_true", 
        help="Use to not create the Sequence Information File in the 'filter' step"
    )
    p.add_argument(
        "--separate", action="store_true", 
        help="Use to separate EFM profiles in the 'filter' step."
    )
    p.add_argument(
        "--nofiltering", action="store_true", 
        help="Use to perform no filtering during the 'filter' step"
    )
