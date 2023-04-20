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
import snakemake
from pkg_resources import resource_filename

## placeholder until I update this

def main(args):
    snakefile = resource_filename("lusSTR", "workflows/snps.smk")
    pretarget = "annotate" if args.filter else "filter"
    result = snakemake.snakemake(
        snakefile, config=args.config, targets=pretarget,
        workdir=args.work_dir
    )
    if result is not True:
        raise NotImplementedError('SNP workflow implementation pending')

def subparser(subparsers):
    p = subparsers.add_parser("snps", description="Running the entire STR pipeline (format, annotate and filter)")
    p.add_argument("--config", default="config.yaml", help="config file used to identify settings.")
    p.add_argument("-w", "--workdir", metavar="W", default=".", help="working directory")
    p.add_argument("--skip-filter", dest="filter", action = "store_true", help="Skip filtering step")
    