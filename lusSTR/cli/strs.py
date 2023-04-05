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
from snakemake import snakemake


def main(args):
    pretarget = args.target if args.target != "all" else "filter"
    workdir = args.workdir
    result = snakemake(
        lusSTR.snakefile(workflow="strs"), targets=[pretarget], workdir=workdir
    )
    if result is not True:
        raise SystemError('Snakemake failed')

def subparser(subparsers):
    p = subparsers.add_parser("strs", description="Running the entire STR pipeline (format, annotate and filter)")
    p.add_argument('target', help="Steps to run; 'all' will run all steps of the STR workflow.")
    p.add_argument("-w", "--workdir", metavar="W", default=".", help="working directory")
    #p.add_argument("--skip-filter", dest="nofilter", action = "store_true", help="Skip filtering step")
    
