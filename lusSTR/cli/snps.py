# -------------------------------------------------------------------------------------------------
# Copyright (c) 2023, DHS.
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
    pretarget = args.target if args.target != "all" else "convert"
    workdir = args.workdir
    result = snakemake(
        lusSTR.snakefile(workflow="snps"), targets=[pretarget], workdir=workdir, verbose=True
    )
    if result is not True:
        raise SystemError("Snakemake failed")


def subparser(subparsers):
    p = subparsers.add_parser("snps", description="Running the SNP pipeline")
    p.add_argument(
        "target",
        choices=["format", "all"],
        help="Steps to run. Specifying 'format' will run only 'format'. Specifying "
        "'all' will run all steps of the SNP workflow ('format' and 'convert').",
    )
    p.add_argument("-w", "--workdir", metavar="W", default=".", help="working directory")
