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

from argparse import ArgumentParser
from importlib.resources import files
import lusSTR
from lusSTR.cli import config, gui, strs, snps
import streamlit.web.cli as stcli
import sys

mains = {"config": config.main, "strs": strs.main, "snps": snps.main}

subparser_funcs = {
    "config": config.subparser,
    "strs": strs.subparser,
    "snps": snps.subparser,
    "gui": gui.subparser,
}


def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    if args.subcmd is None:
        get_parser().parse_args(["-h"])
    elif args.subcmd == "gui":
        gui_path = files("lusSTR") / "cli" / "gui.py"
        sys.argv = ["streamlit", "run", str(gui_path)]
        sys.exit(stcli.main())
    else:
        mainmethod = mains[args.subcmd]
        result = mainmethod(args)
        return result


def get_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "-v", "--version", action="version", version="lusSTR v" + lusSTR.__version__
    )
    subcommandstr = ", ".join(sorted(subparser_funcs.keys()))
    subparsers = parser.add_subparsers(dest="subcmd", metavar="subcmd", help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return parser


if __name__ == "__main__":
    main()
