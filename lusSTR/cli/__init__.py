import argparse
import os
import subprocess
import lusSTR
from lusSTR.cli import config
from lusSTR.cli import strs
from lusSTR.cli import snps
from lusSTR.cli import gui
import snakemake

mains = {
    "config": config.main,
    "strs": strs.main,
    "snps": snps.main,
    "gui": gui.main
}

subparser_funcs = {
    "config": config.subparser,
    "strs": strs.subparser,
    "snps": snps.subparser,
    "gui": gui.subparser
}

def main(args=None):
    if args is None:
        args = get_parser().parse_args()
    if args.subcmd is None:
        get_parser().parse_args(["-h"])
    elif args.subcmd == "gui":  
        # Get the directory containing the script (cli folder)
        script_dir = os.path.dirname(os.path.realpath(__file__))
        # Construct the path to gui.py relative to the script directory
        gui_path = os.path.join(script_dir, "gui.py")
        # Call streamlit run command
        subprocess.run(["streamlit", "run", gui_path])
    else:
        mainmethod = mains[args.subcmd]
        result = mainmethod(args)
        return result

def get_parser():
    parser = argparse.ArgumentParser()
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
