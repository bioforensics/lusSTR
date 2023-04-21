import argparse
import lusSTR
from lusSTR.cli import config
from lusSTR.cli import strs
from lusSTR.cli import snps
import snakemake


mains = {
    "config": config.main,
    "strs": strs.main,
    "snps": snps.main
}

subparser_funcs = {
    "config": config.subparser,
    "strs": strs.subparser,
    "snps": snps.subparser
}


def main(args=None):
    if args is None: 
        args = get_parser().parse_args()
    if args.subcmd is None:
        get_parser().parse_args(["-h"])
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
