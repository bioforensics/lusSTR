import argparse
import lusSTR
from lusSTR.cli import strs
from lusSTR.cli import snps
import snakemake


mains = {
    "strs": strs.main,
    "snps": snps.main
}

subparser_funcs = {
    "strs": strs.subparser,
    "snps": snps.subparser
}


def main(args=None):
    if args is None: 
        args = get_parser().parse_args()
    if args.cmd is None:
        get_parser().parse_args(["-h"])
    mainmethod = mains[args.cmd]
    result = mainmethod(args)
    return result


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--version", action="version", version="lusSTR v" + lusSTR.__version__
    )
    subcommandstr = ", ".join(sorted(subparser_funcs.keys()))
    subparsers = parser.add_subparsers(dest="cmd", metavar="cmd", help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return parser
