from . import all
from . import snps
import snakemake


mains = {
    "all": all.main,
    "snps": snps.main
}

subparser_funcs = {
    "all": all.subparser,
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
    subparsers = parser.add_subparsers(dest="subcmd", metavar="subcmd", help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return parser
