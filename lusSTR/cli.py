#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import lusSTR
from . import format, annot


def format_subparser(subparsers):
    cli = subparsers.add_parser('format')


def annot_subparser(subparsers):
    cli = subparsers.add_parser('annotate')


mains = {
    'format': lusSTR.format.main,
    'annotate': lusSTR.annot.main,
}

subparser_funcs = {
    'format': format_subparser,
    'annotate': annot_subparser,
}


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o', '--out', metavar='FILE',
        help='file to which output will be written; default is terminal (stdout)'
    )
    parser.add_argument(
        'input', help='sample(s) in CSV format; first four columns must be Locus, NumReads, '
        'Sequence, SampleID'
    )
    subcommandstr = ', '.join(sorted(subparser_funcs.keys()))
    subparsers = parser.add_subparsers(dest='subcmd', metavar='subcmd', help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return parser
