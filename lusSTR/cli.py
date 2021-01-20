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
from . import format, annot, snps


def format_subparser(subparsers):
    cli = subparsers.add_parser('format')
    cli.add_argument(
        '-o', '--out', metavar='FILE',
        help='File to which output will be written; default is terminal (stdout)'
    )
    cli.add_argument(
        'input',
        help='Input is either a single file (UAS Sample Details Report, in .xlsx format) or a '
        'directory of STRait Razor output files. If input is the UAS Sample Details Report '
        '(in .xlsx format), use of the --uas flag is required. If STRait Razor output is '
        'used, the name of the provided directory will be used as the Analysis ID in the '
        'final annotation table. Output files within the directory should be named as such: '
        'SampleID_STRaitRazor.txt (e.g. A001_STRaitRazor.txt).'
    )
    cli.add_argument(
        '--uas', action='store_true',
        help='Use if sequences have been previously run through the ForenSeq UAS.'
    )
    cli.add_argument(
        '--include-sex', dest='sex', action='store_true',
        help='Use if including the X and Y STR markers'
    )


def annot_subparser(subparsers):
    cli = subparsers.add_parser('annotate')
    cli.add_argument(
        '-o', '--out', metavar='FILE',
        help='file to which output will be written; default is terminal (stdout)'
    )
    cli.add_argument(
        'input', help='sample(s) in CSV format; first four columns must be Locus, NumReads, '
        'Sequence, SampleID; Optional last two columns can be Project and Analysis.'
    )
    cli.add_argument(
        '--kit', choices=['forenseq', 'powerseq'], default='forenseq',
        help='Kit used to develop sequences; only forenseq or powerseq accepted;'
        'default = forenseq'
    )
    cli.add_argument(
        '--uas', action='store_true',
        help='Use if sequences have been run through the ForenSeq UAS.'
    )
    cli.add_argument(
        '--nocombine', dest='combine', action='store_false',
        help='Do not combine read counts for duplicate sequences within the UAS region. '
        'By default, read counts are combined for sequences not run through the UAS.'
    )
    cli.add_argument(
        '--include-sex', dest='sex', action='store_true',
        help='Use if including the X and Y STR markers. Separate reports for these markers '
        'will be created.'
    )


def snps_subparser(subparsers):
    cli = subparsers.add_parser('snps')
    cli.add_argument(
        '-o', '--out', metavar='FILE',
        help='file to which output will be written; default is terminal (stdout)'
    )
    cli.add_argument(
        'input',
        help='Input is either a directory of either UAS output files (Sample Details Report and '
        'Phenotype Report) or of STRait Razor output files. If input is the UAS output file(s) '
        '(in .xlsx format), use of the --uas flag is required. If STRait Razor output is '
        'used, the name of the provided directory will be used as the Analysis ID in the '
        'final annotation table.'
    )
    cli.add_argument(
        '--type', choices=['all', 'p', 'i'], default='i',
        help='Specify the type of SNPs to include in the final report. "p" will include only the '
        'Phenotype and Ancestry SNPs; "i" will include only the Identity SNPs; and "all" will '
        'include all SNPs. Default is Identity SNPs only (i).'
    )
    cli.add_argument(
        '--uas', action='store_true',
        help='Use if sequences have been run through the ForenSeq UAS.'
    )


mains = {
    'format': lusSTR.format.main,
    'annotate': lusSTR.annot.main,
    'snps': lusSTR.snps.main,
}

subparser_funcs = {
    'format': format_subparser,
    'annotate': annot_subparser,
    'snps': snps_subparser,
}


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--version', action='version', version='lusSTR v' + lusSTR.__version__
    )
    subcommandstr = ', '.join(sorted(subparser_funcs.keys()))
    subparsers = parser.add_subparsers(dest='subcmd', metavar='subcmd', help=subcommandstr)
    for func in subparser_funcs.values():
        func(subparsers)
    return parser
