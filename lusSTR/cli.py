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
from . import format, annot, snps, filter


def format_subparser(subparsers):
    cli = subparsers.add_parser('format')
    cli.add_argument(
        '-o', '--out', metavar='FILE',
        help='File to which output will be written; default is terminal (stdout)'
    )
    cli.add_argument(
        'input',
        help='Input is either a UAS Sample Details Report (in .xlsx format) or a STRait Razor '
        'output file (.txt format). Both single files and directories containing multiple UAS or '
        'STRait Razor files are accepted. If input is the UAS Sample Details Report, use of the '
        '--uas flag is required. If a directory of STRait Razor files is provided, the name of '
        'the directory will be used as the Project and Analysis IDs in the final annotation '
        'table. If a single file is provided, the Project and Analysis IDs will be NA. '
        'STRaitRazor files should be named as the Sample ID, e.g. A001.txt, A002.txt, etc.'
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
        help='file to which output will be written; default is terminal (stdout). If the '
        '--separate flag is used, this will be the name of the directory which the individual '
        'files are written to.'
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
    cli.add_argument(
        '--separate', action='store_true',
        help='This flag will result in the creation of individual output files per sample.'
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
    cli.add_argument(
        '--separate', action='store_true',
        help='This flag will result in the creation of individual output files per sample.'
    )


def filter_subparser(subparsers):
    cli = subparsers.add_parser('filter')
    cli.add_argument(
        'input',
        help='Input is a single lusSTR output file (.txt format)'
    )
    cli.add_argument(
        '--separate', action='store_true',
        help='Used to create separate final output files for each Sample. If not used, a single '
        'file containing all samples will be created.'
    )
    cli.add_argument(
        '--info', action='store_true',
        help='Use to create a text document containing additional information on filtered '
        'sequences and stutter.'
    )
    cli.add_argument(
        '--output-type', dest='output', choices=['efm', 'strmix'], default='efm',
        help='Choose the file format of the output file, either "efm" or "strmix". '
        'Default is efm.'
    )
    cli.add_argument(
        '--no-filters', dest='nofilters', action='store_true',
        help='Used to skip all filtering steps. All input alleles will be included in the output.'
    )
    cli.add_argument(
        '--out', '-o', metavar='FILE',
        help='Name of output file containing all samples for EFM or name/path of directory for '
        'STRmix. If separate files are specified for EFM, the sample ID will be used as the '
        'filename. Output files are in CSV format.'
    )
    cli.add_argument(
        '--profile-type', dest='profile', choices=['evidence', 'reference'], default='evidence',
        help='Choose the type of profile, either evidence or reference. Default is evidence.'
    )
    cli.add_argument(
        '--data-type', dest='data', choices=['ngs', 'ce'], default='ce',
        help='Choose the type of data, either ngs or ce. Default is ce.'
    )


mains = {
    'format': lusSTR.format.main,
    'annotate': lusSTR.annot.main,
    'snps': lusSTR.snps.main,
    'filter': lusSTR.filter.main
}

subparser_funcs = {
    'format': format_subparser,
    'annotate': annot_subparser,
    'snps': snps_subparser,
    'filter': filter_subparser
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
