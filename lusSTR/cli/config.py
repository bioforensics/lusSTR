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
import os
from pathlib import Path
from pkg_resources import resource_filename
import yaml


def main(args):
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    if args.snps:
        final_dest = f"{args.workdir}/snp_config.yaml"
        config = resource_filename("lusSTR", "data/snp_config.yaml")
        final_config = edit_snp_config(config, args)
    else:
        final_dest = f"{args.workdir}/config.yaml"
        config = resource_filename("lusSTR", "data/config.yaml")
        final_config = edit_str_config(config, args)
    with open(final_dest, "w") as file:
        yaml.dump(final_config, file)


def edit_snp_config(config, args):
    with open(config, "r") as file:
        data = yaml.safe_load(file)
        if args.straitrazor:
            data["uas"] = False
        if args.input:
            data["samp_input"] = args.input
        if args.out:
            data["output"] = args.out
        if args.snptype:
            data["types"] = args.snptype
        if args.kintelligence:
            data["kit"] = "kintelligence"
        if args.separate:
            data["separate"] = True
        if args.nofiltering:
            data["nofilter"] = True
        if args.ref:
            data["references"] = args.ref
        else:
            data["references"] = None
        if args.strand:
            data["strand"] = args.strand
        if args.input:
            data["samp_input"] = args.input
        else:
            data["samp_input"] = os.getcwd()
        return data


def edit_str_config(config, args):
    with open(config, "r") as file:
        data = yaml.safe_load(file)
    if args.straitrazor:
        data["uas"] = False
    if args.powerseq:
        data["kit"] = "powerseq"
    if args.input:
        data["samp_input"] = args.input
    else:
        data["samp_input"] = os.getcwd()
    if args.out:
        data["output"] = args.out
    if args.sex:
        data["sex"] = True
    if args.separate:
        data["separate"] = True
    if args.nocombine:
        data["nocombine"] = True
    if args.nofiltering:
        data["nofilters"] = True
    if args.noinfo:
        data["info"] = False
    if args.reference:
        data["profile_type"] = "reference"
    if args.efm:
        data["output_type"] = "efm"
    if args.strand:
        data["strand"] = args.strand
    data["data_type"] = args.allele
    return data


def subparser(subparsers):
    p = subparsers.add_parser("config", description="Create config file for running lusSTR")
    all_args = p.add_argument_group("General Settings")
    all_args.add_argument(
        "-w", "--workdir", metavar="W", default=".",
        help="directory to add config file; default is current working directory")
    all_args.add_argument(
        "--straitrazor", action="store_true",
        help="Use if sequences have been previously run through STRait Razor."
    )
    all_args.add_argument("--input", help="Input file or directory")
    all_args.add_argument("--out", "-o", help="Output file/directory name")
    all_args.add_argument(
        "--nofiltering", action="store_true", 
        help="For STRs, use to perform no filtering during the 'filter' step. For SNPs, "
        "only alleles specified as 'Typed' by the UAS will be included at the 'format' step."
    )
    all_args.add_argument(
        "--strand", choices=["uas", "forward"],
        help="Specify the strand orientation for the final output files. UAS orientation is "
        "default for STRs; forward strand is default for SNPs."
    )

    str_args = p.add_argument_group("Setting Specific for STR Workflow")
    str_args.add_argument(
        "--powerseq", action="store_true",
        help="Use to indicate sequences were created using the PowerSeq Kit."
    )
    str_args.add_argument(
        "--sex", action="store_true",
        help="Use if including the X and Y STR markers. Separate reports for these markers "
        "will be created.",
    )
    str_args.add_argument(
        "--nocombine", action="store_true",
        help="Do not combine read counts for duplicate sequences within the UAS region "
        "during the 'convert' step. By default, read counts are combined for sequences "
        "not run through the UAS.",
    )
    str_args.add_argument(
        "--reference", action="store_true", 
        help="Use for creating Reference profiles for STR workflow"
    )
    str_args.add_argument("--efm", action="store_true",help="Use to create EuroForMix profiles")
    str_args.add_argument(
        "--allele", choices=["ce", "ngs", "lusplus"], default="ce",
        help="Specify the allele type for the files generated for use in STRmix or EuroForMix. "
        "Options for STRmix: 'ce' or 'ngs'. Options for EuroForMix: 'ce' or 'lusplus'."
    )
    str_args.add_argument(
        "--noinfo", action="store_true", 
        help="Use to not create the Sequence Information File in the 'filter' step"
    )
    str_args.add_argument(
        "--separate", action="store_true", 
        help="Use to separate EFM profiles in the 'filter' step."
    )

    snp_args = p.add_argument_group("Setting Specific for SNP Workflow")
    snp_args.add_argument(
        "--snp-type", default="all", dest="snptype",
        help="Specify the type of SNPs to include in the final report. 'p' will include only the "
        "Phenotype SNPs; 'a' will include only the Ancestry SNPs; 'i' will include only the "
        "Identity SNPs; and 'all' will include all SNPs. More than one type can be specified (e.g. "
        " 'p, a'). Default is all."
    )
    snp_args.add_argument(
        "--snps", action="store_true",
        help="Use to create a config file for the SNP workflow"
    )
    snp_args.add_argument(
        "--kintelligence", action="store_true",
        help="Use if processing Kintelligence SNPs within a Kintellience Report(s)"
    )
    snp_args.add_argument(
        "--snp-reference", dest="ref",
        help="Specify any references for SNP data for use in EFM."
    )
