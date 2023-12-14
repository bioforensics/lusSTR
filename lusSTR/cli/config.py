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
import importlib.resources
import lusSTR
import os
from pathlib import Path
import yaml


def main(args):
    Path(args.workdir).mkdir(parents=True, exist_ok=True)
    if args.snps:
        final_dest = f"{args.workdir}/snp_config.yaml"
        config = importlib.resources.files("lusSTR") / "data/snp_config.yaml"
        final_config = edit_snp_config(config, args)
    else:
        final_dest = f"{args.workdir}/config.yaml"
        config = importlib.resources.files("lusSTR") / "data/config.yaml"
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
    if args.datatype:
        data["data_type"] = args.datatype
    if args.software:
        data["output_type"] = args.software
    if args.strand:
        data["strand"] = args.strand
    return data


def subparser(subparsers):
    p = subparsers.add_parser("config", description="Create config file for running STR pipeline")
    p.add_argument(
        "-w", "--workdir", metavar="W", default=".",
        help="directory to add config file; default is current working directory")
    p.add_argument(
        "--straitrazor", action="store_true",
        help="Use if sequences have been previously run through STRait Razor."
    )
    p.add_argument("--input", help="Input file or directory")
    p.add_argument("--out", "-o", help="Output file/directory name")
    p.add_argument(
        "--powerseq", action="store_true",
        help="Use to indicate sequences were created using the PowerSeq Kit."
    )
    p.add_argument(
        "--sex", action="store_true",
        help="Use if including the X and Y STR markers. Separate reports for these markers "
        "will be created.",
    )
    p.add_argument(
        "--nocombine", action="store_true",
        help="Do not combine read counts for duplicate sequences within the UAS region "
        "during the 'convert' step. By default, read counts are combined for sequences "
        "not run through the UAS.",
    )
    p.add_argument(
        "--reference", action="store_true", 
        help="Use for creating Reference profiles for STR workflow"
    )
    p.add_argument(
        "--software", choices=["efm", "mpsproto", "strmix"], default="strmix",
        help="Specify the probabilistic genotyping software package of choice. The final output"
        " files will be in the correct format for direct use. Default is strmix."
    )
    p.add_argument(
        "--str-type", choices=["ce", "ngs", "lusplus"], default="ngs",
        dest="datatype", help="Data type for STRs. Options are: CE allele ('ce'), sequence "
        "or bracketed sequence form('ngs'), or LUS+ allele ('lusplus'). Default is 'ngs'.",
    )
    p.add_argument(
        "--noinfo", action="store_true", 
        help="Use to not create the Sequence Information File in the 'filter' step"
    )
    p.add_argument(
        "--separate", action="store_true", 
        help="Use to separate EFM profiles in the 'filter' step. If specifying for SNPs, "
        "each sample will also be separated into 10 different bins for mixture deconvolution."
    )
    p.add_argument(
        "--nofiltering", action="store_true", 
        help="For STRs, use to perform no filtering during the 'filter' step. For SNPs, "
        "only alleles specified as 'Typed' by the UAS will be included at the 'format' step."
    )
    p.add_argument(
        "--snps", action="store_true",
        help="Use to create a config file for the SNP workflow"
    )
    p.add_argument(
        "--snp-type", default="all", dest="snptype",
        help="Specify the type of SNPs to include in the final report. 'p' will include only the "
        "Phenotype SNPs; 'a' will include only the Ancestry SNPs; 'i' will include only the "
        "Identity SNPs; and 'all' will include all SNPs. More than one type can be specified (e.g. "
        " 'p, a'). Default is all."
    )
    p.add_argument(
        "--kintelligence", action="store_true",
        help="Use if processing Kintelligence SNPs within a Kintellience Report(s)"
    )
    p.add_argument(
        "--snp-reference", dest="ref",
        help="Specify any references for SNP data for use in EFM."
    )
    p.add_argument(
        "--strand", choices=["uas", "forward"],
        help="Specify the strand orientation for the final output files. UAS orientation is "
        "default for STRs; forward strand is default for SNPs."
    )
