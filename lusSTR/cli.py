#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import argparse
import csv
import os
import sys
import lusSTR
from lusSTR.annot import str_dict
from lusSTR.annot import split_string, rev_complement_anno, rev_comp_forward_strand_bracket
from lusSTR.annot import rev_comp_uas_output_bracket, loci_need_split_anno, traditional_str_allele
from lusSTR.annot import lus_anno, D21_bracket, TH01_annotation, PentaD_annotation


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--version', action='version', version='lusSTR v' + lusSTR.__version__
    )
    parser.add_argument(
        '-o', '--out', metavar='FILE',
        help='file to which output will be written; default is terminal (stdout)'
    )
    parser.add_argument(
        'input', help='sample(s) in CSV format; first four columns must be Locus, NumReads, '
        'Sequence, SampleID'
    )
    return parser


def main(args=None):
    if args is None:  # pragma: no cover
        args = get_parser().parse_args()

    cannot_split = [
        "D19S433", "D6S1043", "TH01", "D21S11", "D1S1656", "D7S820", "D5S818", "D12S391",
        "D9S1122", "D1S1656"
    ]
    must_split = ["D13S317", "D18S51"]

    sample_file = open(args.input, "r")
    data = csv.reader(sample_file)
    next(data)
    for row in data:
        locus = row[0]
        reads = row[1]
        sequence = row[2]
        sampleid = row[3]
        try:
            project = row[4]
            analysis = row[5]
        except IndexError:
            project = "NA"
            analysis = "NA"
        repeats = str_dict[locus]['Repeats']
        no_of_repeat_bases = len(str_dict[locus]['LUS'])
        no_of_sub_bases = str_dict[locus]['BasesToSubtract']
        no_of_split_bases = str_dict[locus]['NumBasesToSeparate']
        lus = str_dict[locus]['LUS']
        sec = str_dict[locus]['Sec']
        tert = str_dict[locus]['Tert']
        str_allele = traditional_str_allele(sequence, no_of_repeat_bases, no_of_sub_bases)
        if (
            locus in cannot_split or
            ((len(sequence) % no_of_repeat_bases != 0) and locus not in must_split)
           ):
            if str_dict[locus]['ReverseCompNeeded'] == "Yes":
                reverse_comp_sequence = rev_complement_anno(sequence)
                forward_strand_bracketed_form = rev_comp_forward_strand_bracket(
                    reverse_comp_sequence, no_of_repeat_bases, repeats, locus, cannot_split
                )
                reverse_strand_bracketed_form = rev_comp_uas_output_bracket(
                    forward_strand_bracketed_form, no_of_repeat_bases
                )
            elif locus == "D21S11":
                forward_strand_bracketed_form = D21_bracket(
                    sequence, no_of_split_bases, repeats
                )
            elif locus == "TH01" and (len(sequence) % no_of_repeat_bases != 0):
                forward_strand_bracketed_form = TH01_annotation(sequence, repeats)
            elif locus == "PentaD":
                forward_strand_bracketed_form = PentaD_annotation(
                    sequence, no_of_repeat_bases, repeats
                )
            else:
                forward_strand_bracketed_form = split_string(
                    sequence, no_of_repeat_bases, repeats
                )
            lus_final, sec_final, tert_final = lus_anno(
                forward_strand_bracketed_form, lus, sec, tert, locus, str_allele
            )
        else:
            if locus == "D18S51":
                if type(str_allele) == str:
                    forward_strand_bracketed_form = split_string(
                        sequence, no_of_repeat_bases, repeats
                    )
                else:
                    forward_strand_bracketed_form = loci_need_split_anno(
                        sequence, no_of_repeat_bases
                    )
            elif str_dict[locus]['ReverseCompNeeded'] == "Yes":
                reverse_comp_sequence = rev_complement_anno(sequence)
                forward_strand_bracketed_form = rev_comp_forward_strand_bracket(
                    reverse_comp_sequence, no_of_repeat_bases, repeats, locus, cannot_split
                )
                reverse_strand_bracketed_form = rev_comp_uas_output_bracket(
                    forward_strand_bracketed_form, no_of_repeat_bases
                )
            elif locus == "PentaD":
                forward_strand_bracketed_form = PentaD_annotation(
                    sequence, no_of_repeat_bases, repeats
                )
            else:
                forward_strand_bracketed_form = loci_need_split_anno(sequence, no_of_repeat_bases)
            lus_final, sec_final, tert_final = lus_anno(
                forward_strand_bracketed_form, lus, sec, tert, locus, str_allele
            )

        if locus == "PentaD":
            if str_allele == "2.2":
                lus_final = 5
            elif str_allele == "3.2":
                lus_final = 6
        lus_final_output = f"{str_allele}_{lus_final}"
        if sec_final == "":
            lus_plus = lus_final_output
        else:
            if tert_final == "":
                lus_plus = f"{str_allele}_{lus_final}_{sec_final}"
            else:
                lus_plus = f"{str_allele}_{lus_final}_{sec_final}_{tert_final}"
        if str_dict[locus]['ReverseCompNeeded'] == "Yes":
            summary = [
                sampleid, project, analysis, locus, sequence, reverse_comp_sequence, str_allele,
                forward_strand_bracketed_form, reverse_strand_bracketed_form, lus_final_output,
                lus_plus, reads
            ]
            summary = '\t'.join(str(i) for i in summary)
        else:
            summary = [
                sampleid, project, analysis, locus, sequence, sequence, str_allele,
                forward_strand_bracketed_form, forward_strand_bracketed_form, lus_final_output,
                lus_plus, reads
            ]
            summary = '\t'.join(str(i) for i in summary)

        output_file = sys.stdout
        if args.out is not None:
            if os.path.exists(args.out):
                output_file = open(args.out, "a")
            else:
                output_file = open(args.out, "w")
                header = [
                    'SampleID', 'Project', 'Analysis', 'Locus', 'UAS_Output_Sequence',
                    'Forward_Strand_Sequence', 'Traditional_STR_Allele',
                    'Forward_Strand_Bracketed_form', 'UAS_Output_Bracketed_Form', 'LUS',
                    'LUS_Plus', 'Reads'
                ]
                header = '\t'.join(header)
                print(header, file=output_file)
        print(summary, file=output_file)
        if args.out is not None:
            output_file.close()
