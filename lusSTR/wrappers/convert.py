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

import csv
import json
import numpy as np
import os
import pandas as pd
import re
import lusSTR
from lusSTR.scripts.marker import get_str_metadata_file, STRMarkerObject
from lusSTR.scripts.repeat import collapse_all_repeats, collapse_repeats_by_length
from lusSTR.scripts.repeat import sequence_to_bracketed_form, split_by_n
from lusSTR.scripts.repeat import reverse_complement, reverse_complement_bracketed
from pathlib import Path
from warnings import warn


with open(get_str_metadata_file(), "r") as fh:
    str_marker_data = json.load(fh)


def format_table(input, software, kit="forenseq", custom=False):
    """
    Function to format final output table and the flanking report (if necessary).
    """
    data = pd.read_csv(input, keep_default_na=False)
    data.iloc[:, 3] = data.iloc[:, 3].astype(str)
    list_of_lists = []
    flanks_list = []
    check_sr = 0
    for i, row in data.iterrows():
        locus = data.iloc[i, 0].upper()
        reads = data.iloc[i, 1]
        sequence = data.iloc[i, 2]
        sampleid = re.sub(" ", "_", data.iloc[i, 3])
        try:
            project = re.sub(" ", "_", data.iloc[i, 4])
            analysis = re.sub(" ", "_", data.iloc[i, 5])
        except IndexError:
            project = "NA"
            analysis = "NA"
        except TypeError:
            project = data.iloc[i, 4]
            analysis = data.iloc[i, 5]
        if locus == "PENTAD" or locus == "PENTA_D":
            locus = "PENTA D"
        if locus == "PENTAE" or locus == "PENTA_E":
            locus = "PENTA E"
        if locus == "DYS385A/B" or locus == "DYS385":
            locus = "DYS385A-B"
        if locus == "AMELOGENIN":
            continue
        metadata = str_marker_data[locus]
        if kit == "forenseq":
            remove_5p = metadata["Foren_5"]
            remove_3p = metadata["Foren_3"]
        else:
            remove_5p = metadata["Power_5"]
            remove_3p = metadata["Power_3"]
        if len(sequence) <= (remove_5p + remove_3p) and software != "uas":
            flank_summary = [
                sampleid,
                project,
                analysis,
                locus,
                reads,
                "NA",
                sequence,
                "NA",
                "NA",
                "NA",
                "Partial sequence",
            ]
            flanks_list.append(flank_summary)
            continue
        elif "N" in sequence:
            flank_summary = [
                sampleid,
                project,
                analysis,
                locus,
                reads,
                "NA",
                sequence,
                "NA",
                "NA",
                "NA",
                "Sequence contains Ns",
            ]
            flanks_list.append(flank_summary)
            continue
        marker = STRMarkerObject(locus, sequence, software, custom=custom, kit=kit)
        if locus == "D12S391" and kit == "powerseq" and software == "straitrazor":
            if "." in str(marker.canonical):
                check_sr += 1
                if check_sr > 10:
                    msg = (
                        "Multiple microvariants identified at D12 locus. "
                        "Please check STRait Razor version!!"
                    )
                    warn(msg)
        indel_flag = marker.indel_flag
        if indel_flag == "Possible indel or partial sequence":
            if locus == "PENTA D" and kit == "powerseq":
                marker = check_pentad(marker, sequence, software, custom)
            elif locus == "D7S820" and kit == "powerseq":
                marker = check_D7(marker, sequence, software, custom)
            elif locus == "VWA" and kit == "powerseq":
                marker = check_vwa(marker, sequence, software, custom)
        summary = [sampleid, project, analysis, locus] + marker.summary + [reads]
        list_of_lists.append(summary)
        if software != "uas":
            flank_summary = [
                sampleid,
                project,
                analysis,
                locus,
                reads,
                marker.canonical,
                marker.sequence,
                marker.flank_5p,
                marker.convert,
                marker.flank_3p,
                indel_flag,
            ]
            flanks_list.append(flank_summary)
    columns = [
        "SampleID",
        "Project",
        "Analysis",
        "Locus",
        "UAS_Output_Sequence",
        "Forward_Strand_Sequence",
        "Custom_Range_Sequence",
        "UAS_Output_Bracketed_Notation",
        "Forward_Strand_Bracketed_Notation",
        "Custom_Bracketed_Notation",
        "CE_Allele",
        "LUS",
        "LUS_Plus",
        "Reads",
    ]
    if not list_of_lists:
        final_output = pd.DataFrame(list_of_lists, columns=columns)
    else:
        final_output = sort_table(pd.DataFrame(list_of_lists, columns=columns))
    if software != "uas":
        flanks_columns = [
            "SampleID",
            "Project",
            "Analysis",
            "Locus",
            "Reads",
            "CE_Allele",
            "Full_Sequence",
            "5_Flank_Bracketed_Notation",
            "UAS_Region_Bracketed_Notation",
            "3_Flank_Bracketed_Notation",
            "Potential_Issues",
        ]
        if not flanks_list:
            final_flank_output = pd.DataFrame(flanks_list, columns=flanks_columns)
        else:
            final_flank_output = sort_table(pd.DataFrame(flanks_list, columns=flanks_columns))
    else:
        final_flank_output = ""
    if not custom:
        final_output = final_output.drop(
            ["Custom_Range_Sequence", "Custom_Bracketed_Notation"], axis=1
        )
    return final_output, final_flank_output, columns


def check_pentad(marker, sequence, software, custom):
    new_marker = marker
    if marker.summary[1][:4] == "AAAG" and marker.flank_5p[-5:] == "GAAAA":
        new_sequence = f"{sequence[:66]}-{sequence[66:]}"
        new_marker = STRMarkerObject(
            "PENTA D", new_sequence, software, custom=custom, kit="powerseq"
        )
    elif marker.summary[1][-4:] == "AAAG" and marker.flank_3p[:7] == "AAAAA A":
        new_sequence = f"{sequence}-"
        new_marker = STRMarkerObject(
            "PENTA D", new_sequence, software, custom=custom, kit="powerseq"
        )
    else:
        return marker
    return new_marker


def check_D7(marker, sequence, software, custom):
    if marker.summary[1][:3] == "AAC" and marker.flank_5p[-8:] == "T AAAAAA":
        new_sequence = f"{sequence[:60]}-{sequence[60:]}"
        new_marker = STRMarkerObject(
            "D7S820", new_sequence, software, custom=custom, kit="powerseq"
        )
    else:
        return marker
    return new_marker


def check_vwa(marker, sequence, software, custom):
    if sequence[:4] == "GATA":
        new_sequence = f"{sequence[:1]}-{sequence[1:]}"
        new_marker = STRMarkerObject("VWA", new_sequence, software, custom=custom, kit="powerseq")
    else:
        return marker
    return new_marker


def combine_reads(table, columns):
    comb_table = table.groupby(columns[:-1], as_index=False)["Reads"].sum()
    sorted = sort_table(comb_table)
    return sorted


def sort_table(table):
    sorted_table = table.sort_values(
        by=["SampleID", "Project", "Analysis", "Locus", "Reads", "CE_Allele"], ascending=False
    )
    return sorted_table


def main(input, out, kit, software, sex, nocombine, custom):
    input = str(input)
    out = str(out)
    output_name = os.path.splitext(out)[0]
    input_name = os.path.splitext(input)[0]
    autosomal_final_table, autosomal_flank_table, columns = format_table(
        input, software, kit, custom
    )
    if sex:
        sex_final_table, sex_flank_table, columns = format_table(
            f"{input_name}_sexloci.csv", software, kit, custom
        )
        if software != "uas":
            if not sex_final_table.empty:
                sex_flank_table.to_csv(f"{output_name}_sexloci_flanks.txt", sep="\t", index=False)
                if nocombine:
                    sex_final_table.to_csv(
                        f"{output_name}_sexloci_no_combined_reads.txt", index=False
                    )
                sex_final_table = combine_reads(sex_final_table, columns, custom)
                sex_final_table.to_csv(f"{output_name}_sexloci.txt", sep="\t", index=False)
        else:
            sex_final_table.to_csv(f"{output_name}_sexloci.txt", sep="\t", index=False)
    if software != "uas":
        if not autosomal_final_table.empty:
            autosomal_flank_table.to_csv(f"{output_name}_flanks.txt", sep="\t", index=False)
            if nocombine:
                autosomal_final_table.to_csv(
                    f"{output_name}_no_combined_reads.txt", sep="\t", index=False
                )
            if custom:
                custom_columns = [
                    "SampleID",
                    "Project",
                    "Analysis",
                    "Locus",
                    "Custom_Range_Sequence",
                    "Custom_Bracketed_Notation",
                    "CE_Allele",
                    "LUS",
                    "LUS_Plus",
                    "Reads",
                ]
                custom_table = autosomal_final_table[custom_columns]
                custom_table_comb = combine_reads(custom_table, custom_columns)
                custom_table_comb.to_csv(out, sep="\t", index=False)
            autosomal_final_table = combine_reads(autosomal_final_table, columns)
            full_table_name = re.sub(r"_custom_range", "", output_name)
            autosomal_final_table.to_csv(f"{full_table_name}.txt", sep="\t", index=False)
    else:
        autosomal_final_table.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    main(
        snakemake.input,
        snakemake.output,
        kit=snakemake.params.kit,
        software=snakemake.params.a_software,
        sex=snakemake.params.sex,
        nocombine=snakemake.params.nocombine,
        custom=snakemake.params.custom,
    )
