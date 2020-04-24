#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import csv
import json
import os
from pkg_resources import resource_filename
import re
import sys
import lusSTR


def get_str_metadata_file():
    return resource_filename('lusSTR', 'str_markers.json')


with open(get_str_metadata_file(), 'r') as fh:
    str_dict = json.load(fh)


def get_annotation(sequence, repeat_list):
    '''
    Function used to (in part) create bracketed annotation.

    This function, called as a part of the split_string() function, identifies and extracts
    specified repeats in order defined in the dict. The function returns the repeats in the
    bracketed form with any remaining sequence returned intact. This function is used to
    identify the specified repeats in the prioritized order within the sequence.
    '''
    for k in range(len(repeat_list)):
        repeat = repeat_list[k]
        if k == 0:
            count = 0
            final = list()
            for element in re.split(repeat, sequence):
                parts = element.split(',')
                for i in parts:
                    if i == "":
                        count += 1
                    else:
                        if count == 1:
                            final.append(repeat)
                        elif count >= 2:
                            final.append(f"[{repeat}]{count}")
                        count = 1
                        final.append(i)
            if parts[-1] == "" and count > 2:
                final.append(f"[{repeat}]{count-1}")
            elif parts[-1] == "" and count <= 2:
                final.append(repeat)
            tmp = ' '.join(final)
        elif k > 0:
            final = list()
            if re.match(repeat, sequence):
                if sequence[:4] == repeat:
                    count = 0
                else:
                    count = 1
            else:
                count = 0
            for element in re.split(repeat, tmp):
                parts = element.split(',')
                for i in parts:
                    if i == "":
                        count += 1
                    else:
                        if count == 1:
                            final.append(repeat)
                        elif count >= 2:
                            final.append(f"[{repeat}]{count}")
                        count = 1
                        final.append(i)
            if parts[-1] == "" and count > 2:
                final.append(f"[{repeat}]{count-1}")
            elif parts[-1] == "" and count <= 2:
                final.append(repeat)
            tmp = ' '.join(final)
    return re.sub("  ", " ", tmp)


def split_by_n(sequence, n):
    '''
    Function to divide sequence into chunks of n
    '''
    while sequence:
        yield sequence[:n]
        sequence = sequence[n:]


def split_string(sequence, n, repeat_list):
    '''
    Function to create bracketed annotation.

    Function creates bracketed annotation using a prioritized ordered list of repeats (repeat
    units are marker specific and are specified in str_markers.json file).
    '''
    strings = get_annotation(sequence, repeat_list)
    final_string = list()
    for unit in strings.split(' '):
        if len(unit) > n and "[" not in unit:
            for x in split_by_n(unit, n):
                final_string.append(x)
        else:
            final_string.append(unit)
    final_string_formatted = ' '.join(final_string)
    return re.sub("  ", " ", final_string_formatted)


def rev_complement_anno(sequence):
    '''
    Function creates reverse complement of sequence

    Sequences in which the UAS software output contains the sequence on the reverse strand
    require translation of the sequence to the forward strand. This allows for consistency
    between both loci and any outside analyses in which comparisons may be made.
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(sequence)
    bases = [complement[base] for base in bases]
    comp = ''.join(bases)
    reverse_comp_sequence = comp[::-1]
    return reverse_comp_sequence


def rev_comp_forward_strand_bracket(rev_sequence, n, repeat_list, locusid, cannot_split_list):
    '''
    Function creates bracketed annotation for reverse complement sequences

    Function is used to create the bracketed annotation for reverse complement sequences (i.e
    the forward strand). It calls additional functions depending on the locus and/or if the
    sequence is a microvariant or not.
    '''
    if locusid in cannot_split_list:
        if locusid == "D19S433":
            forward_strand_bracketed_form = D19_annotation(rev_sequence, repeat_list, "CCTT")
        elif locusid == "D1S1656":
            forward_strand_bracketed_form = D1_annotation(rev_sequence, repeat_list, "CACA")
        else:
            forward_strand_bracketed_form = split_string(rev_sequence, n, repeat_list)
    elif locusid == "FGA":
        if len(rev_sequence) % n != 0:
            forward_strand_bracketed_form = FGA_anno(rev_sequence, repeat_list)
        else:
            forward_strand_bracketed_form = loci_need_split_anno(rev_sequence, n)
    else:
        forward_strand_bracketed_form = loci_need_split_anno(rev_sequence, n)
    return re.sub("  ", " ", forward_strand_bracketed_form)


def rev_comp_uas_output_bracket(forward_bracket, n):
    '''
    Function creates bracketed annotation for the UAS output sequence.

    There may be instances where the bracketed annotation is required for the UAS output sequence.
    This function performs the reverse complement of the forward strand bracketed annotation to
    create the bracketed annotation for the UAS output (for those loci in which the UAS output
    reports the reverse strand).
    '''
    ind = list(forward_bracket)
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_strand_form = list()
    for j in ind:
        if j.isalpha():
            reverse_strand_form.append(complement[j])
        elif j == "[":
            reverse_strand_form.append("]")
        elif j == "]":
            reverse_strand_form.append("[")
        else:
            reverse_strand_form.append(j)
    reverse_form_anno_tmp = ''.join(reversed(reverse_strand_form))
    reverse_form_anno_final = list()
    for unit in reverse_form_anno_tmp.split(' '):
        if "[" in unit:
            if len(unit) == (n+3):
                final_string = f"{unit[1:(len(unit))]}{unit[0]}"
            else:
                final_string = f"{unit[2:(len(unit))]}{unit[1]}{unit[0]}"
            reverse_form_anno_final.append(final_string)
        else:
            reverse_form_anno_final.append(unit)
    reverse_strand_bracketed_form = ' '.join(reverse_form_anno_final)
    return re.sub("  ", " ", reverse_strand_bracketed_form)


def get_blocks(sequence, n):
    '''
    Function to split a sequence into blocks of size n

    This function is used as a part of the loci_need_split_anno() function. It splits the sequence
    into blocks of size n bases (as specified in the str_markers.json file).
    '''
    count = 0
    prev = None
    for unit in split_by_n(sequence, n):
        if unit != prev:
            if prev is not None:
                yield prev, count
            prev = unit
            count = 0
        count += 1
    yield prev, count


def loci_need_split_anno(sequence, n):
    '''
    Function creates bracketed annotation form using the sequence split into blocks of size n.

    This is the simplest way to create bracketed annotation form. It splits sequence by specified
    # of bases then checks for repeats of units in sequential order.
    '''
    alleles = list()
    for unit, count in get_blocks(sequence, n):
        if count == 1:
            alleles.append(unit)
        else:
            alleles.append(f"[{unit}]{count}")
    alleles_final = '  '.join(alleles)
    return re.sub("  ", " ", alleles_final)


def traditional_str_allele(sequence, n, n_sub_out):
    '''
    Function used to calculate the traditional STR allele designation
    '''
    new_seq = sequence[:(len(sequence)-n_sub_out)]
    if (len(new_seq) % n) == 0:
        trad_allele = int(len(new_seq)/n)
    else:
        allele_tmp = int(len(new_seq)/n)
        allele_dec = int(len(new_seq) % n)
        trad_allele = f"{allele_tmp}.{allele_dec}"
    return trad_allele


def extract(s, single_repeat):
    '''
    Function to identify the # of alleles for a specified repeat unit from the bracketed
    annotation form.
    '''
    finalcount = 0
    for m in re.finditer(single_repeat, s):
        count = s[m.end()+1:m.end()+3]
        if count == "" or count[0] == "[" or count[0] == " " or count.isalpha():
            count = 1
        try:
            if float(count) > float(finalcount):
                finalcount = count
                try:
                    if str(finalcount)[1] == " ":
                        finalcount = finalcount[0]
                except IndexError:
                    count = count
        except ValueError:
            count = count
    return finalcount


def lus_anno(sequence, lus, sec, tert, locusid, str_allele):
    '''
    Function to identify LUS alleles, secondary motif alleles (if applicable) and tertiary motif
    alleles (if applicable).

    This function identifies the specified repeat unit (from the STR_markers.json file) and finds
    it in the forward strand bracketed annotation. These alleles are used in the LUS+ annotation.

    The D7S820 locus requires a different script due to it's designated tertiary motif (+T at end
    of sequence).

    The D21S11 locus also requires a separate function because both the lus and secondary motif
    are the same, but differ based on the location of the repeat.
    '''
    if locusid == "D21S11":
        lus_allele, sec_allele, tert_allele = lus_D21_anno(sequence, lus, sec, tert)
    else:
        lus_allele = extract(sequence, lus)
        if sec != "":
            sec_allele = extract(sequence, sec)
            if tert != "":
                tert_allele = extract(sequence, tert)
            elif locusid == "D7S820":
                if str(sequence)[-1] == "T" and isinstance(str_allele, str):
                    tert_allele = 1
                else:
                    tert_allele = 0
            else:
                tert_allele = ""
        else:
            sec_allele = ""
            tert_allele = ""
    return lus_allele, sec_allele, tert_allele


def D21_lus_sec(sequence, repeat, tert):
    '''
    Function to identify the number of LUS and secondary motif alleles for the D21S11 locus.

    A separate function is required because the LUS repeat motif is the last "TCTA" repeat set and
    the secondary repeat motif is the first set of "TCTA" repeats in the sequence.
    '''
    remaining = list()
    lus_sec = list()
    lus_allele = None
    for element in re.split(tert, sequence):
        if element == sequence:
            lus_sec = extract(element, repeat)
            sec_allele = lus_sec[0]
            if len(lus_sec) == 1 and element[-4:] == repeat:
                lus_allele = 1
            elif len(lus_sec) == 1 and element[-4:] != repeat:
                lus_allele = 0
            else:
                lus_allele = lus_sec[1]
        else:
            parts = element.split('[,]')
            for i in parts:
                if i != "":
                    repeats = extract(i, repeat)
                    lus_sec.append(repeats)
    if lus_allele is None:
        lus_allele = lus_sec[1]
        sec_allele = lus_sec[0]
    return lus_allele, sec_allele


def extract_D21_tert(s, single_repeat):
    '''
    Function to identify the number of tertiary motif alleles for the D21S11 locus.

    A separate function from the LUS and secondary motif is required to identify the tertiary
    motif for this locus.
    '''
    finalcount = 0
    for m in re.finditer(single_repeat, s):
        count = s[m.end()+1:m.end()+3]
        if count == "" or count[0] == "[" or count[0] == " " or count.isalpha():
            count = 1
        try:
            if float(count) > float(finalcount):
                finalcount = count
                try:
                    if str(finalcount)[1] == " ":
                        finalcount = finalcount[0]
                except IndexError:
                    count = count
        except ValueError:
            count = count
    return finalcount


def lus_D21_anno(sequence, lus, sec, tert):
    '''
    Function to compile the number of LUS alleles, secondary motif alleles and tertiary motif
    alleles for the D21S11 locus.
    '''
    lus_allele, sec_allele = D21_lus_sec(sequence, lus, tert)
    tert_allele = extract_D21_tert(sequence, tert)
    return lus_allele, sec_allele, tert_allele


def D21_bracket(sequence, no_of_split_bases, repeats):
    '''
    Function to create bracketed annotation for the D21S11 locus.

    A specialized function is required for this locus due to the potential end of the sequence
    containing 'TA TCTA' and other variants. This sequence needs to remain intact to conform with
    the conventional annotation for this particular locus. However, if the 'TATCTA' is included in
    a repeat unit, the repeat unit needs to be reported (i.e. [TCTA]2).
    '''
    forward_strand_bracketed_form = split_string(sequence, no_of_split_bases, repeats)
    prev = 0
    for m in re.finditer("]", forward_strand_bracketed_form):
        prev = m.end()
    if (
        prev == (len(forward_strand_bracketed_form) - 1) or
        prev == (len(forward_strand_bracketed_form) - 2)
       ):
        return forward_strand_bracketed_form
    else:
        first_string = forward_strand_bracketed_form[:prev+2]
        second_string = forward_strand_bracketed_form[prev+2:]
        second_string_final = re.sub(" ", "", second_string)
        if len(second_string_final) % 4 == 0:
            split_second_string = loci_need_split_anno(second_string_final, 4)
            final_string = f"{first_string} {second_string}"
        elif len(second_string_final) == 6:
            third_string = second_string_final[-6:-4]
            fourth_string = second_string_final[-4:]
            final_string = f"{first_string} {third_string} {fourth_string}"
        elif len(second_string_final) % 4 == 2:
            third_string = second_string_final[:-6]
            fourth_string = second_string_final[-6:-4]
            last_string = second_string_final[-4:]
            third_string_final = loci_need_split_anno(third_string, 4)
            final_string = f"{first_string} {third_string_final} {fourth_string} {last_string}"
        else:
            third_string = loci_need_split_anno(second_string_final, 4)
            final_string = f"{first_string} {third_string}"
        return re.sub("  ", " ", final_string)


def split_sequence_into_two_strings(sequence, repeat_for_split):
    '''
    Function to split a sequence into two separate strings at a specified repeat unit.
    '''
    last = 0
    prev = 0
    for m in re.finditer(repeat_for_split, sequence):
        if m.start() == prev or m.start() == last or prev == 0:
            prev = m.end()
        else:
            last = m.end()
    first_string = sequence[:prev]
    second_string = sequence[prev:]
    return first_string, second_string


def D1_annotation(sequence, repeat_list, repeat_for_split):
    '''
    Function to create bracketed annotation for the D1S1656 locus.

    This function identifies if the sequence is a microvariant in order to call different
    functions to create the bracketed annotation.
    '''
    sequence_filt = sequence[2:]
    final = list()
    first_string, second_string = split_sequence_into_two_strings(sequence_filt, repeat_for_split)
    final.append(sequence[:2])
    if first_string == "":
        final.append(repeat_for_split)
    else:
        final.append(loci_need_split_anno(first_string, 4))
    if (len(second_string) % 4 != 0):
        final.append(split_string(second_string, 4, repeat_list))
    else:
        final.append(loci_need_split_anno(second_string, 4))
    final_string = ' '.join(final)
    return re.sub("  ", " ", final_string)


def D19_annotation(sequence, repeat_list, repeat_for_split):
    '''
    Function to create bracketed annotation for the D19S433 locus

    A specialized function is required for this locus. The sequence is first broken into two
    different strings. The two sets of sequences are processed separately in order to identify the
    potential presence of a deletion in either sequence.

    Simply identifying repeat units in a specified order does not result in the final annotation
    which is consistent with previously published annotation for this locus.
    '''
    final = list()
    last = 0
    prev = 0
    for m in re.finditer(repeat_for_split, sequence):
        if m.start() == prev or m.start() == last:
            prev = m.end()
        else:
            last = m.end()
    final.append(sequence[:2])
    first_string = sequence[2:prev]
    second_string = sequence[prev:]
    if (len(first_string) % 4 != 0):
        final.append(split_string(first_string, 4, repeat_list))
    else:
        final.append(loci_need_split_anno(first_string, 4))
    if (len(second_string) % 4 != 0):
        third_string = second_string[:-6]
        final.append(loci_need_split_anno(third_string, 4))
        final.append(second_string[-6:-4])
        final.append(second_string[-4:])
    else:
        final.append(loci_need_split_anno(second_string, 4))
    final_string = ' '.join(final)
    return re.sub("  ", " ", final_string)


def FGA_anno(sequence, repeat_list):
    '''
    Function to create the forward strand bracketed annotation for FGA locus.

    A specialized function is required because which repeat unit should be identified differs
    based on its location in the sequence. For example, the "GGAA" repeat should be identified at
    the beginning of the sequence; the "GAAA" repeat should be identified at the end of the
    sequence; and the repeat "AAAG" should be identified within the two end repeats.

    Simply identifying repeat units in a specified order does not result in the final annotation
    which is consistent with previously published annotation for this locus.
    '''
    final = list()
    prev = 0
    if (len(sequence) % 4 == 0):
        final_string = loci_need_split_anno(sequence, 4)
    else:
        for m in re.finditer("GGAA", sequence):
            if prev == 0 or m.start() == prev:
                prev = m.end()
            else:
                break
        first_string = sequence[:prev]
        second_string = sequence[prev:]
        prev = 0
        for m in re.finditer("AAAA", second_string):
            prev = m.start()
            break
        if second_string[prev:(prev+6)] == "AAAAAA":
            third_string = second_string[:prev+2]
            fourth_string = second_string[prev+2:]
        elif prev == 0:
            third_string = second_string[:-6]
            fourth_string = second_string[-6:]
        else:
            third_string = second_string[:prev]
            fourth_string = second_string[prev:]
        final.append(loci_need_split_anno(first_string, 4))
        final.append(split_string(third_string, 4, repeat_list))
        count = 0
        tmp = list()
        for element in re.split("GAAA", fourth_string):
            parts = element.split(',')
            for i in parts:
                if i == "":
                    count += 1
                else:
                    if count == 1:
                        tmp.append("GAAA")
                    elif count >= 2:
                        tmp.append("[GAAA]" + str(count))
                    count = 1
                    if i == "AAAAAA":
                        tmp.append("AA AAAA")
                    elif len(i) > 4:
                        for x in split_by_n(i, 4):
                            tmp.append(x)
                    else:
                        tmp.append(i)
        if parts[-1] == "" and count > 2:
            tmp.append("[GAAA]" + str(count-1))
        elif parts[-1] == "" and count <= 2:
            tmp.append("GAAA")
        last_string_final = ' '.join(tmp)
        final.append(last_string_final)
        final_string = ' '.join(final)
        return re.sub("  ", " ", final_string)


def TH01_annotation(sequence, repeat_list):
    '''
    Function to create bracketed annotation for the TH01 locus.

    A separate function is required for the microvariants of the TH01 locus because of the
    insertion of the "ATG" between the repeat units "AATG".
    '''
    strings = get_annotation(sequence, repeat_list)
    final_string = list()
    for unit in strings.split(' '):
        if "[" not in unit and len(unit) > 3 and (len(unit) % 4 != 0) and unit[:3] == "ATG":
            group1 = unit[:3]
            final_string.append(group1)
            for x in split_by_n(unit[3:], n=4):
                final_string.append(x)
        else:
            final_string.append(unit)
    final_form = ' '.join(final_string)
    return final_form


def PentaD_annotation(sequence, no_of_repeat_bases, repeat_list):
    '''
    Function to create bracketed annotation for the PentaD locus.

    A separate function is required for the PentaD locus in order to first identify smaller
    sequences (i.e. those with length <18bp) to annotate the bracketed form correctly.
    If the sequence is >= 18bp, the flanking region (first 5 bases) is first split off in the
    sequence to preserve that sequence. Then the repeat units are identified and bracketed.
    '''
    if len(sequence) < 18:
        final_string = split_string(sequence, no_of_repeat_bases, repeat_list)
        return final_string
    else:
        first_string = sequence[:5]
        second_string = sequence[5:]
        second_string_anno = split_string(second_string, no_of_repeat_bases, repeat_list)
        final_string = f"{first_string} {second_string_anno}"
        return re.sub("  ", " ", final_string)


def full_seq_to_uas(full_seq, front, back):
    '''
    Function to trim full sequences to the UAS region.

    It identifies the number of bases to remove from the 5' and 3' ends of the sequence to
    trim to the UAS region. The downstream annotation, including length-based allele
    designations, LUS, LUS+ and bracketed annotation is based on this region in the sequence.
    '''
    if front == 0:
        seq_uas = full_seq[:-back]
    elif back == 0:
        seq_uas = full_seq[front:]
    else:
        seq_uas = full_seq[front:-back]
    return seq_uas


def main(args):
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
        foren_5 = str_dict[locus]['Foren_5']
        foren_3 = str_dict[locus]['Foren_3']
        power_5 = str_dict[locus]['Power_5']
        power_3 = str_dict[locus]['Power_3']
        if args.uas:
            uas_sequence = sequence
        else:
            if args.kit == "forenseq":
                if str_dict[locus]['ReverseCompNeeded'] == "No":
                    uas_sequence = full_seq_to_uas(sequence, foren_5, foren_3)
                else:
                    uas_from_full = full_seq_to_uas(sequence, foren_5, foren_3)
                    uas_sequence = rev_complement_anno(uas_from_full)
            elif args.kit == "powerseq":
                if str_dict[locus]['ReverseCompNeeded'] == "No":
                    uas_sequence = full_seq_to_uas(sequence, power_5, power_3)
                else:
                    uas_from_full = full_seq_to_uas(sequence, power_5, power_3)
                    uas_sequence = rev_complement_anno(uas_from_full)
        str_allele = traditional_str_allele(uas_sequence, no_of_repeat_bases, no_of_sub_bases)
        cantsplit = locus in cannot_split
        havetosplit = locus in must_split
        split_incompatible = len(uas_sequence) % no_of_repeat_bases != 0 and not havetosplit
        if cantsplit or split_incompatible:
            if str_dict[locus]['ReverseCompNeeded'] == "Yes":
                reverse_comp_sequence = rev_complement_anno(uas_sequence)
                print(reverse_comp_sequence)
                forward_strand_bracketed_form = rev_comp_forward_strand_bracket(
                    reverse_comp_sequence, no_of_repeat_bases, repeats, locus, cannot_split
                )
                print(forward_strand_bracketed_form)
                reverse_strand_bracketed_form = rev_comp_uas_output_bracket(
                    forward_strand_bracketed_form, no_of_repeat_bases
                )
            elif locus == "D21S11":
                forward_strand_bracketed_form = D21_bracket(
                    uas_sequence, no_of_split_bases, repeats
                )
            elif locus == "TH01" and (len(uas_sequence) % no_of_repeat_bases != 0):
                forward_strand_bracketed_form = TH01_annotation(uas_sequence, repeats)
            elif locus == "PentaD":
                forward_strand_bracketed_form = PentaD_annotation(
                    uas_sequence, no_of_repeat_bases, repeats
                )
            else:
                forward_strand_bracketed_form = split_string(
                    uas_sequence, no_of_repeat_bases, repeats
                )
            lus_final, sec_final, tert_final = lus_anno(
                forward_strand_bracketed_form, lus, sec, tert, locus, str_allele
            )
        else:
            if locus == "D18S51":
                if type(str_allele) == str:
                    forward_strand_bracketed_form = split_string(
                        uas_sequence, no_of_repeat_bases, repeats
                    )
                else:
                    forward_strand_bracketed_form = loci_need_split_anno(
                        uas_sequence, no_of_repeat_bases
                    )
            elif str_dict[locus]['ReverseCompNeeded'] == "Yes":
                reverse_comp_sequence = rev_complement_anno(uas_sequence)
                forward_strand_bracketed_form = rev_comp_forward_strand_bracket(
                    reverse_comp_sequence, no_of_repeat_bases, repeats, locus, cannot_split
                )
                reverse_strand_bracketed_form = rev_comp_uas_output_bracket(
                    forward_strand_bracketed_form, no_of_repeat_bases
                )
            elif locus == "PentaD":
                forward_strand_bracketed_form = PentaD_annotation(
                    uas_sequence, no_of_repeat_bases, repeats
                )
            else:
                forward_strand_bracketed_form = loci_need_split_anno(
                    uas_sequence, no_of_repeat_bases
                )
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
                sampleid, project, analysis, locus, uas_sequence, reverse_comp_sequence,
                str_allele, forward_strand_bracketed_form, reverse_strand_bracketed_form,
                lus_final_output, lus_plus, reads
            ]
            summary = '\t'.join(str(i) for i in summary)
        else:
            summary = [
                sampleid, project, analysis, locus, uas_sequence, uas_sequence, str_allele,
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
