#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
from pkg_resources import resource_filename
import re


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

    Sequences in which the UAS software output contains the sequence on the reverse strand require
    translation of the sequence to the forward strand. This allows for consistency between both
    loci and any outside analyses in which comparisons may be made.
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

    Function is used to create the bracketed annotation for reverse complement sequences (i.e. the
    forward strand). It calls additional functions depending on the locus and/or if the sequence
    is a microvariant or not.
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