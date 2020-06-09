#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import json
import lusSTR
from lusSTR.annot import split_sequence_into_two_strings
from lusSTR.repeat import collapse_repeats_by_length, sequence_to_bracketed_form
from lusSTR.repeat import reverse_complement, reverse_complement_bracketed
from lusSTR.repeat import repeat_copy_number, collapse_all_repeats, split_by_n
from pkg_resources import resource_filename
import re


def get_str_metadata_file():
    return resource_filename('lusSTR', 'str_markers.json')


with open(get_str_metadata_file(), 'r') as fh:
    str_marker_data = json.load(fh)


class InvalidLocusError(ValueError):
    pass


class UnsupportedKitError(ValueError):
    pass


class STRMarker():
    def __init__(self, locus, sequence, uas=False, kit='forenseq'):
        self.locus = locus
        self.sequence = sequence
        if locus not in str_marker_data:
            raise InvalidLocusError(locus)
        self.data = str_marker_data[locus]
        self.uas = uas
        if kit.lower() not in ('forenseq', 'powerseq'):
            raise UnsupportedKitError(kit)
        self.kit = kit.lower()
        if uas and self.data['ReverseCompNeeded'] == "Yes":
            self.sequence = reverse_complement(sequence)

    @property
    def repeat_size(self):
        return len(self.data['LUS'])

    @property
    def repeats(self):
        return self.data['Repeats']

    def _uas_bases_to_trim(self):
        '''Number of bases to trim off each side to get the UAS sequence.

        ForenSeq and PowerSeq amplicons extend beyond the core UAS sequence for some loci. This
        function determines the number of bases that need to be trimmed from the full amplicon
        sequence to recover the UAS core sequence.
        '''
        if self.uas:
            return 0, 0
        elif self.kit == 'forenseq':
            return self.data['Foren_5'], self.data['Foren_3']
        elif self.kit == 'powerseq':
            return self.data['Power_5'], self.data['Power_3']
        else:
            raise UnsupportedKitError(self.kit)

    @property
    def forward_sequence(self):
        '''Sequence from the UAS region, in the forward orientation

        If it is the full ForenSeq or PowerSeq amplicon sequence, trim the full forward sequence
        back to the UAS region. If the sequence has already been run through UAS, no trimming is
        required.
        '''
        if self.uas:
            return self.sequence
        front, back = self._uas_bases_to_trim()
        if back == 0:
            back = None
        else:
            back *= -1
        return self.sequence[front:back]

    @property
    def uas_sequence(self):
        '''Sequence from the UAS region, matching the orientation of the UAS output

        The UAS software outputs the reverse complement of the forward sequence for some loci.
        '''
        if self.data['ReverseCompNeeded'] == "Yes":
            return reverse_complement(self.forward_sequence)
        return self.forward_sequence

    @property
    def flankseq_5p(self):
        if self.uas:
            return None
        front, back = self._uas_bases_to_trim()
        return self.sequence[:front]

    @property
    def flank_5p(self):
        if self.uas:
            return None
        flank_rev = collapse_repeats_by_length(self.flankseq_5p[::-1], self.repeat_size)
        flank = flank_rev[::-1]
        return flank

    @property
    def flankseq_3p(self):
        if self.uas:
            return None
        front, back = self._uas_bases_to_trim()
        if back == 0:
            return ''
        return self.sequence[-back:]

    @property
    def flank_3p(self):
        if self.uas:
            return None
        return collapse_repeats_by_length(self.flankseq_3p, self.repeat_size)

    @property
    def canonical(self):
        '''Canonical STR allele designation'''
        n = self.repeat_size
        nsubout = self.data['BasesToSubtract']
        if nsubout == 0:
            nsubout = None
        else:
            nsubout *= -1
        new_seq = self.uas_sequence[:nsubout]
        if len(new_seq) % n == 0:
            canon_allele = int(len(new_seq) / n)
        else:
            allele_int = int(len(new_seq) / n)
            allele_dec = int(len(new_seq) % n)
            canon_allele = f'{allele_int}.{allele_dec}'
        return canon_allele

    @property
    def cannot_split(self):
        return self.locus in [
            'D19S433', 'D6S1043', 'TH01', 'D21S11', 'D1S1656', 'D7S820', 'D5S818', 'D12S391',
            'D9S1122', 'PentaE'
        ]

    @property
    def must_split(self):
        return self.locus in ['D13S317', 'D18S51']

    @property
    def split_compatible(self):
        return self.must_split or len(self.uas_sequence) % self.repeat_size == 0

    @property
    def do_split(self):
        return not self.cannot_split or self.split_compatible

    @property
    def annotation(self):
        bylength = (
            (self.data['ReverseCompNeeded'] == 'Yes' and self.split_compatible)
            or (self.locus == 'D3S1358' and self.split_compatible)
            or self.locus == 'D16S539'
        )
        if bylength:
            collapseseq = collapse_repeats_by_length(self.forward_sequence, self.repeat_size)
        else:
            collapseseq = sequence_to_bracketed_form(
                self.forward_sequence, self.repeat_size, self.repeats
            )
        return collapseseq

    @property
    def annotation_uas(self):
        if self.data['ReverseCompNeeded'] == 'Yes':
            return reverse_complement_bracketed(self.annotation)
        return self.annotation

    @property
    def annotation_with_flanks(self):
        full_annot = f'{self.flank_5p} {self.annotation} {self.flank_3p}'
        return full_annot.strip()

    @property
    def designation(self):
        lus, sec, ter = None, None, None
        lus = repeat_copy_number(self.annotation, self.data['LUS'])
        if self.data['Sec'] != '':
            sec = repeat_copy_number(self.annotation, self.data['Sec'])
        if self.data['Tert'] != '':
            ter = repeat_copy_number(self.annotation, self.data['Tert'])
        return lus, sec, ter

    @property
    def summary(self):
        lus, sec, ter = self.designation
        canon = self.canonical
        lus_final_output = f'{canon}_{lus}'
        if sec is None:
            lus_plus = lus_final_output
        else:
            if ter is None:
                lus_plus = f'{canon}_{lus}_{sec}'
            else:
                lus_plus = f'{canon}_{lus}_{sec}_{ter}'
        return [
            self.uas_sequence, self.forward_sequence, canon, self.annotation,
            self.annotation_uas, lus_final_output, lus_plus
        ]


class STRMarker_D8S1179(STRMarker):
    @property
    def flank_5p(self):
        return ''


class STRMarker_D13S317(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return (
            f'{flank_seq[:2]} {collapse_repeats_by_length(flank_seq[2:14], 4)} {flank_seq[14]} '
            f'{flank_seq[15]} {flank_seq[16:19]} {collapse_repeats_by_length(flank_seq[19:], 4)}'
        )

    @property
    def annotation(self):
        if len(self.uas_sequence) < 110:
            bracketed_form = collapse_repeats_by_length(self.uas_sequence, 4)
        else:
            for m in re.finditer('GGGC', self.uas_sequence):
                break_point = m.end()
            bracketed_form = (
                f'{collapse_repeats_by_length(self.uas_sequence[:break_point], 4)} '
                f'{sequence_to_bracketed_form(self.uas_sequence[break_point:], 4, self.repeats)}'
            )
        return bracketed_form


class STRMarker_D20S482(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return (
            f'{flank_seq[:2]} {flank_seq[2:6]} {flank_seq[6]} {flank_seq[7:10]} '
            f'{flank_seq[10:]} {flank_seq[14:18]}'
        )


class STRMarker_D2S441(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return f'{flank_seq[:4]} {flank_seq[4]} {collapse_repeats_by_length(flank_seq[5:], 4)}'


class STRMarker_D7S820(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return f'{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:13], 4)} {flank_seq[13:]}'

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        return collapse_repeats_by_length(flank_seq, 4)

    @property
    def annotation(self):
        '''
        Function to correctly bracket microvariants in the D7S820 locus.

        Microvariants N.1 and N.3 separate the first base (or first three bases) of the sequence
        in order to correctly bracket the sequence. Microvariant N.2 has a third repeat motif
        added in order to bracket the sequence correctly.
        '''
        sequence = self.forward_sequence
        if type(self.canonical) == int:
            forward_strand_brack_form = sequence_to_bracketed_form(
                sequence, self.repeat_size, self.repeats
            )
        else:
            if re.search(r'\d{1,2}.1', self.canonical):
                if sequence[-1] == 'T':
                    forward_strand_brack_form = sequence_to_bracketed_form(
                        sequence, self.repeat_size, self.repeats
                    )
                else:
                    bf = sequence_to_bracketed_form(
                        sequence[1:], self.repeat_size, self.repeats
                    )
                    forward_strand_brack_form = f'{sequence[0]} {bf}'
            elif re.search(r'\d{1,2}.2', self.canonical):
                new_repeat_list = ['TATC', 'TGTC', 'AATC']
                forward_strand_brack_form = sequence_to_bracketed_form(
                    sequence, self.repeat_size, new_repeat_list
                )
            else:
                bf = sequence_to_bracketed_form(
                    sequence[3:], self.repeat_size, self.repeats
                )
                forward_strand_brack_form = f'{sequence[:3]} {bf}'
        return forward_strand_brack_form

    @property
    def designation(self):
        lus, sec, ter = None, None, None
        lus = repeat_copy_number(self.annotation, self.data['LUS'])
        sec = repeat_copy_number(self.annotation, self.data['Sec'])
        if str(self.annotation)[-1] == 'T' and isinstance(self.canonical, str):
            ter = 1
        else:
            ter = 0
        return lus, sec, ter


class STRMarker_D16S539(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return (
            f'{flank_seq[:2]} {flank_seq[2:6]} {flank_seq[6]} '
            f'{collapse_repeats_by_length(flank_seq[7:], 4)}'
        )

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        return (
            f'{collapse_repeats_by_length(flank_seq[:12], 4)} {flank_seq[12:15]} {flank_seq[15]} '
            f'{collapse_repeats_by_length(flank_seq[16:28], 4)} {flank_seq[28:31]} '
            f'{flank_seq[31:33]} {flank_seq[33]} {flank_seq[-2:]}'
        )


class STRMarker_D1S1656(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return f'{flank_seq[:3]} {collapse_repeats_by_length(flank_seq[3:], 4)}'

    @property
    def flank_3p(self):
        return ''

    @property
    def annotation(self):
        '''Bracketed annotation for D1S1656

        This function identifies if the sequence is a microvariant in order to call different
        functions to create the bracketed annotation.
        '''
        sequence = self.forward_sequence
        sequence_filt = sequence[2:]
        final = list()
        first_string, second_string = split_sequence_into_two_strings(sequence_filt, 'CACA')
        final.append(sequence[:2])
        if first_string == '':
            final.append('CACA')
        else:
            final.append(collapse_repeats_by_length(first_string, 4))
        if (len(second_string) % 4 != 0):
            final.append(sequence_to_bracketed_form(second_string, 4, self.repeats))
        else:
            final.append(collapse_repeats_by_length(second_string, 4))
        final_string = ' '.join(final)
        final_string = re.sub(r' +', ' ', final_string)
        return final_string


class STRMarker_PentaD(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return (
            f'{collapse_repeats_by_length(flank_seq[:20], 5)} {flank_seq[20]} {flank_seq[21:25]} '
            f'{collapse_repeats_by_length(flank_seq[25:], 5)}'
        )

    @property
    def designation(self):
        lus, sec, ter = super(STRMarker_PentaD, self).designation
        if self.canonical == '2.2':
            lus = 5
        elif self.canonical == '3.2':
            lus = 6
        return lus, sec, ter

    @property
    def annotation(self):
        '''Bracketed annotation for PentaD

        If the sequence is >= 18bp, the flanking region (first 5 bases) is first split off in the
        sequence to preserve that sequence. Then the repeat units are identified and bracketed.
        '''
        if len(self.uas_sequence) < 18:
            return sequence_to_bracketed_form(self.uas_sequence, self.repeat_size, self.repeats)
        else:
            prefix = self.uas_sequence[:5]
            suffix = self.uas_sequence[5:]
            brack_form = sequence_to_bracketed_form(suffix, self.repeat_size, self.repeats)
            result = f'{prefix} {brack_form}'
            result = re.sub(r' +', ' ', result)
            return result


class STRMarker_vWA(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return f'{flank_seq[:3]} {collapse_repeats_by_length(flank_seq[3:], 4)}'


class STRMarker_D10S1248(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return f'{flank_seq[:2]} {collapse_repeats_by_length(flank_seq[2:], 4)}'


class STRMarker_D22S1045(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return f'{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:], 3)}'


class STRMarker_FGA(STRMarker):
    @property
    def flank_3p(self):
        return ''

    @property
    def annotation(self):
        '''Bracketed annotation for FGA

        Specialized handling is required because which repeat unit should be identified differs
        based on its location in the sequence. For example, the 'GGAA' repeat should be identified
        at the beginning of the sequence; the 'GAAA' repeat should be identified at the end of the
        sequence; and the repeat 'AAAG' should be identified within the two end repeats.

        Simply identifying repeat units in a specified order does not result in the final
        annotation which is consistent with previously published annotation for this locus.
        '''
        sequence = self.forward_sequence
        if len(sequence) % self.repeat_size == 0 or (not ('GGAA') in sequence):
            return collapse_repeats_by_length(sequence, self.repeat_size)
        else:
            final = list()
            prev = 0
            if len(sequence) % 4 == 0:
                final_string = collapse_repeats_by_length(sequence, 4)
            else:
                for m in re.finditer('GGAA', sequence):
                    if prev == 0 or m.start() == prev:
                        prev = m.end()
                    else:
                        break
                first_string = sequence[:prev]
                second_string = sequence[prev:]
                prev = 0
                for m in re.finditer('AAAA', second_string):
                    prev = m.start()
                    break
                if second_string[prev:(prev+6)] == 'AAAAAA':
                    third_string = second_string[:prev+2]
                    fourth_string = second_string[prev+2:]
                elif prev == 0:
                    third_string = second_string[:-6]
                    fourth_string = second_string[-6:]
                else:
                    third_string = second_string[:prev]
                    fourth_string = second_string[prev:]
                final.append(collapse_repeats_by_length(first_string, 4))
                final.append(sequence_to_bracketed_form(third_string, 4, self.repeats))
                count = 0
                tmp = list()
                for element in re.split('GAAA', fourth_string):
                    parts = element.split(',')
                    for i in parts:
                        if i == '':
                            count += 1
                        else:
                            if count == 1:
                                tmp.append('GAAA')
                            elif count >= 2:
                                tmp.append('[GAAA]' + str(count))
                            count = 1
                            if i == 'AAAAAA':
                                tmp.append('AA AAAA')
                            elif len(i) > 4:
                                for x in split_by_n(i, 4):
                                    tmp.append(x)
                            else:
                                tmp.append(i)
                if parts[-1] == '' and count > 2:
                    tmp.append('[GAAA]' + str(count-1))
                elif parts[-1] == '' and count <= 2:
                    tmp.append('GAAA')
                last_string_final = ' '.join(tmp)
                final.append(last_string_final)
                final_string = ' '.join(final)
            return re.sub(r' +', ' ', final_string)


class STRMarker_CSF1PO(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        return f'{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:-1], 4)} {flank_seq[-1]}'


class STRMarker_D18S51(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        return (
            f'{flank_seq[:2]} {collapse_repeats_by_length(flank_seq[2:30], 4)} {flank_seq[30:33]} '
            f'{flank_seq[33]} {collapse_repeats_by_length(flank_seq[34:42], 4)} '
            f'{flank_seq[42:44]} {flank_seq[44:]}'
        )

    @property
    def annotation(self):
        if isinstance(self.canonical, str):
            return sequence_to_bracketed_form(self.uas_sequence, self.repeat_size, self.repeats)
        elif isinstance(self.canonical, int):
            return collapse_repeats_by_length(self.uas_sequence, self.repeat_size)


class STRMarker_D21S11(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        return (
            f'{flank_seq[:2]} {flank_seq[2]} {collapse_repeats_by_length(flank_seq[3:11], 4)} '
            f'{flank_seq[-1]}'
        )

    @property
    def annotation(self):
        '''Bracketed annotation for D21

        A specialized function is required for this locus due to the potential end of the sequence
        containing 'TA TCTA' and other variants. This sequence needs to remain intact to conform
        with the conventional annotation for this particular locus. However, if the 'TATCTA' is
        included in a repeat unit, the repeat unit needs to be reported (i.e. [TCTA]2).
        '''
        forward_strand_brack_form = sequence_to_bracketed_form(
            self.uas_sequence, self.data['NumBasesToSeparate'], self.repeats
        )
        prev = 0
        for m in re.finditer(']', forward_strand_brack_form):
            prev = m.end()
        if (
            prev == (len(forward_strand_brack_form) - 1) or
            prev == (len(forward_strand_brack_form) - 2) or
            prev == (len(forward_strand_brack_form) - 4) or
            prev == (len(forward_strand_brack_form) - 5)
           ):
            return forward_strand_brack_form
        else:
            first_string = forward_strand_brack_form[:prev+2]
            second_string = forward_strand_brack_form[prev+2:]
            second_string_final = re.sub(' ', '', second_string)
            if len(second_string_final) % 4 == 0:
                split_second_string = collapse_repeats_by_length(second_string_final, 4)
                final_string = f'{first_string} {second_string}'
            elif len(second_string_final) == 6:
                third_string = second_string_final[-6:-4]
                fourth_string = second_string_final[-4:]
                final_string = f'{first_string} {third_string} {fourth_string}'
            elif len(second_string_final) % 4 == 2:
                third_string = second_string_final[:-6]
                fourth_string = second_string_final[-6:-4]
                last_string = second_string_final[-4:]
                third_string_final = collapse_repeats_by_length(third_string, 4)
                final_string = f'{first_string} {third_string_final} {fourth_string} {last_string}'
            else:
                third_string = collapse_repeats_by_length(second_string_final, 4)
                final_string = f'{first_string} {third_string}'
            return re.sub('  ', ' ', final_string)

    @property
    def designation(self):
        '''Primary, secondary, and tertiary motif alleles for the D21S11 locus

        Special handling is required because the LUS repeat motif is the last 'TCTA' repeat set and
        the secondary repeat motif is the first set of 'TCTA' repeats in the sequence.
        '''
        sequence = self.annotation
        repeat = self.data['LUS']
        remaining = list()
        lus_sec = list()
        lus_allele = None
        for element in re.split(self.data['Tert'], sequence):
            if element == sequence:
                lus_sec = repeat_copy_number(element, repeat)
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
                    if i != '':
                        repeats = repeat_copy_number(i, repeat)
                        lus_sec.append(repeats)
        if lus_allele is None:
            lus_allele = lus_sec[1]
            sec_allele = lus_sec[0]

        finalcount = 0
        for m in re.finditer(self.data['Tert'], self.annotation):
            count = self.annotation[m.end()+1:m.end()+3]
            if count == '' or count[0] == '[' or count[0] == ' ' or count.isalpha():
                count = 1
            try:
                if float(count) > float(finalcount):
                    finalcount = count
                    try:
                        if str(finalcount)[1] == ' ':
                            finalcount = finalcount[0]
                    except IndexError:
                        count = count
            except ValueError:
                count = count
        return lus_allele, sec_allele, finalcount


class STRMarker_TH01(STRMarker):
    @property
    def annotation(self):
        '''Bracketed annotation for TH01

        Special handling is required for the microvariants of the TH01 locus because of the
        insertion of the 'ATG' between the repeat units 'AATG'.
        '''
        strings = collapse_all_repeats(self.uas_sequence, self.repeats)
        final_string = list()
        for unit in strings.split(' '):
            if '[' not in unit and len(unit) > 3 and (len(unit) % 4 != 0) and unit[:3] == 'ATG':
                group1 = unit[:3]
                final_string.append(group1)
                for x in split_by_n(unit[3:], n=4):
                    final_string.append(x)
            else:
                final_string.append(unit)
        final_form = ' '.join(final_string)
        return final_form


class STRMarker_D19S433(STRMarker):
    @property
    def annotation(self):
        '''Bracketed annotation for D19S433

        A specialized function is required for this locus. The sequence is first broken into two
        different strings. The two sets of sequences are processed separately in order to identify
        the potential presence of a deletion in either sequence.

        Simply identifying repeat units in a specified order does not result in the final
        annotation which is consistent with previously published annotation for this locus.
        '''
        sequence = self.forward_sequence
        final = list()
        last = 0
        prev = 0
        for m in re.finditer('CCTT', sequence):
            if m.start() == prev or m.start() == last:
                prev = m.end()
            else:
                last = m.end()
        final.append(sequence[:2])
        first_string = sequence[2:prev]
        second_string = sequence[prev:]
        if (len(first_string) % 4 != 0):
            final.append(sequence_to_bracketed_form(first_string, 4, self.repeats))
        else:
            final.append(collapse_repeats_by_length(first_string, 4))
        if (second_string != ""):
            if (len(second_string) % 4 != 0):
                if (len(second_string) > 6):
                    third_string = second_string[:-6]
                    final.append(collapse_repeats_by_length(third_string, 4))
                final.append(second_string[-6:-4])
                final.append(second_string[-4:])
            else:
                final.append(collapse_repeats_by_length(second_string, 4))
        final_string = ' '.join(final)
        return re.sub(r' +', ' ', final_string)


def STRMarkerObject(locus, sequence, uas=False, kit='forenseq'):
    constructors = {
        'D8S1179': STRMarker_D8S1179,
        'D13S317': STRMarker_D13S317,
        'D20S482': STRMarker_D20S482,
        'D2S441': STRMarker_D2S441,
        'D7S820': STRMarker_D7S820,
        'D16S539': STRMarker_D16S539,
        'D1S1656': STRMarker_D1S1656,
        'PentaD': STRMarker_PentaD,
        'vWA': STRMarker_vWA,
        'D10S1248': STRMarker_D10S1248,
        'D22S1045': STRMarker_D22S1045,
        'FGA': STRMarker_FGA,
        'CSF1PO': STRMarker_CSF1PO,
        'D18S51': STRMarker_D18S51,
        'D21S11': STRMarker_D21S11,
        'TH01': STRMarker_TH01,
        'D19S433': STRMarker_D19S433,
    }
    if locus in constructors:
        constructor = constructors[locus]
        return constructor(locus, sequence, uas=uas, kit=kit)
    else:
        return STRMarker(locus, sequence, uas=uas, kit=kit)
