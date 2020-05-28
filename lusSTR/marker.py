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
from lusSTR.repeat import collapse_repeats_by_length, sequence_to_bracketed_form
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
    def __init__(self, locus, sequence, uas=True, kit='forenseq'):
        self.locus = locus
        self.sequence = sequence
        if locus not in str_marker_data:
            raise InvalidLocusError(locus)
        self.data = str_marker_data[locus]
        self.uas = uas
        if kit.lower() not in ('forenseq', 'powerseq'):
            raise UnsupportedKitError(kit)
        self.kit = kit.lower()

    @property
    def repeat_size(self):
        return len(self.data['LUS'])

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
    def uas_sequence(self):
        '''Determine the UAS core sequence.'''
        if self.uas:
            return self.sequence
        front, back = self._uas_bases_to_trim()
        if back == 0:
            back = None
        else:
            back *= -1
        return self.sequence[front:back]

    @property
    def forward_sequence(self):
        if self.data['ReverseCompNeeded']:
            return lusSTR.annot.reverse_complement(self.uas_sequence)
        return self.uas_sequence

    @property
    def flankseq_5p(self):
        front, back = self._uas_bases_to_trim()
        return self.sequence[:front]

    @property
    def flank_5p(self):
        flank_rev = collapse_repeats_by_length(self.flankseq_5p, self.repeat_size)
        flank = flank_rev[::-1]
        return flank

    @property
    def flankseq_3p(self):
        front, back = self._uas_bases_to_trim()
        if back == 0:
            return ''
        return self.sequence[-back:]

    @property
    def flank_3p(self):
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
        lus, sec, ter = None, None, None
        if self.data['ReverseCompNeeded']:
            collapsed = collapse_repeats_by_length(self.forward_sequence, self.repeat_size)
        else:
            collapsed = sequence_to_bracketed_form(
                self.forward_sequence, self.repeat_size, self.data['Repeats']
            )
        lus = repeat_copy_number(collapsed, self.data['LUS'])
        if self.data['Sec'] != '':
            sec = repeat_copy_number(collapsed, self.data['Sec'])
        if self.data['Tert'] != '':
            ter = repeat_copy_number(collapsed, self.data['Tert'])
        return lus, sec, ter


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


class STRMarker_PentaD(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return (
            f'{collapse_repeats_by_length(flank_seq[:20], 5)} {flank_seq[20]} {flank_seq[21:25]} '
            f'{collapse_repeats_by_length(flank_seq[25:], 5)}'
        )


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


class STRMarker_D21S11(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        return (
            f'{flank_seq[:2]} {flank_seq[2]} {collapse_repeats_by_length(flank_seq[3:11], 4)} '
            f'{flank_seq[-1]}'
        )


def init_str_marker(locus, sequence, uas=True, kit='forenseq'):
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
    }
    if locus in constructors:
        constructor = constructors[locus]
        return constructor(locus, sequence, uas=uas, kit=kit)
    else:
        return STRMarker(locus, sequence, uas=uas, kit=kit)
