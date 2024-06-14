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

import importlib.resources
import json
import lusSTR
from lusSTR.scripts.repeat import collapse_repeats_by_length, collapse_repeats_by_length_flanks
from lusSTR.scripts.repeat import sequence_to_bracketed_form, split_sequence_into_two_strings
from lusSTR.scripts.repeat import reverse_complement, reverse_complement_bracketed
from lusSTR.scripts.repeat import repeat_copy_number, collapse_all_repeats, split_by_n
import re


def get_str_metadata_file():
    return importlib.resources.files("lusSTR") / "data/str_markers.json"


with open(get_str_metadata_file(), "r") as fh:
    str_marker_data = json.load(fh)


class InvalidLocusError(ValueError):
    pass


class UnsupportedKitError(ValueError):
    pass


class InvalidSequenceError(ValueError):
    pass


class UnsupportedSoftwareError(ValueError):
    pass


class STRMarker:
    def __init__(self, locus, sequence, software, custom=False, kit="forenseq"):
        self.locus = locus
        self.sequence = sequence
        if locus not in str_marker_data:
            raise InvalidLocusError(locus)
        self.data = str_marker_data[locus]
        if software.lower() not in ("uas", "straitrazor", "genemarker"):
            raise UnsupportedSoftwareError(software)
        self.software = software
        if kit.lower() not in ("forenseq", "powerseq"):
            raise UnsupportedKitError(kit)
        self.kit = kit.lower()
        if software == "uas" and self.data["ReverseCompNeeded"] == "Yes":
            self.sequence = reverse_complement(sequence)
        self.custom = custom

    @property
    def repeat_size(self):
        return len(self.data["LUS"])

    @property
    def repeats(self):
        return self.data["Repeats"]

    def _uas_bases_to_trim(self):
        """Number of bases to trim off each side to get the UAS sequence.

        ForenSeq and PowerSeq amplicons extend beyond the core UAS sequence for some loci. This
        function determines the number of bases that need to be trimmed from the full amplicon
        sequence to recover the UAS core sequence.
        """
        if self.software == "uas":
            return 0, 0
        elif self.kit == "forenseq":
            return self.data["Foren_5"], self.data["Foren_3"]
        elif self.kit == "powerseq":
            if self.locus == "D16S539" and self.software == "genemarker":
                return self.data["Power_5"], (self.data["Power_3"] - 3)
            elif self.locus == "D8S1179" and self.software == "genemarker":
                return (self.data["Power_5"] - 5), (self.data["Power_3"] - 5)
            else:
                return self.data["Power_5"], self.data["Power_3"]
        else:
            raise UnsupportedKitError(self.kit)

    @property
    def forward_sequence(self):
        """Sequence from the UAS region, in the forward orientation

        If it is the full ForenSeq or PowerSeq amplicon sequence, trim the full forward sequence
        back to the UAS region. If the sequence has already been run through UAS, no trimming is
        required.
        """
        if self.software == "uas":
            return self.sequence
        front, back = self._uas_bases_to_trim()
        if back == 0:
            back = None
        else:
            back *= -1
        return self.sequence[front:back]

    @property
    def uas_sequence(self):
        """Sequence from the UAS region, matching the orientation of the UAS output

        The UAS software outputs the reverse complement of the forward sequence for some loci.
        """
        if self.data["ReverseCompNeeded"] == "Yes":
            return reverse_complement(self.forward_sequence)
        return self.forward_sequence

    @property
    def custom_sequence(self):
        """Custom range for sequences; PowerSeq sequences only"""
        if self.custom:
            front, back = self._uas_bases_to_trim()
            custom_front = front - self.data["Custom_5"]
            custom_back = back - self.data["Custom_3"]
            if custom_back == 0:
                custom_back = None
            else:
                custom_back *= -1
            return self.sequence[custom_front:custom_back]
        else:
            return None

    @property
    def flankseq_5p(self):
        if self.software == "uas":
            return None
        front, back = self._uas_bases_to_trim()
        if front == 0:
            return ""
        return self.sequence[:front]

    @property
    def flank_5p(self):
        if self.software == "uas" or self.flankseq_5p == "":
            return None
        elif (
            self.kit == "powerseq"
            and self.data["Power_5"] > self.data["Foren_5"]
            and self.data["Foren_5"] > 0
        ):
            power_seq_flank = collapse_repeats_by_length_flanks(
                self.flankseq_5p[: -self.data["Foren_5"]], self.repeat_size
            )
            foren_seq_flank = collapse_repeats_by_length_flanks(
                self.flankseq_5p[-self.data["Foren_5"] :], self.repeat_size
            )
            flank = f"{power_seq_flank} {foren_seq_flank}"
        else:
            flank = collapse_repeats_by_length_flanks(self.flankseq_5p, self.repeat_size)
        return flank

    @property
    def flankseq_3p(self):
        if self.software == "uas":
            return None
        front, back = self._uas_bases_to_trim()
        if back == 0:
            return ""
        return self.sequence[-back:]

    @property
    def flank_3p(self):
        if self.software == "uas" or self.flankseq_3p == "":
            return None
        elif (
            self.kit == "powerseq"
            and self.data["Power_3"] > self.data["Foren_3"]
            and self.data["Foren_3"] > 0
        ):
            foren_seq_flank = collapse_repeats_by_length(
                self.flankseq_3p[: self.data["Foren_3"]], self.repeat_size
            )
            power_seq_flank = collapse_repeats_by_length(
                self.flankseq_3p[self.data["Foren_3"] :], self.repeat_size
            )
            flank = f"{foren_seq_flank} {power_seq_flank}"
        else:
            flank = collapse_repeats_by_length(self.flankseq_3p, self.repeat_size)
        return flank

    @property
    def canonical(self):
        """Canonical STR allele designation"""
        n = self.repeat_size
        nsubout = self.data["BasesToSubtract"]
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
            canon_allele = f"{allele_int}.{allele_dec}"
        return canon_allele

    @property
    def indel_flag(self):
        powerseq_loci = ["DYS393", "DYS458", "DYS456"]
        partial_loci = ["DYS439", "DYS391", "DYS385A-B"]
        """Check for potential indels within flanking regions"""
        if str(self.canonical) not in self.data["Alleles"]:
            if self.locus in powerseq_loci:
                flag = "UAS region indicates entire sequence; Possible indel or partial sequence"
            else:
                flag = "Possible indel or partial sequence"
        else:
            if self.locus in powerseq_loci:
                flag = "UAS region indicates entire sequence"
            elif (
                self.locus in partial_loci and self.kit == "powerseq"
            ) or self.locus == "Y-GATA-H4":
                flag = "Partial UAS region sequence"
            else:
                flag = " "
        return flag

    @property
    def cannot_split(self):
        return self.locus in [
            "D19S433",
            "D6S1043",
            "TH01",
            "D21S11",
            "D1S1656",
            "D7S820",
            "D5S818",
            "D12S391",
            "D9S1122",
            "PENTA E",
            "DXS7132",
        ]

    @property
    def must_split(self):
        return self.locus in [
            "D13S317",
            "D18S51",
            "DYS643",
            "DYS635",
            "DYS635",
            "DYS612",
            "DYS576",
            "DYS570",
            "DYS549",
            "DYS533",
            "DYS505",
            "DYS481",
            "DYS460",
            "DYS439",
            "DYS438",
            "DYS437",
            "DYS392",
            "DYS391",
            "DYS390",
            "DYS389II",
            "DYS389I",
            "DYS385A-B",
            "DYS19",
            "DYF387S1",
            "DYS393",
            "DYS456",
            "HPRTB",
            "DXS8378",
            "DXS7423",
            "DXS10103",
        ]

    @property
    def split_compatible(self):
        return self.must_split or len(self.uas_sequence) % self.repeat_size == 0

    @property
    def do_split(self):
        return not self.cannot_split or self.split_compatible

    @property
    def convert(self):
        bylength = (
            self.split_compatible
            or (self.data["ReverseCompNeeded"] == "Yes" and self.split_compatible)
            or (self.locus == "D3S1358" and self.split_compatible)
            or self.locus == "D16S539"
        )
        if bylength:
            collapseseq = collapse_repeats_by_length(self.forward_sequence, self.repeat_size)
        else:
            collapseseq = sequence_to_bracketed_form(
                self.forward_sequence, self.repeat_size, self.repeats
            )
        return collapseseq

    @property
    def convert_uas(self):
        if self.data["ReverseCompNeeded"] == "Yes":
            return reverse_complement_bracketed(self.convert)
        return self.convert

    @property
    def custom_brack(self):
        return self.convert

    @property
    def designation(self):
        lus, sec, ter = None, None, None
        lus = repeat_copy_number(self.convert, self.data["LUS"])
        if self.data["Sec"] != "":
            sec = repeat_copy_number(self.convert, self.data["Sec"])
        if self.data["Tert"] != "":
            ter = repeat_copy_number(self.convert, self.data["Tert"])
        return lus, sec, ter

    @property
    def summary(self):
        lus, sec, ter = self.designation
        canon = self.canonical
        lus_final_output = f"{canon}_{lus}"
        if sec is None:
            lus_plus = lus_final_output
        else:
            if ter is None:
                lus_plus = f"{canon}_{lus}_{sec}"
            else:
                lus_plus = f"{canon}_{lus}_{sec}_{ter}"
        return [
            self.uas_sequence,
            self.forward_sequence,
            self.custom_sequence,
            self.convert_uas,
            self.convert,
            self.custom_brack,
            canon,
            lus_final_output,
            lus_plus,
        ]


class STRMarker_D8S1179(STRMarker):
    @property
    def flank_5p(self):
        if self.kit == "powerseq":
            flank = collapse_repeats_by_length_flanks(self.flankseq_5p, 4)
        else:
            flank = ""
        return flank

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[: self.data["Foren_3"]]
            foren_flank_anno = collapse_repeats_by_length(foren_seq_flank, 4)
            power_seq_flank = flank_seq[self.data["Foren_3"] :]
            power_flank_anno = (
                f"{collapse_repeats_by_length(power_seq_flank[:41], 4)} "
                f"{collapse_repeats_by_length(power_seq_flank[41:62], 4)} "
                f"{collapse_repeats_by_length(power_seq_flank[62:], 4)}"
            )
            flank = f"{foren_flank_anno} {power_flank_anno}"
        else:
            flank = collapse_repeats_by_length(flank_seq, 4)
        return flank


class STRMarker_D13S317(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[-self.data["Foren_5"] :]
            foren_flank_anno = (
                f"{foren_seq_flank[:2]} {collapse_repeats_by_length(foren_seq_flank[2:14], 4)} "
                f"{foren_seq_flank[14]} {foren_seq_flank[15]} {foren_seq_flank[16:19]} "
                f"{collapse_repeats_by_length(foren_seq_flank[19:], 4)}"
            )
            power_seq_flank = flank_seq[: -self.data["Foren_5"]]
            power_flank_anno = (
                f"{power_seq_flank[:2]} "
                f"{collapse_repeats_by_length_flanks(power_seq_flank[2:17], 4)} "
                f"{power_seq_flank[17]} {power_seq_flank[18:22]} {power_seq_flank[22:27]} "
                f"{power_seq_flank[27:]}"
            )
            flank = f"{power_flank_anno} {foren_flank_anno}"
        else:
            flank = (
                f"{flank_seq[:2]} {collapse_repeats_by_length(flank_seq[2:14], 4)} "
                f"{flank_seq[14]} {flank_seq[15]} {flank_seq[16:19]} "
                f"{collapse_repeats_by_length(flank_seq[19:], 4)}"
            )
        return flank

    @property
    def convert(self):
        if len(self.uas_sequence) < 110:
            bracketed_form = collapse_repeats_by_length(self.uas_sequence, 4)
        else:
            if "GGGCTGCCTA" in self.uas_sequence:
                break_point = self.uas_sequence.index("GGGCTGCCTA") + 10
                bracketed_form = (
                    f"{collapse_repeats_by_length(self.uas_sequence[:break_point], 4)} "
                    f"{collapse_repeats_by_length(self.uas_sequence[break_point:], 4)}"
                )
            elif "TTTT" in self.uas_sequence:
                break_point = self.uas_sequence.index("TTTT") + 14
                bracketed_form = (
                    f"{collapse_repeats_by_length(self.uas_sequence[:break_point], 4)} "
                    f"{collapse_repeats_by_length(self.uas_sequence[break_point:], 4)}"
                )
            else:
                bracketed_form = ""
        return bracketed_form


class STRMarker_D20S482(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        return (
            f"{flank_seq[:2]} {flank_seq[2:6]} {flank_seq[6]} {flank_seq[7:10]} "
            f"{flank_seq[10:14]} {flank_seq[14:18]}"
        )


class STRMarker_D2S441(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[-self.data["Foren_5"] :]
            foren_flank_anno = (
                f"{foren_seq_flank[:4]} {foren_seq_flank[4]} "
                f"{collapse_repeats_by_length(foren_seq_flank[5:], 4)}"
            )
            power_seq_flank = flank_seq[: -self.data["Foren_5"]]
            power_flank_anno = (
                f"{power_seq_flank[:2]} {collapse_repeats_by_length(power_seq_flank[2:45], 4)}"
                f" {power_seq_flank[45]} {collapse_repeats_by_length(power_seq_flank[-6:], 4)}"
            )
            flank = f"{power_flank_anno} {foren_flank_anno}"
        else:
            flank = (
                f"{flank_seq[:4]} {flank_seq[4]} {collapse_repeats_by_length(flank_seq[5:], 4)}"
            )
        return flank


class STRMarker_D5S818(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:7], 4)} "
                f"{collapse_repeats_by_length(flank_seq[7:13], 4)} "
                f"{collapse_repeats_by_length(flank_seq[13:58], 4)} "
                f"{collapse_repeats_by_length(flank_seq[58:], 4)}"
            )
        else:
            flank = collapse_repeats_by_length(flank_seq, 4)
        return flank


class STRMarker_D7S820(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[-self.data["Foren_5"] :]
            foren_flank_anno = (
                f"{foren_seq_flank[0]} {collapse_repeats_by_length(foren_seq_flank[1:13], 4)} "
                f"{foren_seq_flank[13:]}"
            )
            power_seq_flank = flank_seq[: -self.data["Foren_5"]]
            power_flank_anno = (
                f"{collapse_repeats_by_length(power_seq_flank[:13], 4)} "
                f"{collapse_repeats_by_length(power_seq_flank[13:38], 4)} {power_seq_flank[-4:]}"
            )
            flank = f"{power_flank_anno} {foren_flank_anno}"
        else:
            flank = (
                f"{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:13], 4)} "
                f"{flank_seq[13:]}"
            )
        return flank

    @property
    def convert(self):
        """
        Function to correctly bracket microvariants in the D7S820 locus.

        Microvariants N.1 and N.3 separate the first base (or first three bases) of the sequence
        in order to correctly bracket the sequence. Microvariant N.2 has a third repeat motif
        added in order to bracket the sequence correctly.
        """
        sequence = self.forward_sequence
        if type(self.canonical) == int:
            forward_strand_brack_form = sequence_to_bracketed_form(
                sequence, self.repeat_size, self.repeats
            )
        else:
            if re.search(r"\d{1,2}.1", self.canonical):
                if sequence[-1] == "T":
                    forward_strand_brack_form = sequence_to_bracketed_form(
                        sequence, self.repeat_size, self.repeats
                    )
                else:
                    bf = sequence_to_bracketed_form(sequence[1:], self.repeat_size, self.repeats)
                    forward_strand_brack_form = f"{sequence[0]} {bf}"
            elif re.search(r"\d{1,2}.2", self.canonical):
                new_repeat_list = ["TATC", "TGTC", "AATC"]
                forward_strand_brack_form = sequence_to_bracketed_form(
                    sequence, self.repeat_size, new_repeat_list
                )
            else:
                bf = sequence_to_bracketed_form(sequence[3:], self.repeat_size, self.repeats)
                forward_strand_brack_form = f"{sequence[:3]} {bf}"
        return forward_strand_brack_form

    @property
    def designation(self):
        lus, sec, ter = None, None, None
        lus = repeat_copy_number(self.convert, self.data["LUS"])
        sec = repeat_copy_number(self.convert, self.data["Sec"])
        if str(self.convert)[-1] == "T" and isinstance(self.canonical, str):
            ter = 1
        else:
            ter = 0
        return lus, sec, ter


class STRMarker_D16S539(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[-self.data["Foren_5"] :]
            foren_flank_anno = (
                f"{foren_seq_flank[:2]} {foren_seq_flank[2:6]} {foren_seq_flank[6]} "
                f"{collapse_repeats_by_length(foren_seq_flank[7:], 4)}"
            )
            power_seq_flank = flank_seq[: -self.data["Foren_5"]]
            power_flank_anno = (
                f"{power_seq_flank[:2]} {collapse_repeats_by_length(power_seq_flank[2:44], 4)} "
                f"{collapse_repeats_by_length(power_seq_flank[44:74], 4)} "
                f"{collapse_repeats_by_length(power_seq_flank[74:95], 4)} {power_seq_flank[-4:]}"
            )
            flank = f"{power_flank_anno} {foren_flank_anno}"
        else:
            flank = (
                f"{flank_seq[:2]} {flank_seq[2:6]} {flank_seq[6]} "
                f"{collapse_repeats_by_length(flank_seq[7:], 4)}"
            )
        return flank

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = flank_seq
        else:
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:12], 4)} {flank_seq[12:15]} "
                f"{flank_seq[15]} {collapse_repeats_by_length(flank_seq[16:28], 4)} "
                f"{flank_seq[28:31]} {flank_seq[31:33]} {flank_seq[33]} {flank_seq[-2:]}"
            )
        return flank


class STRMarker_D1S1656(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = f"{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:], 4)}"
        else:
            flank = f"{flank_seq[:3]} {collapse_repeats_by_length(flank_seq[3:], 4)}"
        return flank

    @property
    def flank_3p(self):
        if self.kit == "powerseq":
            flank_seq = self.flankseq_3p
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:6], 4)} {flank_seq[6:9]} "
                f"{collapse_repeats_by_length(flank_seq[9:], 4)}"
            )
        else:
            flank = ""
        return flank

    @property
    def convert(self):
        """Bracketed sequence form for D1S1656

        This function identifies if the sequence is a microvariant in order to call different
        functions to create the bracketed form.
        """
        sequence = self.forward_sequence
        sequence_filt = sequence[2:]
        final = list()
        first_string, second_string = split_sequence_into_two_strings(sequence_filt, "CACA")
        final.append(sequence[:2])
        if first_string == "":
            final.append("CACA")
        else:
            final.append(collapse_repeats_by_length(first_string, 4))
        if len(second_string) % 4 != 0:
            final.append(sequence_to_bracketed_form(second_string, 4, self.repeats))
        else:
            final.append(collapse_repeats_by_length(second_string, 4))
        final_string = " ".join(final)
        final_string = re.sub(r" +", " ", final_string)
        return final_string


class STRMarker_PentaD(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = collapse_repeats_by_length_flanks(flank_seq, 5)
        else:
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:20], 5)} {flank_seq[20]} "
                f"{flank_seq[21:25]} {collapse_repeats_by_length(flank_seq[25:], 5)}"
            )
        return flank

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[: self.data["Foren_3"]]
            foren_flank_anno = collapse_repeats_by_length(foren_seq_flank, 5)
            power_seq_flank = flank_seq[self.data["Foren_3"] :]
            power_flank_anno = f"{collapse_repeats_by_length(power_seq_flank, 5)}"
            flank = f"{foren_flank_anno} {power_flank_anno}"
        else:
            flank = collapse_repeats_by_length(flank_seq, 5)
        return flank

    @property
    def designation(self):
        lus, sec, ter = super(STRMarker_PentaD, self).designation
        if self.canonical == "2.2":
            lus = 5
        elif self.canonical == "3.2":
            lus = 6
        return lus, sec, ter

    @property
    def convert(self):
        """Bracketed sequence form for PentaD

        If the sequence is >= 18bp, the flanking region (first 5 bases) is first split off in the
        sequence to preserve that sequence. Then the repeat units are identified and bracketed.
        """
        if len(self.uas_sequence) < 18:
            return sequence_to_bracketed_form(self.uas_sequence, self.repeat_size, self.repeats)
        else:
            prefix = self.uas_sequence[:5]
            suffix = self.uas_sequence[5:]
            brack_form = sequence_to_bracketed_form(suffix, self.repeat_size, self.repeats)
            result = f"{prefix} {brack_form}"
            result = re.sub(r" +", " ", result)
            return result


class STRMarker_PentaE(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[-self.data["Foren_5"] :]
            foren_flank_anno = collapse_repeats_by_length_flanks(foren_seq_flank, 5)
            power_seq_flank = flank_seq[: -self.data["Foren_5"]]
            power_flank_anno = collapse_repeats_by_length_flanks(power_seq_flank, 5)
            flank = f"{power_flank_anno} {foren_flank_anno}"
        else:
            flank = collapse_repeats_by_length_flanks(flank_seq, 5)
        return flank

    @property
    def flank_3p(self):
        if self.kit == "powerseq":
            flank = self.flankseq_3p
        else:
            flank = collapse_repeats_by_length(self.flankseq_3p, 5)
        return flank


class STRMarker_vWA(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = collapse_repeats_by_length_flanks(flank_seq, 4)
        else:
            flank = f"{flank_seq[:3]} {collapse_repeats_by_length(flank_seq[3:], 4)}"
        return flank

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:5], 4)} "
                f"{collapse_repeats_by_length(flank_seq[5:], 4)}"
            )
        else:
            flank = collapse_repeats_by_length(flank_seq, 4)
        return flank


class STRMarker_D10S1248(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            power_seq_flank = collapse_repeats_by_length_flanks(
                flank_seq[: -self.data["Foren_5"]], 4
            )
            foren_seq_flank = collapse_repeats_by_length_flanks(
                flank_seq[-self.data["Foren_5"] :], 4
            )
            flank = f"{power_seq_flank} {foren_seq_flank}"
        else:
            flank = collapse_repeats_by_length_flanks(flank_seq, 4)
        return flank

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "forenseq":
            flank = collapse_repeats_by_length(flank_seq, 4)
        else:
            flank = ""
        return flank


class STRMarker_D22S1045(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = collapse_repeats_by_length_flanks(flank_seq, 3)
        else:
            flank = f"{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:], 3)}"
        return flank


class STRMarker_D2S1338(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length_flanks(flank_seq[:30], 4)} {flank_seq[30]} "
                f"{collapse_repeats_by_length(flank_seq[31:-3], 4)} {flank_seq[-3:]}"
            )
        else:
            flank = flank_seq
        return flank


class STRMarker_FGA(STRMarker):
    @property
    def flank_3p(self):
        return ""

    @property
    def convert(self):
        """Bracketed sequence form for FGA

        Specialized handling is required because which repeat unit should be identified differs
        based on its location in the sequence. For example, the 'GGAA' repeat should be identified
        at the beginning of the sequence; the 'GAAA' repeat should be identified at the end of the
        sequence; and the repeat 'AAAG' should be identified within the two end repeats.

        Simply identifying repeat units in a specified order does not result in the final
        form which is consistent with the previously published sequence form for this locus.
        """
        sequence = self.forward_sequence
        if len(sequence) % self.repeat_size == 0 or (not ("GGAA") in sequence):
            return collapse_repeats_by_length(sequence, self.repeat_size)
        else:
            final = list()
            prev = 0
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
            if second_string[prev : (prev + 6)] == "AAAAAA":
                third_string = second_string[: prev + 2]
                fourth_string = second_string[prev + 2 :]
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
            for element in re.split("GAAA", fourth_string):
                parts = element.split(",")
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
                            for x in split_by_n(i, 4, False):
                                tmp.append(x)
                        else:
                            tmp.append(i)
            if parts[-1] == "" and count > 2:
                tmp.append("[GAAA]" + str(count - 1))
            elif parts[-1] == "" and count <= 2:
                tmp.append("GAAA")
            last_string_final = " ".join(tmp)
            final.append(last_string_final)
            final_string = " ".join(final)
            return re.sub(r" +", " ", final_string)


class STRMarker_CSF1PO(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[: self.data["Foren_3"]]
            foren_flank_anno = (
                f"{foren_seq_flank[0]} {collapse_repeats_by_length(foren_seq_flank[1:-1], 4)}"
                f" {foren_seq_flank[-1]}"
            )
            power_seq_flank = flank_seq[self.data["Foren_3"] :]
            power_flank_anno = (
                f"{collapse_repeats_by_length(power_seq_flank[:71], 4)} {power_seq_flank[71:73]}"
                f" {collapse_repeats_by_length(power_seq_flank[73:], 4)}"
            )
            flank = f"{foren_flank_anno} {power_flank_anno}"
        else:
            flank = (
                f"{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:-1], 4)} {flank_seq[-1]}"
            )
        return flank


class STRMarker_D18S51(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{flank_seq[:2]} {collapse_repeats_by_length(flank_seq[2:33], 4)} "
                f"{flank_seq[33]} {collapse_repeats_by_length(flank_seq[33:], 4)}"
            )
        else:
            flank = (
                f"{flank_seq[:2]} {collapse_repeats_by_length(flank_seq[2:30], 4)} "
                f"{flank_seq[30:33]} {flank_seq[33]} "
                f"{collapse_repeats_by_length(flank_seq[34:42], 4)} {flank_seq[42:44]} "
                f"{flank_seq[44:]}"
            )
        return flank

    @property
    def convert(self):
        if isinstance(self.canonical, str):
            return sequence_to_bracketed_form(self.uas_sequence, self.repeat_size, self.repeats)
        elif isinstance(self.canonical, int):
            return collapse_repeats_by_length(self.uas_sequence, self.repeat_size)


class STRMarker_D21S11(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            foren_seq_flank = flank_seq[: self.data["Foren_3"]]
            foren_flank_anno = (
                f"{foren_seq_flank[:2]} {foren_seq_flank[2]} "
                f"{collapse_repeats_by_length(foren_seq_flank[3:11], 4)} {flank_seq[-1]}"
            )
            power_seq_flank = flank_seq[self.data["Foren_3"] :]
            power_flank_anno = (
                f"{collapse_repeats_by_length(power_seq_flank[:21], 4)} "
                f"{collapse_repeats_by_length(power_seq_flank[21:], 4)}"
            )
            flank = f"{foren_flank_anno} {power_flank_anno}"
        else:
            flank = (
                f"{flank_seq[:2]} {flank_seq[2]} "
                f"{collapse_repeats_by_length(flank_seq[3:11], 4)} {flank_seq[-1]}"
            )
        return flank

    @property
    def convert(self):
        """Bracketed sequence form for D21

        A specialized function is required for this locus due to the potential end of the sequence
        containing 'TA TCTA' and other variants. This sequence needs to remain intact to conform
        with the conventional bracketed form for this particular locus. However, if the 'TATCTA'
        is included in a repeat unit, the repeat unit needs to be reported (i.e. [TCTA]2).
        """
        forward_strand_brack_form = sequence_to_bracketed_form(
            self.uas_sequence, self.data["NumBasesToSeparate"], self.repeats
        )
        prev = 0
        for m in re.finditer("]", forward_strand_brack_form):
            prev = m.end()
        if (
            prev == (len(forward_strand_brack_form) - 1)
            or prev == (len(forward_strand_brack_form) - 2)
            or prev == (len(forward_strand_brack_form) - 4)
            or prev == (len(forward_strand_brack_form) - 5)
        ):
            return forward_strand_brack_form
        else:
            first_string = forward_strand_brack_form[: prev + 2]
            second_string = forward_strand_brack_form[prev + 2 :]
            second_string_final = re.sub(" ", "", second_string)
            if len(second_string_final) % 4 == 0:
                split_second_string = collapse_repeats_by_length(second_string_final, 4)
                final_string = f"{first_string} {split_second_string}"
            elif len(second_string_final) == 6:
                third_string = second_string_final[-6:-4]
                fourth_string = second_string_final[-4:]
                final_string = f"{first_string} {third_string} {fourth_string}"
            elif len(second_string_final) % 4 == 2:
                third_string = second_string_final[:-6]
                fourth_string = second_string_final[-6:-4]
                last_string = second_string_final[-4:]
                third_string_final = collapse_repeats_by_length(third_string, 4)
                final_string = f"{first_string} {third_string_final} {fourth_string} {last_string}"
            else:
                third_string = collapse_repeats_by_length(second_string_final, 4)
                final_string = f"{first_string} {third_string}"
            return re.sub("  ", " ", final_string)

    @property
    def designation(self):
        """Primary, secondary, and tertiary motif alleles for the D21S11 locus

        Special handling is required because the LUS repeat motif is the last 'TCTA' repeat set and
        the secondary repeat motif is the first set of 'TCTA' repeats in the sequence.
        """
        sequence = self.convert
        repeat = self.data["LUS"]
        remaining = list()
        lus_sec = list()
        lus_allele = None
        for element in re.split(self.data["Tert"], sequence):
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
                parts = element.split("[,]")
                for i in parts:
                    if i != "":
                        repeats = repeat_copy_number(i, repeat)
                        lus_sec.append(repeats)
        if lus_allele is None:
            if len(lus_sec) == 2:
                lus_allele = lus_sec[1]
                sec_allele = lus_sec[0]
            elif len(lus_sec) > 2:
                lus_allele = lus_sec[-1]
                sec_allele = lus_sec[0]
            else:
                lus_allele = 0
                sec_allele = lus_sec[0]
        finalcount = 0
        for m in re.finditer(self.data["Tert"], self.convert):
            count = self.convert[m.end() + 1 : m.end() + 3]
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
        return lus_allele, sec_allele, finalcount


class STRMarker_TH01(STRMarker):
    @property
    def convert(self):
        """Bracketed sequence form for TH01

        Special handling is required for the microvariants of the TH01 locus because of the
        insertion of the 'ATG' between the repeat units 'AATG'.
        """
        strings = collapse_all_repeats(self.uas_sequence, self.repeats)
        final_string = list()
        for unit in strings.split(" "):
            if "[" not in unit and len(unit) > 3 and len(unit) % 4 != 0 and unit[:3] == "ATG":
                group1 = unit[:3]
                final_string.append(group1)
                for x in split_by_n(unit[3:], n=4, rev=False):
                    final_string.append(x)
            elif "[" not in unit and len(unit) > 4:
                final_string.append(collapse_repeats_by_length(unit, 4))
            else:
                final_string.append(unit)
        final_form = " ".join(final_string)
        return final_form

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:129], 4)} "
                f"{collapse_repeats_by_length(flank_seq[129:], 4)}"
            )
        else:
            flank = collapse_repeats_by_length(flank_seq, 4)
        return flank


class STRMarker_TPOX(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = (
                f"{flank_seq[:3]} {collapse_repeats_by_length(flank_seq[3:16], 4)} "
                f"{collapse_repeats_by_length(flank_seq[16:23], 4)} "
                f"{collapse_repeats_by_length(flank_seq[23:29], 4)} "
                f"{collapse_repeats_by_length(flank_seq[29:72], 4)} "
                f"{flank_seq[72:74]} {collapse_repeats_by_length(flank_seq[74:-2], 4)} "
                f"{flank_seq[-2:]}"
            )
        else:
            flank = flank_seq
        return flank

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:5], 4)} "
                f"{collapse_repeats_by_length(flank_seq[5:], 4)}"
            )
        else:
            flank = collapse_repeats_by_length(flank_seq, 4)
        return flank


class STRMarker_D19S433(STRMarker):
    @property
    def convert(self):
        """Bracketed sequence form for D19S433

        A specialized function is required for this locus. The sequence is first broken into two
        different strings. The two sets of sequences are processed separately in order to identify
        the potential presence of a deletion in either sequence.

        Simply identifying repeat units in a specified order does not result in the final
        bracketed form which is consistent with the previously published bracketed form for this
        locus.
        """
        sequence = self.forward_sequence
        final = list()
        last = 0
        prev = 0
        for m in re.finditer("CCTT", sequence):
            if m.start() == prev or m.start() == last:
                prev = m.end()
            else:
                last = m.end()
        final.append(sequence[:2])
        if prev != 0:
            first_string = sequence[2:prev]
            second_string = sequence[prev:]
            if len(first_string) % 4 != 0:
                final.append(sequence_to_bracketed_form(first_string, 4, self.repeats))
            else:
                final.append(collapse_repeats_by_length(first_string, 4))
        else:
            second_string = sequence[2:]
        if second_string != "":
            if len(second_string) % 4 != 0:
                if len(second_string) > 6:
                    third_string = second_string[:-6]
                    final.append(collapse_repeats_by_length(third_string, 4))
                final.append(second_string[-6:-4])
                final.append(second_string[-4:])
            else:
                final.append(collapse_repeats_by_length(second_string, 4))
        final_string = " ".join(final)
        return re.sub(r" +", " ", final_string)


class STRMarker_DYS643(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:8], 5)} "
                f"{collapse_repeats_by_length(flank_seq[8:], 5)}"
            )
        else:
            flank = collapse_repeats_by_length(flank_seq, 5)
        return flank


class STRMarker_DYS635(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length_flanks(flank_seq[:11], 4)} {flank_seq[11]} "
                f"{collapse_repeats_by_length(flank_seq[12:42], 4)}"
            )
        else:
            flank = f"{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:31], 4)}"
        return flank


class STRMarker_DYS612(STRMarker):
    @property
    def designation(self):
        """
        The LUS and Secondary motif allele are the same, "TCT".
        The LUS motif is always the last "TCT" in the sequence.
        The secondary motif, however, is not as easily identified given there could be
        multiple instances of "TCT" within the sequence. After the LUS motif is identified,
        the LUS allele is identified as the "TCT" repeat sequence with the largest number of
        repeats.
        """
        lus, sec, ter = None, None, None
        anno = self.convert
        repeat = "TCT"
        match_list = []
        for block in anno.split(" "):
            if block == repeat:
                match_list.append(1)
            match = re.match(r"\[" + repeat + r"\](\d+)", block)
            if match:
                length = int(match.group(1))
                match_list.append(length)
        if len(match_list) == 1:
            lus = match_list[0]
            sec = 0
        elif len(match_list) > 1:
            lus = match_list[-1]
            sec = match_list[-2]
        else:
            pass
        return lus, sec, ter

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        flank = f"{flank_seq[0]} {collapse_repeats_by_length(flank_seq[1:], 3)}"
        return flank


class STRMarker_DYS576(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length(flank_seq[-62:-23], 4)} "
                f"{collapse_repeats_by_length(flank_seq[-23:], 4)}"
            )
        else:
            flank = collapse_repeats_by_length(flank_seq, 4)
        return flank


class STRMarker_DYS549(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length_flanks(flank_seq[:49], 4)} {flank_seq[49:51]} "
                f"{collapse_repeats_by_length(flank_seq[51:], 4)}"
            )
        else:
            flank = f"{flank_seq[:2]} {collapse_repeats_by_length(flank_seq[2:], 4)}"
        return flank


class STRMarker_DYS533(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:21], 4)} "
                f"{collapse_repeats_by_length(flank_seq[21:], 4)}"
            )
        else:
            flank = collapse_repeats_by_length(flank_seq, 4)
        return flank


class STRMarker_DYS522(STRMarker):
    @property
    def convert(self):
        sequence = self.forward_sequence
        final_seq = f"{sequence[:3]} {collapse_repeats_by_length(sequence[3:], 4)}"
        return final_seq


class STRMarker_DYS439(STRMarker):
    @property
    def canonical(self):
        """Canonical STR allele designation"""
        n = self.repeat_size
        if self.kit == "powerseq":
            nsubout = self.data["BasesToSubtract"] - 48
        else:
            nsubout = self.data["BasesToSubtract"]
        nsubout *= -1
        new_seq = self.uas_sequence[:nsubout]
        if len(new_seq) % n == 0:
            canon_allele = int(len(new_seq) / n)
        else:
            allele_int = int(len(new_seq) / n)
            allele_dec = int(len(new_seq) % n)
            canon_allele = f"{allele_int}.{allele_dec}"
        return canon_allele

    @property
    def convert(self):
        sequence = self.forward_sequence
        if self.kit == "powerseq" or (len(sequence) % 4 != 0):
            final_seq = sequence_to_bracketed_form(sequence, self.repeat_size, self.repeats)
        else:
            final_seq = collapse_repeats_by_length(sequence, self.repeat_size)
        return final_seq


class STRMarker_DYS437(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        flank = f"{flank_seq[:3]} {collapse_repeats_by_length(flank_seq[3:], 4)}"
        return flank


class STRMarker_DYS392(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        flank = (
            f"{collapse_repeats_by_length(flank_seq[:72], 3)} "
            f"{collapse_repeats_by_length(flank_seq[72:], 3)}"
        )
        return flank


class STRMarker_DYS391(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:7], 4)} "
                f"{collapse_repeats_by_length(flank_seq[7:29], 4)} "
                f"{collapse_repeats_by_length(flank_seq[29:], 4)}"
            )
        else:
            flank = (
                f"{collapse_repeats_by_length(flank_seq[:7], 4)} "
                f"{collapse_repeats_by_length(flank_seq[7:], 4)}"
            )
        return flank

    @property
    def canonical(self):
        """Canonical STR allele designation"""
        n = self.repeat_size
        if self.kit == "powerseq":
            nsubout = self.data["BasesToSubtract"] - 6
        else:
            nsubout = self.data["BasesToSubtract"]
        nsubout *= -1
        new_seq = self.uas_sequence[:nsubout]
        if len(new_seq) % n == 0:
            canon_allele = int(len(new_seq) / n)
        else:
            allele_int = int(len(new_seq) / n)
            allele_dec = int(len(new_seq) % n)
            canon_allele = f"{allele_int}.{allele_dec}"
        return canon_allele

    @property
    def convert(self):
        sequence = self.forward_sequence
        if self.kit == "powerseq":
            final_seq = (
                f"{collapse_repeats_by_length_flanks(sequence[:6], 4)} "
                f"{collapse_repeats_by_length(sequence[6:], 4)}"
            )
        elif len(sequence) % 4 != 0:
            final_seq = sequence_to_bracketed_form(sequence, self.repeat_size, self.repeats)
        else:
            final_seq = collapse_repeats_by_length(sequence, self.repeat_size)
        return final_seq


class STRMarker_DYS19(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = (
                f"{collapse_repeats_by_length_flanks(flank_seq[:30], 4)} {flank_seq[30:32]} "
                f"{collapse_repeats_by_length(flank_seq[32:], 4)}"
            )
        else:
            flank = ""
        return flank


class STRMarker_DYS458(STRMarker):
    @property
    def convert(self):
        sequence = self.forward_sequence
        final_string = (
            f"{collapse_repeats_by_length(sequence[:14], 4)} "
            f"{collapse_repeats_by_length(sequence[14:], 4)}"
        )
        return final_string


class STRMarker_HPRTB(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        flank = f"{collapse_repeats_by_length(flank_seq[:-4], 4)} {flank_seq[-4]} {flank_seq[-3:]}"
        return flank


class STRMarker_DXS8378(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        flank = (
            f"{collapse_repeats_by_length(flank_seq[:21], 4)} "
            f"{collapse_repeats_by_length(flank_seq[21:47], 4)} {flank_seq[47]} "
            f"{collapse_repeats_by_length_flanks(flank_seq[48:-1], 4)} {flank_seq[-1]}"
        )
        return flank

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        flank = (
            f"{collapse_repeats_by_length(flank_seq[:26], 4)} {flank_seq[26]} "
            f"{collapse_repeats_by_length(flank_seq[27:], 4)}"
        )
        return flank


class STRMarker_DXS7132(STRMarker):
    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        flank = (
            f"{collapse_repeats_by_length(flank_seq[:5], 4)} "
            f"{collapse_repeats_by_length(flank_seq[5:36], 4)} {flank_seq[36]} "
            f"{collapse_repeats_by_length(flank_seq[37:], 4)}"
        )
        return flank


class STRMarker_DXS10135(STRMarker):
    @property
    def convert(self):
        sequence = self.forward_sequence
        final_string = (
            f"{collapse_repeats_by_length(sequence[:12], 4)} "
            f"{sequence[12:19].lower()} {collapse_repeats_by_length(sequence[19:], 4)}"
        )
        return final_string


class STRMarker_DXS10074(STRMarker):
    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        flank = (
            f"{collapse_repeats_by_length(flank_seq[:22], 4)} {flank_seq[22]} "
            f"{collapse_repeats_by_length(flank_seq[23:], 4)}"
        )
        return flank

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        flank = (
            f"{collapse_repeats_by_length(flank_seq[:5], 4)} "
            f"{collapse_repeats_by_length(flank_seq[:-3], 4)} {flank_seq[-2:]}"
        )
        return flank


class STRMarker_Y_GATA_H4(STRMarker):
    @property
    def convert(self):
        sequence = self.forward_sequence
        if self.kit == "powerseq":
            final_string = collapse_repeats_by_length(sequence, self.repeat_size)
        else:
            final_string = f"{sequence[0]} {collapse_repeats_by_length(sequence[1:], 4)}"
        return final_string

    @property
    def canonical(self):
        """Canonical STR allele designation"""
        n = self.repeat_size
        if self.software == "uas":
            nsubout = self.data["BasesToSubtract"]
        elif self.kit == "forenseq":
            nsubout = self.data["BasesToSubtract"] - 12
        else:
            nsubout = self.data["BasesToSubtract"] - 49
        nsubout *= -1
        new_seq = self.uas_sequence[:nsubout]
        if len(new_seq) % n == 0:
            canon_allele = int(len(new_seq) / n)
        else:
            allele_int = int(len(new_seq) / n)
            allele_dec = int(len(new_seq) % n)
            canon_allele = f"{allele_int}.{allele_dec}"
        return canon_allele


class STRMarker_DYS390(STRMarker):
    @property
    def canonical(self):
        """Canonical STR allele designation"""
        n = self.repeat_size
        if self.software == "uas" or self.kit == "powerseq":
            nsubout = self.data["BasesToSubtract"]
        else:
            nsubout = self.data["BasesToSubtract"] - 3
        nsubout *= -1
        new_seq = self.uas_sequence[:nsubout]
        if len(new_seq) % n == 0:
            canon_allele = int(len(new_seq) / n)
        else:
            allele_int = int(len(new_seq) / n)
            allele_dec = int(len(new_seq) % n)
            canon_allele = f"{allele_int}.{allele_dec}"
        return canon_allele

    @property
    def designation(self):
        lus, sec, ter = None, None, None
        lus = repeat_copy_number(self.convert, self.data["LUS"])
        sec = repeat_copy_number(self.convert, self.data["Sec"])
        if self.software == "uas" or self.kit == "powerseq":
            ter = repeat_copy_number(self.convert, self.data["Tert"])
        else:
            if self.convert[-1] == "G":
                ter = "1"
            else:
                ter = "0"
        return lus, sec, ter

    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        flank = f"{flank_seq[:2]} {collapse_repeats_by_length(flank_seq[2:], 4)}"
        return flank


class STRMarker_DYS385(STRMarker):
    @property
    def canonical(self):
        """Canonical STR allele designation"""
        n = self.repeat_size
        if self.software == "uas" or self.kit == "forenseq":
            nsubout = self.data["BasesToSubtract"]
        else:
            nsubout = self.data["BasesToSubtract"] - 4
        nsubout *= -1
        new_seq = self.uas_sequence[:nsubout]
        if len(new_seq) % n == 0:
            canon_allele = int(len(new_seq) / n)
        else:
            allele_int = int(len(new_seq) / n)
            allele_dec = int(len(new_seq) % n)
            canon_allele = f"{allele_int}.{allele_dec}"
        return canon_allele


class STRMarker_DYS448(STRMarker):
    @property
    def designation(self):
        lus, sec, ter = None, None, None
        anno = self.convert
        repeat = "AGAGAT"
        match_list = []
        for block in anno.split(" "):
            if block == repeat:
                match_list.append(1)
            match = re.match(r"\[" + repeat + r"\](\d+)", block)
            if match:
                length = int(match.group(1))
                match_list.append(length)
        if len(match_list) == 0:
            lus, sec, ter = 0, 0, 0
        elif len(match_list) == 1:
            lus = match_list[0]
            sec = 0
        else:
            lus = match_list[0]
            sec = match_list[-1]
        return lus, sec, ter

    @property
    def flank_3p(self):
        flank_seq = self.flankseq_3p
        if self.kit == "powerseq":
            flank = f"{flank_seq[:5]} {collapse_repeats_by_length(flank_seq[5:], 6)}"
        else:
            flank = flank_seq
        return flank


class STRMarker_DXS10103(STRMarker):
    @property
    def designation(self):
        """
        The LUS and Secondary motif allele are the same, "TAGA".
        The Secondary motif is always the first "TAGA" in the sequence.
        The LUS, however, is not as easily identified given there could be multiple instances
        of "TAGA" within the sequence. After the Secondary motif is identified, the LUS allele
        is identified as the "TAGA" repeat sequence with the largest number of repeats.
        """
        lus, sec, ter = None, None, None
        anno = self.convert
        repeat = "TAGA"
        match_list = []
        for block in anno.split(" "):
            if block == repeat:
                match_list.append(1)
            match = re.match(r"\[" + repeat + r"\](\d+)", block)
            if match:
                length = int(match.group(1))
                match_list.append(length)
        if len(match_list) == 1:
            lus = match_list[0]
            sec = 0
        elif len(match_list) == 2:
            lus = match_list[-1]
            sec = match_list[0]
        elif len(match_list) > 2:
            lus = max(match_list)
            sec = match_list[0]
        else:
            pass
        return lus, sec, ter


class STRMarker_DYS389II(STRMarker):
    @property
    def designation(self):
        """
        The LUS and Secondary motif allele are the same, "TAGA".
        The Secondary motif is always the first "TAGA" in the sequence.
        The LUS, however, is not as easily identified given there could be multiple instances
        of "TAGA" within the sequence. After the Secondary motif is identified, the LUS allele
        is identified as the "TAGA" repeat sequence with the largest number of repeats.
        """
        lus, sec, ter = None, None, None
        anno = self.convert
        repeat = "TAGA"
        match_list = []
        for block in anno.split(" "):
            if block == repeat:
                match_list.append(1)
            match = re.match(r"\[" + repeat + r"\](\d+)", block)
            if match:
                length = int(match.group(1))
                match_list.append(length)
        if len(match_list) == 1:
            lus = match_list[0]
            sec = 0
        elif len(match_list) == 2:
            lus = match_list[-1]
            sec = match_list[0]
        elif len(match_list) > 2:
            lus = max(match_list)
            sec = match_list[0]
        else:
            pass
        return lus, sec, ter

    @property
    def flank_5p(self):
        flank_seq = self.flankseq_5p
        if self.kit == "powerseq":
            flank = f"{flank_seq[:3]} {collapse_repeats_by_length(flank_seq[3:], 4)}"
        else:
            flank = ""
        return flank


def STRMarkerObject(locus, sequence, software, custom=False, kit="forenseq"):
    constructors = {
        "D8S1179": STRMarker_D8S1179,
        "D13S317": STRMarker_D13S317,
        "D20S482": STRMarker_D20S482,
        "D2S441": STRMarker_D2S441,
        "D7S820": STRMarker_D7S820,
        "D16S539": STRMarker_D16S539,
        "D1S1656": STRMarker_D1S1656,
        "PENTA D": STRMarker_PentaD,
        "VWA": STRMarker_vWA,
        "D10S1248": STRMarker_D10S1248,
        "D22S1045": STRMarker_D22S1045,
        "FGA": STRMarker_FGA,
        "CSF1PO": STRMarker_CSF1PO,
        "D18S51": STRMarker_D18S51,
        "D21S11": STRMarker_D21S11,
        "TH01": STRMarker_TH01,
        "D19S433": STRMarker_D19S433,
        "D2S1338": STRMarker_D2S1338,
        "D5S818": STRMarker_D5S818,
        "PENTA E": STRMarker_PentaE,
        "TPOX": STRMarker_TPOX,
        "DYS643": STRMarker_DYS643,
        "DYS635": STRMarker_DYS635,
        "DYS612": STRMarker_DYS612,
        "DYS576": STRMarker_DYS576,
        "DYS549": STRMarker_DYS549,
        "DYS533": STRMarker_DYS533,
        "DYS522": STRMarker_DYS522,
        "DYS439": STRMarker_DYS439,
        "DYS437": STRMarker_DYS437,
        "DYS392": STRMarker_DYS392,
        "DYS391": STRMarker_DYS391,
        "DYS19": STRMarker_DYS19,
        "DYS458": STRMarker_DYS458,
        "HPRTB": STRMarker_HPRTB,
        "DXS8378": STRMarker_DXS8378,
        "DXS7132": STRMarker_DXS7132,
        "DXS10135": STRMarker_DXS10135,
        "DXS10074": STRMarker_DXS10074,
        "Y-GATA-H4": STRMarker_Y_GATA_H4,
        "DYS390": STRMarker_DYS390,
        "DYS385A-B": STRMarker_DYS385,
        "DYS448": STRMarker_DYS448,
        "DXS10103": STRMarker_DXS10103,
        "DYS389II": STRMarker_DYS389II,
    }
    if locus in constructors:
        constructor = constructors[locus]
        return constructor(locus, sequence, software=software, custom=custom, kit=kit)
    else:
        return STRMarker(locus, sequence, software=software, custom=custom, kit=kit)
