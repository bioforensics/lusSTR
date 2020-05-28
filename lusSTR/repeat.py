#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2020, Battelle National Biodefense Institute.
#
# This file is part of lusSTR (http://github.com/bioforensics/lusSTR)
# and is licensed under the BSD license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import re


def collapse_tandem_repeat(fullseq, repeat):
    '''Collapse tandem stretches of the specified repeat sequence in a larger sequence.

    >>> collapse_tandem_repeat('TAGATTATTATTTAGTAGATTTAGTAG', 'ATT')
    'TAG [ATT]3 TAGTAG ATT TAGTAG'
    >>> collapse_tandem_repeat('TAGATTATTATTTAGTAGATTTAGTAG', 'TAG')
    'TAG ATTATTATT [TAG]2 ATT [TAG]2'
    '''
    if repeat not in fullseq:
        return fullseq
    i = fullseq.find(repeat)
    prefix = fullseq[:i]
    suffix = fullseq[i:]
    count = 0
    while suffix.startswith(repeat):
        count += 1
        suffix = suffix[len(repeat):]
    if count == 1:
        formatted = f' {repeat} '
    else:
        formatted = f' [{repeat}]{count} '
    final = prefix + formatted + collapse_tandem_repeat(suffix, repeat)
    final = final.strip()
    final = re.sub(r' +', ' ', final)
    return final


def collapse_all_repeats(sequence, repeats):
    '''Convert a sequence to bracketed form by collapsing stretches of tandem repeats.

    >>> collapse_all_repeats('TAGATTATTATTTAGTAGATTTAGTAG', ['ATT', 'TAG'])
    'TAG [ATT]3 [TAG]2 ATT [TAG]2'
    '''
    collapsed_seq = sequence
    for repeat in repeats:
        collapsed_seq = collapse_tandem_repeat(collapsed_seq, repeat)
    return collapsed_seq


def split_by_n(sequence, n):
    '''Split a sequence into non-overlapping chunks of length n.'''
    while sequence:
        yield sequence[:n]
        sequence = sequence[n:]


def get_blocks(sequence, n):
    '''Split a sequence into chunks of length n, and count adjacent repeated chunks.'''
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


def collapse_repeats_by_length(sequence, n):
    '''Convert to bracketed annotation form by splitting the sequence into blocks of size n.'''
    units = list()
    for unit, count in get_blocks(sequence, n):
        assert unit is not None, (sequence, n)
        if count == 1:
            units.append(unit)
        else:
            units.append(f'[{unit}]{count}')
    result = '  '.join(units)
    result = re.sub(r' +', ' ', result)
    return result


def sequence_to_bracketed_form(sequence, n, repeats):
    '''Convert sequence to bracketed annotation.

    Uses a combination of repeat-based and length-based methods to convert a sequence containing
    tandem repeats into a concise bracketed representation.
    '''
    collapsed = collapse_all_repeats(sequence, repeats)
    blocks = list()
    for unit in collapsed.split(' '):
        if len(unit) > n and '[' not in unit:
            for x in split_by_n(unit, n):
                blocks.append(x)
        else:
            blocks.append(unit)
    result = ' '.join(blocks)
    result = re.sub(r' +', ' ', result)
    return result


def reverse_complement(sequence):
    '''
    Function creates reverse complement of sequence

    Sequences in which the UAS software output contains the sequence on the reverse strand
    require translation of the sequence to the forward strand. This allows for consistency
    between both loci and any outside analyses in which comparisons may be made.
    '''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    rclist = [complement[base] for base in sequence[::-1]]
    rc = ''.join(rclist)
    return rc


def reverse_complement_bracketed(forward_bracket):
    '''Compute reverse complement of a bracketed form annotation.'''
    inblocks = forward_bracket.split(' ')
    outblocks = list()
    for block in reversed(inblocks):
        match = re.match(r'\[([ACGT]+)\](\d+)', block)
        if match:
            rcrep = reverse_complement(match.group(1))
            count = match.group(2)
            rcblock = f'[{rcrep}]{count}'
        else:
            if re.match(r'[^ACGT]', block):
                raise ValueError(f'annotation block "{block}" includes invalid characters')
            rcblock = reverse_complement(block)
        outblocks.append(rcblock)
    return ' '.join(outblocks)


def repeat_copy_number(bf, repeat):
    '''Determine the longest uninterrupted stretch of the specified repeat.

    The input is a sequence string collapsed to bracketed annotation form.
    '''
    longest = 0
    for block in bf.split(' '):
        if block == repeat:
            if 1 > longest:
                longest = 1
        match = re.match(r'\[' + repeat + r'\](\d+)', block)
        if match:
            length = int(match.group(1))
            if length > longest:
                longest = length
    return str(longest)
