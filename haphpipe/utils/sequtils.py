# -*- coding: utf-8 -*-
"""Utilities for working with seqeunces or sets of sequences
"""

import sys
import re

from helpers import merge_interval_list

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def fastagen(fh):
    lines = (l.strip() for l in fh)
    name, seq = None, ''
    for l in lines:
        if not l: continue
        if l.startswith('>'):
            if name is not None:
                yield name, seq
            name = l.lstrip('>')
            seq = ''
        else:
            seq = seq + l
    yield name, seq

def clean_seqnames(fh):
    for i, tup in enumerate(fastagen(fh)):
        yield ('%s' % tup[0].split()[0], tup[1])

def wrap(s, wraplen=60):
    return '\n'.join(s[i:(i+wraplen)] for i in range(0, len(s), wraplen))

def N50(l):
    slen = sorted(l)
    cumsum = [sum(slen[:i+1]) for i in range(len(slen))]
    i = len([c for c in cumsum if c < sum(slen)*0.5])
    return len(slen)-i, slen[i]

def assembly_stats(fh, outh=sys.stdout):
    contig_lengths = sorted(len(s) for n,s in fastagen(fh))
    print >>outh, 'num_contigs\t%d' % len(contig_lengths)
    print >>outh, 'max_contig\t%d' % contig_lengths[-1]
    print >>outh, 'contig>1kb\t%d' % sum(l>=1000 for l in contig_lengths)
    l50,n50 = N50(contig_lengths)
    print >>outh, 'contig_N50\t%d' % n50
    print >>outh, 'contig_L50\t%d' % l50
    contig_lengths.sort(reverse=True)
    for i,l in enumerate(contig_lengths[1:5]):
        print >>outh, 'rank%d\t%d' % (i+2, l)

def unambig_intervals(s, maxgap=30):
    """ Return unambiguous intervals of a sequence string
    
    """
    tmp = [(m.start(), m.end()) for m in re.finditer('[ACGTacgt]+', s)]
    return merge_interval_list(tmp, maxgap)

def extract_amplicons(name, seq, maxgap=30):
    ivs = unambig_intervals(seq, maxgap)
    for iv in ivs:
        amp_name = '%s:%d-%d' % (name, iv[0]+1, iv[1])
        amp = seq[iv[0]:iv[1]]
        yield amp_name, amp

'''
def extract_amplicons_fh(fh, mingap=30):
    amplicons = []
    for name, seq in fastagen(fh):
        for aname, aseq in extract_amplicons(name, seq):
            amplicons.append((aname, aseq))
    return amplicons
'''

AMBIG_MAP = {'AC': 'M',
             'AG': 'R',
             'AT': 'W',
             'CG': 'S',
             'CT': 'Y',
             'GT':'K',
             'ACG': 'V',
             'ACT': 'H',
             'AGT': 'D',
             'CGT': 'B',
             'ACGT': 'N'
             }

def get_ambig(bases):
    return AMBIG_MAP[''.join(sorted(set(bases)))]             