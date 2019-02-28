# -*- coding: utf-8 -*-
"""Utilities for working with seqeunces or sets of sequences
"""
from __future__ import print_function
from __future__ import absolute_import

import sys
import re
import os

from haphpipe.utils.helpers import merge_interval_list


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


class HPSeq(object):
    def __init__(self, name, seq):
        self.id = name
        self.seq = seq

    def __len__(self):
        return len(self.seq)


class HPSeqIO(object):
    @staticmethod
    def parse(filename, format):
        assert format == 'fasta'
        with open(filename, 'rU') as fh:
            for n,s in fastagen(fh):
                yield HPSeq(n,s)

    @staticmethod
    def read(filename, format):
        assert format == 'fasta'
        with open(filename, 'rU') as fh:
            n, s = list(fastagen(fh))[0]
            return HPSeq(n, s)


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


def wrap(s, wraplen=60):
    return '\n'.join(s[i:(i+wraplen)] for i in range(0, len(s), wraplen))


def clean_seqnames(fh):
    """ Remove sequence description

    Some programs do not play nicely with FASTA description lines that include
    whitespace. Strip off anything after first whitespace

    Args:
        fh (file): Filehandle of FASTA file

    Yields:
        seq_id, seq (2-tuple of str): The next sequence ID and sequence pair

    """
    for i, tup in enumerate(fastagen(fh)):
        yield ('%s' % tup[0].split()[0], tup[1])


def clean_seqnames_file(infile, outdir):
    """ Remove sequence description from file

    Some programs do not play nicely with FASTA description lines that include
    whitespace. Strip off anything after first whitespace. See `clean_seqnames`

    Args:
        infile (str): Path to input FASTA file
        outdir (str): Path to output directory

    Returns:
        out_clean (str): Path to clean FASTA file

    """
    out_clean = os.path.join(outdir, 'tmp_clean_seqnames.fa')
    with open(out_clean, 'w') as outh, open(infile, 'rU') as fh:
        for n,s in clean_seqnames(fh):
            print('>%s\n%s' % (n, wrap(s)), file=outh)
    return out_clean


def N50(l):
    slen = sorted(l)
    cumsum = [sum(slen[:i+1]) for i in range(len(slen))]
    i = len([c for c in cumsum if c < sum(slen)*0.5])
    return len(slen)-i, slen[i]


def assembly_stats(fh, outh=sys.stdout):
    contig_lengths = sorted(len(s) for n,s in fastagen(fh))
    print('num_contigs\t%d' % len(contig_lengths), file=outh)
    print('max_contig\t%d' % contig_lengths[-1], file=outh)
    print('contig>1kb\t%d' % sum(l>=1000 for l in contig_lengths), file=outh)
    l50,n50 = N50(contig_lengths)
    print('contig_N50\t%d' % n50, file=outh)
    print('contig_L50\t%d' % l50, file=outh)
    contig_lengths.sort(reverse=True)
    for i,l in enumerate(contig_lengths[1:5]):
        print('rank%d\t%d' % (i+2, l), file=outh)

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

def region_to_tuple(regstr):
    chrom = regstr.split(':')[0]
    reg_s, reg_e = map(int, regstr.split(':')[1].split('-'))
    return chrom, reg_s, reg_e

def parse_seq_id(s, delim='|'):
    f = [_ for _ in s.split(delim) if _]
    return {f[i]:f[i+1] for i in range(0, len(f), 2)}


def make_seq_id(**kwargs):
    """ Return sequence ID in consistent format"""
    ret = ''
    keyorder = ['sid', 'ref', 'reg',]
    for k in keyorder:
        if k in kwargs:
            ret += '%s|%s|' % (k, kwargs.pop(k))
    for k,v in kwargs.items():
        ret += '%s|%s|' % (k, v)
    return ret

def update_seq_id(seq_id, sample_id):
    if sample_id == 'sampleXX':
        return seq_id
    d = parse_seq_id(seq_id)
    d['sid'] = sample_id
    return make_seq_id(**d)

"""
class BedLine(object):
    BEDFIELDS = [
        ('chrom', str),
        ('chromStart', int),
        ('chromEnd', int),
        ('name', str),
        ('score', lambda x: None if x=='.' else float(x)),
        ('strand', str),
        ('thickStart', int),
        ('thickEnd', int),
        ('itemRgb', lambda x: tuple(map(int, x.split(',')))),
        ('blockCount', int),
        ('blockSizes', lambda x: map(int, x.split(','))),
        ('blockStarts', lambda x: map(int, x.split(','))),
    ]
    
    def __init__(self, l):
        for v, (n, f) in zip(l, self.BEDFIELDS[:len(l)]):
            setattr(self, n, f(v))
    
    def __str__(self):
        ret = []
        for n,f in self.BEDFIELDS:
           if hasattr(self, n):
               v = getattr(self, n)
               if v is None:
                   ret.append('.')
               elif type(v) is list or type(v) is tuple:
                   ret.append(','.join(str(_) for _ in v))
               else:
                   ret.append(str(v))
        return '\t'.join(ret)
        
def bed_parser(infile):
    fh = open(infile, 'rU') if type(infile) is str else infile
    bedlines = (l.strip('\n').split('\t') for l in fh)
    for bl in bedlines:
        yield BedLine(bl)
"""
