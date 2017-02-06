#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import re
from collections import defaultdict

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

class ReferenceAlignment(object):
    """
        Examples:
        >>> aln = ReferenceAlignment('ac-gtac', 'accg--c')
        >>> aln.rseq()
        'acgtac'
        >>> aln.qseq()
        'accgc'
    """
    gapchar = set('.-')
    unkchar = set('?')
    
    def __init__(self, raln=None, qaln=None):
        if raln is not None and qaln is not None:
            self.load_alignment(raln, qaln)
        else:
            self.aln_positions = []
            self._rpos_to_apos = {}
            self._qpos_to_apos = {}
            self.alen = 0
            self.rstart = self.rend = 0
            self.qstart = self.qend = 0
    
    def load_alignment(self, raln, qaln):
        assert len(raln) == len(qaln), 'Error: Alignments are different lengths'
        self.aln_positions = []
        
        # Convert all gapchars to '.'
        for rb, qb in zip(raln, qaln):
            rb = '.' if rb in self.gapchar else rb
            qb = '.' if qb in self.gapchar else qb
            self.aln_positions.append((rb, qb))
        
        self.alen = len(self.aln_positions)
        self._index_alignment()
    
    def _index_alignment(self):
        if self.alen != len(self.aln_positions):
            # print >>sys.stderr, 'Warning: self.alen != len(self.aln_positions)'
            self.alen = len(self.aln_positions)
        
        # Initialize
        self._rpos_to_apos = {} # one-to-one with some apos not included
        self._qpos_to_apos = {} # one-to-one with some apos not included
        self._apos_to_rpos = defaultdict(list)
        self._apos_to_qpos = defaultdict(list)
        rpos = qpos = 0
        for apos in range(self.alen):
            rb, qb = self.aln_positions[apos]
            if rb is not '.':
                self._rpos_to_apos[rpos] = apos
                rpos += 1
            if qb is not '.':
                self._qpos_to_apos[qpos] = apos
                qpos += 1
        
        self.rstart, self.rend = 0, rpos
        self.qstart, self.qend = 0, qpos
    
    def adjust_ref_start(self, refstart):
        new_rpos_to_apos = {}
        adj = refstart - self.rstart
        for i in range(self.rstart, self.rend):
            new_rpos_to_apos[i + adj] = self._rpos_to_apos[i]
        self._rpos_to_apos = new_rpos_to_apos
        self.rstart += adj
        self.rend += adj
        
    def merge_alignments(self, other):
        newaln = ReferenceAlignment()
        
        # First and last reference positions in other
        other_rs = min(other._rpos_to_apos.keys())
        other_re = max(other._rpos_to_apos.keys())
        
        if other_rs in self._rpos_to_apos:
            left = self.aln_positions[:self._rpos_to_apos[other_rs]]
        if other_re in self._rpos_to_apos:
            right = self.aln_positions[self._rpos_to_apos[other_re] + 1:]
        
        newaln = ReferenceAlignment()
        newaln.aln_positions = left + other.aln_positions + right
        newaln._index_alignment()
        return newaln
    
    def rseq(self):
        return ''.join(self.aln_positions[self._rpos_to_apos[i]][0]
           for i in range(self.rstart, self.rend))
    
    def raln(self):
        return ''.join(t[0] for t in self.aln_positions)
    
    def qseq(self):
        return ''.join(self.aln_positions[self._qpos_to_apos[i]][1]
            for i in range(self.qstart, self.qend))
    def qaln(self):
        return ''.join(t[1] for t in self.aln_positions)

    def imputed(self):
        ret = ''
        for i in range(self.qstart, self.qend):
            rb, qb = self.aln_positions[self._qpos_to_apos[i]]
            ret += rb.lower() if qb == '?' else qb.upper()
        return ret
    
    def scaffold(self):
        return self.qseq().upper().replace('?', 'n')
    
    def padded(self):
        return self.qseq().upper().replace('?', '.')
    
    def convert_rpos(self, rpos):
        if self.rstart <= rpos <= self.rend:
            alntarget = self._rpos_to_apos[rpos]
            _apos_to_qpos = {ap:qp for qp,ap in self._qpos_to_apos.iteritems()}
            while True:
                if alntarget in _apos_to_qpos:
                    return alntarget[rev]
                else:
                    alntarget -= 1
        return 0

class EmptyReferenceAlignment(ReferenceAlignment):
    def __init__(self, refseq):
        super(EmptyReferenceAlignment, self).__init__()
        for p, rb in enumerate(refseq):
            self._qpos_to_apos[p] = p
            self._rpos_to_apos[p] = p
            self.aln_positions.append((rb, '?'))
        
        self.alen = len(self.aln_positions)
        self.rstart, self.rend = 0, len(self.aln_positions)
        self.qstart, self.qend = 0, len(self.aln_positions)

class NucmerReferenceAlignment(ReferenceAlignment):
    def __init__(self, outlines=None):
        super(NucmerReferenceAlignment, self).__init__()
        self.ref_frm = None
        self.ref_s = None
        self.ref_e = None
        self.qry_frm = None
        self.qry_s = None
        self.qry_e = None
        
        if outlines is not None:
            self.parse(outlines)
    
    def parse(self, outlines):
        aln_lines = []
        flag = False
        for l in outlines:
            if l.startswith('-- BEGIN'):
                m = re.search('\[([+\-\s\d]+)\|([+\-\s\d]+)\]', l)
                ref_frm, ref_s, z ,ref_e = m.group(1).strip().split()
                qry_frm, qry_s, z ,qry_e = m.group(2).strip().split()
                flag = True
            elif l.startswith('--   END'):
                flag = False
            if flag and re.match('^\d+', l):
                aln_lines.append(l)
        
        self.ref_frm, self.ref_s, self.ref_e = map(int, [ref_frm, ref_s, ref_e])
        self.qry_frm, self.qry_s, self.qry_e = map(int, [qry_frm, qry_s, qry_e])
        
        self.load_alignment(
            ''.join(aln_lines[i].split()[1] for i in range(0, len(aln_lines), 2)),
            ''.join(aln_lines[i].split()[1] for i in range(1, len(aln_lines), 2))
        )
        # Adjust for reference start
        # Nucmer uses 1-based numbering, we are going to use 0-based numbering
        self.adjust_ref_start(self.ref_s - 1)



class TilingRow:
    """ Alignment specification from show-tiling (mummer)
    """
    def __init__(self, l, seqdict=None):
        fields = l.strip('\n').split('\t')
        self.ref_s, self.ref_e, self.qry_s, self.qry_e = map(int, fields[:4])
        self.ref_alen, self.qry_alen = map(int, fields[4:6])
        self.pid = float(fields[6])
        self.ref_len, self.qry_len = map(int, fields[7:9])
        self.ref_cov, self.qry_cov = map(float, fields[9:11])
        self.ref, self.qry = fields[11:]
        if seqdict is not None:
            self.add_seq(seqdict)
    
    def __str__(self):
        return '\t'.join(map(str, [self.ref_s, self.ref_e, self.qry_s, self.qry_e,
                                   self.ref_alen, self.qry_alen, self.pid,
                                   self.ref_len, self.qry_len,
                                   self.ref_cov, self.qry_cov,
                                   self.ref, self.qry]))

