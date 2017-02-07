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
    
    def __init__(self, raln=None, qaln=None, rname=None, qname=None):
        if raln is not None and qaln is not None:
            self.load_alignment(raln, qaln)
        else:
            self.aln_positions = []
            self._rpos_to_apos = {}
            self._qpos_to_apos = {}
            self.alen = 0
            self.rstart = self.rend = 0
            self.qstart = self.qend = 0
            self.rname = rname
            self.qname = qname
    
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
        print >>sys.stderr, "Indexing alignment"
        if self.alen != len(self.aln_positions):
            # print >>sys.stderr, 'Warning: self.alen != len(self.aln_positions)'
            self.alen = len(self.aln_positions)
        
        # Initialize
        self._rpos_to_apos = {} # one-to-one with some apos not included
        self._qpos_to_apos = {} # one-to-one with some apos not included
        self._apos_to_rpos = {}
        self._apos_to_qpos = {}
        rpos = qpos = 0
        for apos in range(self.alen):
            rb, qb = self.aln_positions[apos]
            if rb is not '.':
                self._rpos_to_apos[rpos] = apos
                self._apos_to_rpos[apos] = rpos
                rpos += 1
            if qb is not '.':
                self._qpos_to_apos[qpos] = apos
                self._apos_to_qpos[apos] = qpos
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
    
    def convert_rpos(self, rp, left=True):
        """ Returns position in query for given reference position
        """
        incr = -1 if left else 1
        if rp < self.rstart:
            print >>sys.stderr, "WARNING: position %d is outside reference boundaries" % rp
            return self.qstart
        elif rp >= self.rend:
            print >>sys.stderr, "WARNING: position %d is outside reference boundaries" % rp
            return self.qend
        else:
            ap = self._rpos_to_apos[rp]
            while 0 <= ap <= self.alen:
                if ap in self._apos_to_qpos:
                    return self._apos_to_qpos[ap]
                ap += incr
            print >>sys.stderr, "WARNING: position %d is outside alignment" % ap
            return None
    
    def merge_alignments(self, other):
        newaln = ReferenceAlignment()
        
        # First and last reference positions in other
        other_rs = min(other._rpos_to_apos.keys())
        assert other_rs == other.rstart
        other_re = max(other._rpos_to_apos.keys())
        assert other_re + 1 == other.rend, "%d != %d" % (other_re, other.rend)
        # print other_rs
        # print other_re
        # print self._rpos_to_apos[other_rs]
        # print self._rpos_to_apos[other_re] + 1
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

# aln = ReferenceAlignment('at.gtacc', 'atcg..cc')
# aln.adjust_ref_start(200)
# aln.convert_rpos(200)

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
    
    def parse(self, aln_report):
        assert re.match('^--\s+BEGIN', aln_report[0]), 'First line of alignment is not "BEGIN"'
        assert re.match('^--\s+END', aln_report[-1]), 'Last line of alignment is not "END"'
        m = re.search('\[([+\-\s\d]+)\|([+\-\s\d]+)\]', aln_report[0])
        ref_frm, ref_s, z ,ref_e = m.group(1).strip().split()
        qry_frm, qry_s, z ,qry_e = m.group(2).strip().split()
        """
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
        """
        self.ref_frm, self.ref_s, self.ref_e = map(int, [ref_frm, ref_s, ref_e])
        self.qry_frm, self.qry_s, self.qry_e = map(int, [qry_frm, qry_s, qry_e])
        aln_lines = [l for l in aln_report if re.match('^\d+', l)]
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

