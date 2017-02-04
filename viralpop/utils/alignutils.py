#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
from subprocess import check_output
from collections import defaultdict, namedtuple
from itertools import chain

from Bio import SeqIO

from utils.sysutils import command_runner
from utils.helpers import overlaps

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

class TilingRow:
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
    
    def add_seq(self, seqdict):
        if self.qry_s < self.qry_e:
            self.seq = seqdict[self.qry].seq[(self.qry_s - 1): self.qry_e]
        else:
            self.seq = seqdict[self.qry].seq[(self.qry_e - 1): self.qry_s].reverse_complement()
    
    def ref_align(self, delta_file, quiet=True):
        out1 = show_aligns(self.ref, self.qry, delta_file)
        if not quiet:
            print out1
        self.alignment = parse_show_aligns(out1)
        self.alt_alignment = NucmerAlignment(out1.strip('\n').split('\n'))
    
    def __str__(self):
        return '\t'.join(map(str, [self.ref_s, self.ref_e, self.qry_s, self.qry_e,
                                   self.ref_alen, self.qry_alen, self.pid,
                                   self.ref_len, self.qry_len,
                                   self.ref_cov, self.qry_cov,
                                   self.ref, self.qry]))


class AlignmentPosition:
    def __init__(self, ap, rp, rb, qp, qb):
        self.aln_pos = ap
        self.ref_pos = rp
        self.ref_base = rb
        self.qry_pos = qp
        self.qry_base = qb
    
    def __repr__(self):
        'Return a nicely formatted representation string'
        ret = 'AlignmentPosition(aln_pos=%r, ref_pos=%r, ref_base=%r, qry_pos=%r, qry_base=%r)'
        return ret % (self.aln_pos, self.ref_pos, self.ref_base, self.qry_pos, self.qry_base)

class GenericAlignment(object):
    def __init__(self):
        self.aln_positions = None
        self.gapchar = set('.-')
        self.unkchar = set('?')
        self._mapping_rpos_to_apos = None # reference positions (keys) map to alignment positions (values)
    
    def _make_mapping_rpos_to_apos(self):
        self._mapping_rpos_to_apos = {} # reference positions (keys) map to alignment positions (values)
        for ap in self.aln_positions:
            if ap.ref_pos is None: continue
            assert ap.ref_pos not in self._mapping_rpos_to_apos, 'ERROR %d' % ap.ref_pos
            self._mapping_rpos_to_apos[ap.ref_pos] = ap.aln_pos
    
    def _reindex_alignment(self):
        # Alignment positions are numbered sequentially starting with 0
        # Query positions are number sequentially for non-gap positions        
        cur_qpos = 0
        for i, ap in enumerate(self.aln_positions):
            ap.aln_pos = i
            if ap.qry_base not in self.gapchar:
                ap.qry_pos = cur_qpos
                cur_qpos += 1
            else:
                assert ap.qry_pos is None
        
        self._make_mapping_rpos_to_apos()
    
    def load_alignment_string(self, raln, qaln):
        aln_len = len(raln)
        assert aln_len == len(qaln), 'ERROR: Alignment strings are unequal:\n%s\n%s' % (raln, qaln)
        
        self.aln_positions = []
        rpos = qpos = 0
        for aln_pos, rb, qb in zip(range(aln_len), raln, qaln):
            if rb in self.gapchar and qb in self.gapchar:
                self.aln_positions.append(AlignmentPosition(aln_pos, None, '.', None, '.'))
            elif rb in self.gapchar:
                self.aln_positions.append(AlignmentPosition(aln_pos, None, '.', qpos, qb))
            elif qb in self.gapchar:
                self.aln_positions.append(AlignmentPosition(aln_pos, rpos, rb, None, '.'))
            else:
                self.aln_positions.append(AlignmentPosition(aln_pos, rpos, rb, qpos, qb))
            
            rpos += 0 if rb in self.gapchar else 1
            qpos += 0 if qb in self.gapchar else 1
        
        self._make_mapping_rpos_to_apos()
    
    def combine_alignments(self, other):
        newaln = GenericAlignment()
        
        # Segment before other
        oref_s = other.aln_positions[0].ref_pos
        # Find corresponding alignment position in self
        apos_s = self._mapping_rpos_to_apos[oref_s]
        seg0 = self.aln_positions[:apos_s]
        # print 'ref segment 0: %d-%d' % (seg0[0].ref_pos, seg0[-1].ref_pos)
        
        # Middle segment
        seg1 = other.aln_positions
        # print 'new segment 1: %d-%d' % (seg1[0].ref_pos, seg1[-1].ref_pos)
        
        # Segment after other
        oref_e = other.aln_positions[-1].ref_pos
        # Find corresponding alignment position in self
        apos_e = self._mapping_rpos_to_apos[oref_e]
        seg2 = self.aln_positions[apos_e + 1:]
        # print 'ref segment 2: %d-%d' % (seg2[0].ref_pos, seg2[-1].ref_pos)
        
        newaln.aln_positions = seg0 + seg1 + seg2
        newaln._reindex_alignment()
        return newaln
        
    def get_aligned_str(self):
        raln = ''.join(ap.ref_base for ap in self.aln_positions)
        qaln = ''.join(ap.qry_base for ap in self.aln_positions)
        return raln, qaln
    
    def get_scaffold(self, ref_s=None, ref_e=None):
        if ref_s is None:
            ap_s = self.aln_positions[0].aln_pos
        else:
            ap_s = self._mapping_rpos_to_apos[ref_s]
        
        if ref_e is None:
            ap_e = self.aln_positions[-1].aln_pos
        else:
            ap_e = self._mapping_rpos_to_apos[ref_e]
        
        scaf_str = ''
        for apos in range(ap_s, ap_e+1):
            if self.aln_positions[apos].qry_base in self.unkchar:
                scaf_str += 'n'
            elif self.aln_positions[apos].qry_base in self.gapchar:
                scaf_str += ''
            else:
                scaf_str += self.aln_positions[apos].qry_base.upper()
        return scaf_str

    def get_imputed(self, ref_s=None, ref_e=None):
        if ref_s is None:
            ap_s = self.aln_positions[0].aln_pos
        else:
            ap_s = self._mapping_rpos_to_apos[ref_s]
        
        if ref_e is None:
            ap_e = self.aln_positions[-1].aln_pos
        else:
            ap_e = self._mapping_rpos_to_apos[ref_e]
        
        imp_str = ''
        for apos in range(ap_s, ap_e+1):
            if self.aln_positions[apos].qry_base in self.unkchar:
                imp_str += self.aln_positions[apos].ref_base.lower()
            elif self.aln_positions[apos].qry_base in self.gapchar:
                imp_str += ''
            else:
                imp_str += self.aln_positions[apos].qry_base.upper()
        return imp_str
    
    def as_matrix(self):
        ret = []
        for ap in self.aln_positions:
            ret.append([
                '%d' % ap.aln_pos,
                '' if ap.ref_pos is None else '%d' % ap.ref_pos,
                ap.ref_base,
                '' if ap.qry_pos is None else '%d' % ap.qry_pos,
                ap.qry_base,
            ])
        return ret
    
    def from_matrix(self, matrix):
        self.aln_positions = []
        for row in matrix:
            self.aln_positions.append(        
                AlignmentPosition(
                    int(row[0]),
                    None if row[1] == '' else int(row[1]),
                    row[2],
                    None if row[3] == '' else int(row[3]),
                    row[4],
                )
            )
        self._reindex_alignment()

class EmptyAlignment(GenericAlignment):
    def __init__(self, refseq):
        super(EmptyAlignment, self).__init__()
        self.aln_positions = []
        for i,b in enumerate(refseq):
            self.aln_positions.append(AlignmentPosition(i, i, b, i, '?'))
        # Create mapping dictionaries
        self._make_mapping_rpos_to_apos()

class NucmerAlignment(GenericAlignment):
    def __init__(self, out=None):
        super(NucmerAlignment, self).__init__()
        self.ref_frm = None
        self.ref_s = None
        self.ref_e = None
        self.qry_frm = None
        self.qry_s = None
        self.qry_e = None
        
        self.raln = ''
        self.qaln = ''
        
        if out is not None:
            self.parse(out)
    
    def parse(self, out):
        aln_lines = []
        flag = False
        for l in out:
            if l.startswith('-- BEGIN'):
                m = re.search('\[([+\-\s\d]+)\|([+\-\s\d]+)\]', out[5])
                ref_frm, ref_s, z ,ref_e = m.group(1).strip().split()
                qry_frm, qry_s, z ,qry_e = m.group(2).strip().split()
                flag = True
            elif l.startswith('--   END'):
                flag = False
            if flag and re.match('^\d+', l):
                aln_lines.append(l)
        
        self.ref_frm, self.ref_s, self.ref_e = map(int, [ref_frm, ref_s, ref_e])
        self.qry_frm, self.qry_s, self.qry_e = map(int, [qry_frm, qry_s, qry_e])
        
        self.load_alignment_string(
            ''.join(aln_lines[i].split()[1] for i in range(0, len(aln_lines), 2)),
            ''.join(aln_lines[i].split()[1] for i in range(1, len(aln_lines), 2))
        )
        # Adjust for reference start
        # Nucmer uses 1-based numbering, we are going to use 0-based numbering
        for ap in self.aln_positions:
            if ap.ref_pos is not None:
                ap.ref_pos = ap.ref_pos + self.ref_s - 1

def align_promer(query_fa, ref_fa, outdir, debug=False):
    # Command 1: promer
    cmd1 = ['promer',
        '--prefix', os.path.join(outdir, 'promer'),
        '--extend',
        '--maxgap', '%d' % 70,
        '--minmatch', '%d' % 3,
        ref_fa,
        query_fa,
    ]
    
    # Command 2: delta-filter
    cmd2 = ['delta-filter',
        '-q',
        os.path.join(outdir, 'promer.delta'),
        '>',
        os.path.join(outdir, 'promer.filter')
    ]
    
    # Command 3: show-tiling
    cmd3 = ['show-tiling',
        '-a',
        '-i', '%.1f' % 0.6,
        '-l', '%d' % 200,
        '-v', '%.1f' % 60,
        os.path.join(outdir, 'promer.filter'),
        '>',
        os.path.join(outdir, 'promer.tiling'),
    ]
    command_runner([cmd1, cmd2, cmd3,], 'align_promer', debug)
    return os.path.join(outdir, 'promer.filter'), os.path.join(outdir, 'promer.tiling')

def align_nucmer(query_fa, ref_fa, outdir, debug=False):
    # Command 1: nucmer
    cmd1 = ['nucmer',
        '--prefix', os.path.join(outdir, 'nucmer'),
        '--extend',
        '--maxgap', '%d' % 200,
        '--minmatch', '%d' % 10,
        ref_fa,
        query_fa,
    ]
    
    # Command 2: delta-filter
    cmd2 = ['delta-filter',
        '-q',
        os.path.join(outdir, 'nucmer.delta'),
        '>',
        os.path.join(outdir, 'nucmer.filter')
    ]
    
    # Command 3: show-tiling
    cmd3 = ['show-tiling',
        '-a',
        '-i', '%.1f' % 0.6,
        '-l', '%d' % 200,
        '-v', '%.1f' % 60,
        os.path.join(outdir, 'nucmer.filter'),
        '>',
        os.path.join(outdir, 'nucmer.tiling'),
    ]
    command_runner([cmd1, cmd2, cmd3,], 'align_nucmer', debug)
    return os.path.join(outdir, 'nucmer.filter'), os.path.join(outdir, 'nucmer.tiling')

def show_aligns(ref, qry, delta):
    cmd1 = ['show-aligns', '-r', delta, ref, qry]
    return check_output(cmd1)

def assemble_to_ref(ref_fa, qry_fa, workdir):
    ref_dict = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    '''
    tmp_contigs_fa = os.path.join(workdir, 'query.fna')
    with open(tmp_contigs_fa, 'w') as outh:
        for i,ctg in enumerate(contigs):
            print >>outh, '>seq%02d\n%s' % (i, ctg)
    '''
    # Align contigs
    fil, til = align_nucmer(qry_fa, ref_fa, workdir)
    print >>sys.stderr, '\n'.join(l.strip() for l in open(til, 'rU'))    
    # Parse tiling
    tr_list = [TilingRow(l) for l in open(til, 'rU')]
    
    scaffolds = {}
    refs = sorted(ref_dict.keys())
    for ref in refs:
        ref_tr = [tr for tr in tr_list if tr.ref == ref]
        if not ref_tr:
            continue
        
        ranked = sorted(ref_tr, key=lambda x:x.pid, reverse=True)
        ranked.sort(key=lambda x:x.qry_alen, reverse=True)

        # Extract nucmer alignments
        nuc_alns = []
        for tr in ranked:
            out1 = show_aligns(tr.ref, tr.qry, fil)
            nuc_alns.append(NucmerAlignment(out1.strip('\n').split('\n')))
        
        scaffold = EmptyAlignment(str(ref_dict[ref].seq).lower())
        for na in reversed(nuc_alns):
            scaffold = scaffold.combine_alignments(na)
        
        scaffolds[ref] = scaffold
    
    return scaffolds


def parse_show_aligns(out):
    if type(out) is str:
        out = out.split('\n')
    
    # Get the alignment part
    aln = []
    flag = False
    for l in out:
        if l.startswith('-- BEGIN'):
            m = re.search('\[([+\-\s\d]+)\|([+\-\s\d]+)\]', out[5])
            ref_frm, ref_s, z ,ref_e = m.group(1).strip().split()
            qry_frm, qry_s, z ,qry_e = m.group(2).strip().split()
            flag = True
        elif l.startswith('--   END'):
            flag = False
        if flag and re.match('^\d+', l):
            aln.append(l)
    
    ref_frm, ref_s, ref_e = map(int, [ref_frm, ref_s, ref_e])
    qry_frm, qry_s, qry_e = map(int, [qry_frm, qry_s, qry_e])
    
    # Set query increment (based on frame/direction)
    qry_inc = 1 if qry_frm > 0 else -1
    
    # Init
    rpos = ref_s
    qpos = qry_s
    rpos_to_qseq = defaultdict(str)
    
    # Loop over pairs of alignment lines
    for i in range(0, len(aln), 2):  
        l1, l2 = aln[i], aln[i+1]
        new_rpos, rseq = l1.split()
        new_qpos, qseq = l2.split()
        assert int(new_rpos) == rpos
        assert int(new_qpos) == qpos
        for rb, qb in zip(rseq, qseq):
            if qb != '.':
                rpos_to_qseq[rpos] += qb
                qpos += qry_inc
            if rb != '.':
                rpos += 1
    return rpos_to_qseq


def pad_alignments(ref_dict, til_rows):
    padded = {}
    for ref in ref_dict.keys():
        ref_tr = sorted([tr for tr in til_rows if tr.ref == ref], key=lambda x:x.ref_s)
        if not ref_tr:
            continue 
        ref_len = ref_tr[0].ref_len
        padded[ref] = []
        for tr in ref_tr:
            n = '%s:%d-%d' % (tr.qry, tr.ref_s, tr.ref_e)
            pad = ['X'] + [''] * ref_len
            for p in range(1, ref_len+1):
                if p in tr.alignment:
                    pad[p] = tr.alignment[p].upper()
                else:
                    pad[p] = '.'
            if pad[0] == 'X': pad = pad[1:]
            padded[ref].append((n, pad))
    
    return padded


def print_padded(ref_dict, padded, out_fn='test.txt'):
    with open(out_fn, 'w') as outh:
        for ref in ref_dict.keys():
            if ref not in padded:
                continue
            print >>outh, '%s%s' % (ref.ljust(50), str(ref_dict[ref].seq).upper())
            for n, pad in padded[ref]:
                trunc = ''.join(pad)[1:len(ref_dict[ref])+1]
                print >>outh, '%s%s' % (n.ljust(50), trunc)

def overlapper(ref_dict, tr_list, quiet=True):
    scaffold = {}
    imputed = {}
    for ref in sorted(ref_dict.keys()):
        ref_tr = sorted([tr for tr in tr_list if tr.ref == ref], key=lambda x:x.ref_s)
        ref_len = len(ref_dict[ref])
        # Get regions of reference
        ref_breaks = set(chain.from_iterable([[tr.ref_s, tr.ref_e] for tr in ref_tr]))
        ref_breaks |= set([1, ref_len])
        ref_breaks = sorted(ref_breaks)
        ref_intervals = [(p, ref_breaks[i+1]) for i,p in enumerate(ref_breaks[:-1]) ]
        
        # Initialize
        scaffold[ref] = ['X'] + [''] * ref_len
        imputed[ref]  = ['X'] + [''] * ref_len
        
        for iv in ref_intervals:
            # Rank possible alignments by query length
            ranked = sorted([tr for tr in ref_tr if overlaps(iv, (tr.ref_s, tr.ref_e))],
                             key=lambda x:x.qry_alen, reverse=True)
            if len(ranked) == 0:
                if not quiet:
                    print >>sys.stderr, 'Reference %s:%d-%d is not covered' % (ref, iv[0], iv[1])
                for p in range(iv[0], iv[1]):
                    scaffold[ref][p] = 'n'
                    imputed[ref][p] = str(ref_dict[ref].seq[p-1]).lower()
            else:
                best = ranked[0]
                if not quiet:
                    if len(ranked) == 1:
                        print >>sys.stderr, 'Reference %s:%d-%d matched by %s:%d-%d' % (ref, iv[0], iv[1], best.qry, best.qry_s, best.qry_e)
                    else:
                        print >>sys.stderr, 'WARNING: Multiple matches found for Reference %s:%d-%d'
                        for tr in ranked:
                            print >>sys.stderr, '\tmatch: %s:%d-%d' % (tr.qry, tr.qry_s, tr.qry_e)
                for p in range(iv[0], iv[1]):
                    if p in best.alignment:
                        scaffold[ref][p] = best.alignment[p].upper()
                        imputed[ref][p]  = best.alignment[p].upper()
        # Remove the leading spacer
        if scaffold[ref][0] == 'X': scaffold[ref] = scaffold[ref][1:]
        if imputed[ref][0] == 'X': imputed[ref] = imputed[ref][1:]
    
    return scaffold, imputed

