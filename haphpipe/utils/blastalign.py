# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import absolute_import
import sys
import os
import json
from subprocess import Popen, PIPE
from collections import defaultdict, Counter

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

from .helpers import overlaps

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"

class BlastxAlignment(object):
    def __init__(self, ref, query, cdspos, workdir="."):
        self.ref = ref    
        self.query = query
        self.cds_s, self.cds_e = cdspos
        self.aa_align = []
        self.nuc_align = []
        
        jobj = self.run_blast(showaln=True, workdir=workdir)
        self.parse_alignment(jobj)
    
    def run_blast(self, showaln=True, workdir="."):
        sfile = os.path.join(workdir, 'tmp.faa')
        qfile = os.path.join(workdir, 'tmp.fna')
        
        with open(sfile, 'w') as outh:
            print('>subject', file=outh)
            print(str(self.ref[self.cds_s:self.cds_e].translate()), file=outh)
        
        with open(qfile, 'w') as outh:
            print('>query', file=outh)
            print(str(self.query), file=outh)
        
        if showaln:
            p0 = Popen(['blastx',
                        '-query', qfile,
                        '-subject', sfile,
                       ], stdout=PIPE, stderr=PIPE)
            o,e = p0.communicate()
            print(o, file=sys.stdout)
        
        p = Popen(['blastx',
                   '-query', qfile,
                   '-subject', sfile,
                   '-outfmt', '12',
                   ], stdout=PIPE, stderr=PIPE)
        o, e = p.communicate()
        return json.loads(o)
    
    def parse_alignment(self, jobj, alnidx=0):
        aln = jobj['Seq_annot']['data']['align'][alnidx]
        for i,seg in enumerate(aln['segs']['std']):
            nzip,azip = self.parse_seg(seg)
            self.nuc_align.extend(nzip)
            self.aa_align.extend(azip)
    
    def parse_seg(self, seg_json):
        assert len(seg_json['loc']) == 2
        qjson, sjson = seg_json['loc']
        
        # Process query
        if 'int' in qjson:
            assert qjson['int']['id']['local']['str'].startswith('Query')
            # Nucleotide
            qrynuc_s = qjson['int']['from']
            qrynuc_e = qjson['int']['to'] + 1
            qrynuc_seq = self.query[qrynuc_s:qrynuc_e]
            qrynuc_pos = xrange(qrynuc_s, qrynuc_e)
            # Protein
            qryaa_seq = qrynuc_seq.translate()
            qryaa_pos = xrange(qrynuc_s, qrynuc_e, 3)
        else:
            assert 'empty' in qjson
            qrynuc_seq = None
        
        if 'int' in sjson:
            assert sjson['int']['id']['local']['str'].startswith('Subject')
            # Protein
            refaa_s = sjson['int']['from']
            refaa_e = sjson['int']['to'] + 1
            refaa_seq = self.ref[self.cds_s:self.cds_e].translate()[refaa_s:refaa_e]
            refaa_pos = xrange(refaa_s, refaa_e)
            # Nucleotide
            refnuc_s = self.cds_s + (refaa_s * 3)
            refnuc_e = self.cds_s + (refaa_e * 3)
            refnuc_seq = self.ref[refnuc_s:refnuc_e]
            refnuc_pos = xrange(refnuc_s, refnuc_e)
        else:
            assert 'empty' in sjson
            refnuc_seq = None
        
        if qrynuc_seq is None:
            qrynuc_seq = Seq('-' * (refnuc_e - refnuc_s))
            qrynuc_pos = [-1] * (refnuc_e - refnuc_s)
            qryaa_seq = Seq('-' * (refaa_e - refaa_s))
            qryaa_pos = [-1] * (refaa_e - refaa_s)
        
        if refnuc_seq is None:
            refnuc_seq = Seq('-' * (qrynuc_e - qrynuc_s))
            refnuc_pos = [-1] * (qrynuc_e - qrynuc_s)
            refaa_seq = Seq('-' * len(qryaa_seq))
            refaa_pos = [-1] * len(qryaa_seq)
        
        return zip(refnuc_pos, refnuc_seq, qrynuc_seq, qrynuc_pos), zip(refaa_pos, refaa_seq, qryaa_seq, qryaa_pos)

def discontinuous_query(aln):
    ''' Find locations in alignment where query is not continuous '''
    ret = []
    prev = aln[0][3] - 1
    for i, t in enumerate(aln):
        if t[3] == -1: continue
        if t[3] != prev + 1:
            ret.append((i, prev, t[3]))
        prev = t[3]
    return ret

def aligned_seqs_to_list(s1, s2, pos1, pos2):
    ''' Convert aligned sequence strings to list
            s1 - aligned sequence 1
            s2 - aligned sequence 2
            pos1 - start position of sequence 1
            pos2 - start position of sequence 2
    '''
    ret = []
    cur1, cur2 = pos1, pos2
    for b1, b2 in zip(s1, s2):
        p1 = p2 = None
        if b1 == '-':
            p1 = -1
        else:
            p1 = cur1
            cur1 += 1
        if b2 == '-':
            p2 = -1
        else:
            p2 = cur2
            cur2 += 1
        ret.append((p1, b1, b2, p2))
    return ret

def nuc_align_insert(merged, mr, refseq, qryseq):
    ''' Nucleotide alignment for unaligned positions within merged alignment '''
    ref_s, ref_e = merged[mr-1][0]+1, merged[mr][0]-1
    qry_s, qry_e = merged[mr-1][3]+1, merged[mr][3]-1
    alns = pairwise2.align.globalms(refseq[ref_s:ref_e+1], qryseq[qry_s:qry_e+1], 2, -1, -10, -1)
    best = alns[0]
    print(pairwise2.format_alignment(*best), file=sys.stdout)
    insert = aligned_seqs_to_list(best[0], best[1], ref_s, qry_s)
    return insert


def nuc_align_head(merged, refseq, qryseq):
    ''' Nucleotide alignment for unaligned positions to the left of merged alignment '''
    qry_s = 0
    qry_e = merged[0][3]
    if qry_s >= qry_e:
        return []
    ref_e = merged[0][0]
    ref_s = ref_e - (qry_e - qry_s)
    alns = pairwise2.align.globalms(refseq[ref_s:ref_e], qryseq[qry_s:qry_e], 2, -1, -10, -1)
    best = alns[0]
    print(pairwise2.format_alignment(*best), file=sys.stdout)
    head = aligned_seqs_to_list(best[0], best[1], ref_s, qry_s)
    # Trim off query gaps at the end
    while head[0][2] == '-':
        head = head[1:]
    return head


def nuc_align_tail(merged, refseq, qryseq):
    ''' Nucleotide alignment for unaligned positions to the right of merged alignment '''
    qry_s = merged[-1][3] + 1
    qry_e = len(qryseq)
    if qry_s >= qry_e:
        return []
    ref_s = merged[-1][0] + 1
    ref_e = ref_s + (qry_e - qry_s)
    alns = pairwise2.align.globalms(refseq[ref_s:ref_e+1], qryseq[qry_s:qry_e+1], 2, -1, -10, -1)
    best = alns[0]
    print(pairwise2.format_alignment(*best), file=sys.stdout)
    tail = aligned_seqs_to_list(best[0], best[1], ref_s, qry_s)
    # Trim off query gaps at the end
    while tail[-1][2] == '-':
        _ = tail.pop()
    return tail

def validate_alignment(aln, refseq, qryseq):
    assert aln[0][3] == 0, 'First query position is not 0'
    prev_r = aln[0][0] - 1
    prev_q = aln[0][3] - 1
    for tup in aln:
        if tup[0] == -1:
            assert tup[1] == '-', '%s: reference should be gap ("-")' % tup
        elif tup[0] == -2:
            assert tup[1] == '*', '%s: reference should be pad ("*")' % tup
        else:
            assert tup[0] == prev_r + 1, '%s: reference pos should be %d' % (tup, (prev_r + 1))
            assert refseq[tup[0]] == tup[1], '%s: reference base should be %s' % (tup, refseq[tup[0]])
            prev_r = tup[0]
        
        if tup[3] == -1:
            assert tup[2] == '-', '%s: query should be gap ("-")' % tup
        elif tup[3] == -2:
            assert tup[2] == '*', '%s: query should be pad ("*")' % tup
        else:
            assert tup[3] == prev_q + 1, '%s: query pos should be %d' % (tup, (prev_q + 1))
            assert qryseq[tup[3]] == tup[2], '%s: query base should be %s' % (tup,  qryseq[tup[3]])
            prev_q = tup[3]
    return True

def get_seg_stats(alnseg):
    statkeys = ['match', 'mismatch', 'rgap', 'qgap', 'uncalled',]
    c = Counter()
    for tup in alnseg:
        if tup[0] == -1: c['rgap'] += 1
        elif tup[3] == -1: c['qgap'] += 1
        elif tup[1] == tup[2]: c['match'] += 1
        elif tup[2] == '*': c['uncalled'] += 1
        elif tup[1] != tup[2]: c['mismatch'] += 1
    return {k:c[k] for k in statkeys}

def called_regions(alnseg):
    called = []
    flag = False
    for tup in alnseg:
        if flag:
            if tup[2] == '*':
                called[-1]['end'] = tup
                flag = False
        else:
            if tup[2] != '*':
                called.append({'start': tup, 'end': None})
                flag = True
    if flag:
        assert called[-1]['end'] is None
        called[-1]['end'] = tup
    
    return [(d['start'][3]+1, d['end'][3]) for d in called]


def load_slot_json(infile, slotname='padded_alignments'):
    fh = open(infile, 'rU') if type(infile) is str else infile
    jobj = json.loads(fh.read())
    return jobj[slotname]


def alignAA(refrec, qryrec, cds, altcds=list(), workdir="."):
    ''' Perform blastx alignment
    
        Use blastx to align translated nucleotide query to protein reference. The user
        provides nucleotide sequence for the reference and query as well as coordinates
        for open reading frames in the reference. The query is first aligned to the
        primary ORF given by the cds parameter. The query is then aligned to other
        alternate ORFs if provided. The alignments are merge and the remaining sections
        of the query that did not align to ORF are nucleotide aligned.
    '''
    ref = refrec.seq    
    qry = qryrec.seq
    bxa = BlastxAlignment(ref, qry, cds, workdir)
    
    merged = bxa.nuc_align
    altalns = [BlastxAlignment(ref, qry, acds, workdir) for acds in altcds]
    # print '\nstart: %s\nend: %s\n' % ('%d %s %s %d' % merged[0], '%d %s %s %d' % merged[-1])
    for bx in altalns:
        # print '\nstart: %s\nend: %s\n' % ('%d %s %s %d' % bx.nuc_align[0], '%d %s %s %d' % bx.nuc_align[-1])
        mrange = (merged[0][3], merged[-1][3])               # merged range
        orange = (bx.nuc_align[0][3], bx.nuc_align[-1][3])   # other range
        if not overlaps(mrange, orange):
            if orange[1] < mrange[0]:                        # other comes before merged
                merged = bx.nuc_align + merged
            else:                                            # merged comes before other
                assert mrange[1] < orange[0]
                merged = merged + bx.nuc_align
        else:
            for gap_idx, gap_s, gap_e in discontinuous_query(merged):
                if overlaps(orange, (gap_s, gap_e)):
                    for sidx in xrange(len(bx.nuc_align)):
                        if gap_s < bx.nuc_align[sidx][3] < gap_e:
                            break
                    for eidx in xrange(len(bx.nuc_align)-1, -1, -1):
                        if gap_s < bx.nuc_align[eidx][3] < gap_e:
                            break
                    assert merged[gap_idx-1][3] < bx.nuc_align[sidx][3]
                    assert merged[gap_idx][3] > bx.nuc_align[eidx][3]
                    merged = merged[:gap_idx] + bx.nuc_align[sidx:eidx+1] + merged[gap_idx:]
                    break
    
    # Nucleotide alignment for unaligned positions within merged alignment
    internal_gaps = discontinuous_query(merged)
    while internal_gaps:
        gap_idx, gap_s, gap_e = internal_gaps[0]
        # print '\ngap start: %s\ngap end: %s\n' % ('%d %s %s %d' % merged[gap_idx-1], '%d %s %s %d' % merged[gap_idx])
        ins = nuc_align_insert(merged, gap_idx, ref, qry)
        merged = merged[:gap_idx] + ins + merged[gap_idx:]
        internal_gaps = discontinuous_query(merged)

    # Nucleotide alignment for unaligned positions to the left of merged alignment
    # print '\ngap start: %s\ngap end: %s\n' % ('%d %s %s %d' % merged[0], '%d %s %s %d' % merged[0])
    h = nuc_align_head(merged, ref, qry)
    merged = h + merged
    
    # Nucleotide alignment for unaligned positions to the right of merged alignment
    # print '\ngap start: %s\ngap end: %s\n' % ('%d %s %s %d' % merged[-1], '%d %s %s %d' % merged[-1])
    t = nuc_align_tail(merged, ref, qry)
    merged = merged + t
    
    if validate_alignment(merged, ref, qry):
        print("Alignment validation passed", file=sys.stderr)
    
    return bxa, merged
