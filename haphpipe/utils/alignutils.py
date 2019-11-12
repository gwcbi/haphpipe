# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import str
from builtins import range
from builtins import object
import sys
import os
import re
from subprocess import check_output
from collections import defaultdict

from Bio import SeqIO

from haphpipe.utils import sysutils


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


"""Classes
"""


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
        assert len(raln) == len(
            qaln), 'Error: Alignments are different lengths'
        self.aln_positions = []

        # Convert all gapchars to '.'
        for rb, qb in zip(raln, qaln):
            rb = '.' if rb in self.gapchar else rb
            qb = '.' if qb in self.gapchar else qb
            self.aln_positions.append((rb, qb))

        self.alen = len(self.aln_positions)
        self._index_alignment()

    def _index_alignment(self):
        # print >>sys.stderr, "Indexing alignment"
        if self.alen != len(self.aln_positions):
            # print >>sys.stderr, 'Warning: self.alen != len(self.aln_positions)'
            self.alen = len(self.aln_positions)

        # Initialize
        self._rpos_to_apos = {}  # one-to-one with some apos not included
        self._qpos_to_apos = {}  # one-to-one with some apos not included
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
            print("WARNING: position %d is outside reference boundaries" % rp,
                  file=sys.stderr)
            return self.qstart
        elif rp >= self.rend:
            print("WARNING: position %d is outside reference boundaries" % rp,
                  file=sys.stderr)
            return self.qend
        else:
            ap = self._rpos_to_apos[rp]
            while 0 <= ap <= self.alen:
                if ap in self._apos_to_qpos:
                    return self._apos_to_qpos[ap]
                ap += incr
            print("WARNING: position %d is outside alignment" % ap,
                  file=sys.stderr)
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

    def scaffold2(self, collapse=True):
        for aln_s in range(len(self.aln_positions)):
            if self.aln_positions[aln_s][1] != '?':
                break
        for aln_e in range(len(self.aln_positions) - 1, 0, -1):
            if self.aln_positions[aln_e][1] != '?':
                break

        ref_s = self._apos_to_rpos[aln_s]
        ref_e = self._apos_to_rpos[aln_e]

        # Extract the scaffold string
        ret = self.qseq().upper()
        ret = ret.strip('?')
        if collapse:
            ret = re.sub('[\?\.\-]', '', ret)
        return ret, ref_s, ref_e


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
        # Get beginning and end of report
        for i,rl in enumerate(aln_report):
            if re.match('^--\s+BEGIN', rl):
                l_begin = i
            if re.match('^--\s+END', rl):
                l_end = i
        aln_report = aln_report[l_begin:(l_end+1)]
        # print('\n'.join(aln_report))
        assert re.match('^--\s+BEGIN', aln_report[
            0]), 'First line of alignment is not "BEGIN"'
        assert re.match('^--\s+END',
                        aln_report[-1]), 'Last line of alignment is not "END"'
        m = re.search('\[([+\-\s\d]+)\|([+\-\s\d]+)\]', aln_report[0])
        ref_frm, ref_s, z, ref_e = m.group(1).strip().split()
        qry_frm, qry_s, z, qry_e = m.group(2).strip().split()
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
        self.ref_frm, self.ref_s, self.ref_e = list(map(int,
                                                   [ref_frm, ref_s, ref_e]))
        self.qry_frm, self.qry_s, self.qry_e = list(map(int,
                                                   [qry_frm, qry_s, qry_e]))
        aln_lines = [l for l in aln_report if re.match('^\d+', l)]
        self.load_alignment(
            ''.join(
                aln_lines[i].split()[1] for i in range(0, len(aln_lines), 2)),
            ''.join(
                aln_lines[i].split()[1] for i in range(1, len(aln_lines), 2))
        )
        # Adjust for reference start
        # Nucmer uses 1-based numbering, we are going to use 0-based numbering
        self.adjust_ref_start(self.ref_s - 1)


class TilingRow(object):
    """ Alignment specification from show-tiling (mummer)
    """

    def __init__(self, l, seqdict=None):
        fields = l.strip('\n').split('\t')
        self.ref_s, self.ref_e, self.qry_s, self.qry_e = list(map(int, fields[:4]))
        self.ref_alen, self.qry_alen = list(map(int, fields[4:6]))
        self.pid = float(fields[6])
        self.ref_len, self.qry_len = list(map(int, fields[7:9]))
        self.ref_cov, self.qry_cov = list(map(float, fields[9:11]))
        self.ref, self.qry = fields[11:]
        if seqdict is not None:
            self.add_seq(seqdict)

    def __str__(self):
        return '\t'.join(
            map(str, [self.ref_s, self.ref_e, self.qry_s, self.qry_e,
                      self.ref_alen, self.qry_alen, self.pid,
                      self.ref_len, self.qry_len,
                      self.ref_cov, self.qry_cov,
                      self.ref, self.qry]))


## Utility functions

def align_promer(
        query_fa, ref_fa, outdir, quiet=False, logfile=None, debug=False
    ):
    """

    Args:
        query_fa:
        ref_fa:
        outdir:
        quiet:
        logfile:
        debug:

    Returns:

    """
    # Outputs
    out_del = os.path.join(outdir, 'promer.delta')
    out_fil = os.path.join(outdir, 'promer.filter')
    out_til = os.path.join(outdir, 'promer.tiling')

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
    cmd2 = [
        'delta-filter', '-q', out_del, '>', out_fil
    ]
    
    # Command 3: show-tiling
    cmd3 = ['show-tiling',
        '-a',
        '-i', '%.1f' % 0.6,
        '-l', '%d' % 200,
        '-v', '%.1f' % 60,
        out_fil,
        '>',
        out_til,
    ]
    sysutils.command_runner(
        [cmd1, cmd2, cmd3,], 'align_promer', quiet, logfile, debug
    )
    return out_fil, out_til


def align_nucmer(
        query_fa, ref_fa, outdir, min_contig_len=200,
        quiet=False, logfile=None, debug=False
    ):
    """

    Args:
        query_fa:
        ref_fa:
        outdir:
        quiet:
        logfile:
        debug:

    Returns:

    """
    # Outputs
    out_del = os.path.join(outdir, 'nucmer.delta')
    out_fil = os.path.join(outdir, 'nucmer.filter')
    out_til = os.path.join(outdir, 'nucmer.tiling')

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
    cmd2 = [
        'delta-filter', '-q', out_del, '>', out_fil
    ]
    
    # Command 3: show-tiling
    cmd3 = ['show-tiling',
        '-a',
        '-i', '%.1f' % 0.6,
        '-l', '%d' % min_contig_len,
        '-v', '%.1f' % 60,
        out_fil,
        '>',
        out_til,
    ]
    sysutils.command_runner(
        [cmd1, cmd2, cmd3, ], 'align_nucmer', quiet, logfile, debug
    )
    return out_fil, out_til


def show_aligns(ref, qry, delta, quiet=False, logfile=None, debug=False):
    cmd1 = ['show-aligns', '-r', delta, ref, qry]
    return check_output(cmd1)


def parse_show_aligns(out):
    """ Returns show_aligns """
    flag = False
    cur_report = None
    outlines = out.decode('utf-8').strip('\n').split('\n')
    for l in outlines:
        if re.match('^--\s+BEGIN', l):
            cur_report = []
        if cur_report is not None:
            cur_report.append(l)
        if re.match('^--\s+END', l):
            yield NucmerReferenceAlignment(cur_report)
            cur_report = None


def assemble_to_ref(
        qry_fa, ref_fa, outdir, pad_fh=None,
        quiet=False, logfile=None, debug=False
    ):
    """

    Args:
        qry_fa:
        ref_fa:
        outdir:
        pad_fh:
        quiet:
        logfile:
        debug:

    Returns:

    """
    # Align query to reference
    fil, til = align_nucmer(
        qry_fa, ref_fa, outdir,
        quiet=quiet, logfile=logfile, debug=debug
    )
    if debug:
        return None

    # Parse tiling rows
    tr_byref = defaultdict(list)
    for l in open(til, 'rU'):
        tr = TilingRow(l)
        tr_byref[tr.ref].append(tr)

    # Load reference(s)
    refs = sorted(tr_byref.keys())
    ref_dict = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    sysutils.log_message(
        '\nReferences: %s\n' % ', '.join(refs), quiet, logfile
    )
    
    scaffolds = {}
    for ref in refs:
        if pad_fh is not None:
            empty = EmptyReferenceAlignment(str(ref_dict[ref].seq).lower())
            print('%s%s' % (ref.ljust(40), empty.rseq().upper()), file=pad_fh)
        scaffolds[ref] = EmptyReferenceAlignment(str(ref_dict[ref].seq).lower())
        # Rank hits so that worst hit is in index 0 (best at the end)
        ranked = sorted(tr_byref[ref], key=lambda x:x.pid)
        ranked.sort(key=lambda x:x.qry_alen)
        
        for tr in ranked:
            out = show_aligns(tr.ref, tr.qry, fil)
            # May be multiple alignments
            flag = False
            aln_reports = []
            for l in out.strip('\n').split('\n'):
                if re.match('^--\s+BEGIN', l):
                    aln_reports.append(list())
                    flag = True
                if flag:
                    aln_reports[-1].append(l)
                if re.match('^--\s+END', l):
                    flag = False
            
            for aln_report in aln_reports:
                # if debug:
                #     print("*" * 80, file=sys.stderr)
                #     print('\n'.join(aln_report), file=sys.stderr)
                #     print("*" * 80, file=sys.stderr)
                nucaln = NucmerReferenceAlignment(aln_report)
                # print('%d-%d' % (nucaln.rstart, nucaln.rend))
                if pad_fh is not None:
                    pad = empty.merge_alignments(nucaln)
                    print('%s%s' % (tr.qry.ljust(40), pad.padded()), file=pad_fh)
                scaffolds[ref] = scaffolds[ref].merge_alignments(nucaln)
    
    return scaffolds

