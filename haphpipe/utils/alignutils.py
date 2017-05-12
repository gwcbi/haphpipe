#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
from subprocess import check_output
from collections import defaultdict

from Bio import SeqIO

from sysutils import command_runner
import alignobj

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

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

def align_nucmer(ref_fa, query_fa, outdir, debug=False):
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

def parse_show_aligns(out):
    """ Returns show_aligns """
    flag = False
    cur_report = None
    for l in out.strip('\n').split('\n'):
        if re.match('^--\s+BEGIN', l):
            cur_report = []
        if cur_report is not None:
            cur_report.append(l)
        if re.match('^--\s+END', l):
            yield alignobj.NucmerReferenceAlignment(cur_report)
            cur_report = None

def assemble_to_ref(ref_fa, qry_fa, workdir, pad_fh=None, debug=True):
    fil, til = align_nucmer(ref_fa, qry_fa, workdir)
    tr_byref = defaultdict(list)
    for l in open(til, 'rU'):
        tr = alignobj.TilingRow(l)
        tr_byref[tr.ref].append(tr)
    
    refs = sorted(tr_byref.keys())
    ref_dict = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    print >>sys.stderr, '\nReferences: %s' % ', '.join(refs)
    
    scaffolds = {}
    for ref in refs:
        if pad_fh is not None:
            empty = alignobj.EmptyReferenceAlignment(str(ref_dict[ref].seq).lower())
            print >>pad_fh, '%s%s' % (ref.ljust(40), empty.rseq().upper())
        scaffolds[ref] = alignobj.EmptyReferenceAlignment(str(ref_dict[ref].seq).lower())
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
                if debug:
                    print >>sys.stderr, "*" * 80
                    print >>sys.stderr, '\n'.join(aln_report)
                    print >>sys.stderr, "*" * 80
                nucaln = alignobj.NucmerReferenceAlignment(aln_report)
                print '%d-%d' % (nucaln.rstart, nucaln.rend)
                if pad_fh is not None:
                    pad = empty.merge_alignments(nucaln)
                    print >>pad_fh, '%s%s' % (tr.qry.ljust(40), pad.padded())
                scaffolds[ref] = scaffolds[ref].merge_alignments(nucaln)
    
    return scaffolds

