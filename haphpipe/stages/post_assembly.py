#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from subprocess import check_output
from collections import defaultdict

from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.sequtils import fastagen, unambig_intervals
from ..utils.alignobj import ReferenceAlignment
from ..utils.helpers import merge_interval_list, percentile # pairwise

from ..utils.sysutils import PipelineStepError

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--alignment', type=existing_file, required=True,
                        help='Alignment file (BAM)')
    group1.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Reference sequence used to align reads (fasta)')                        
    group1.add_argument('--depth_txt', type=existing_file,
                        help='File with depth per position (samtools format)')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    group2 = parser.add_argument_group('Post assembly options')
    group2.add_argument('--min_depth', type=int,
                        help='Minimum depth to consider position covered')
    group2.add_argument('--max_ambig', type=int,
                        help='''Maximum size of ambiguous sequence within a reconstruction
                                region''')                    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=post_assembly)

def get_intervals(covlines, ref_fa, max_ambig=None, min_depth=None, debug=False):
    # Load reference fasta
    seqs = {}
    with open(ref_fa, 'rU') as fh:
        for n, s in fastagen(fh):
            seqs[n.split()[0]] = s    
    
    # Load coverage data
    covdict = defaultdict(lambda: defaultdict(int))
    for l in covlines:
        covdict[l[0]][int(l[1]) - 1] = int(l[2])
    
    # Set of shared reference
    refs = set(seqs.keys()) & set(covdict.keys())
    
    refined = {}
    for ref in refs:
        refined[ref] = []
        if min_depth is None:
            allcovs = sorted(covdict[ref].values())
            min_depth = max(20, percentile(allcovs, 0.5) * 1e-2)
            print >>sys.stderr, 'Using cutoff %s for %s' % (min_depth, ref)
        
        for uiv in unambig_intervals(seqs[ref], max_ambig):
            if debug:
                print >>sys.stderr, 'Unambiguous interval: %s:%d-%d' % (ref, uiv[0], uiv[1])
            start = uiv[0]
            while covdict[ref][start] < min_depth and start < uiv[1]:
                start += 1
            if start == uiv[1]:
                continue
            end = uiv[1] - 1
            while covdict[ref][end] < min_depth and end > uiv[0]:
                end -= 1
            if debug:
                print >>sys.stderr, 'Covered interval: %s:%d-%d' % (ref, start, end)
            refined[ref].append((start,end))
    
    return refined          

def get_covered_intervals(covlines, cutoff=None):# , indel=10, debug=False):
    """ Identify covered intervals """
    covdict = defaultdict(dict)
    for l in covlines:
        covdict[l[0]][int(l[1])] = int(l[2])
    
    ret = {}
    for ref in covdict.keys():
        # Calculate a cutoff value
        if cutoff is None:
            allcovs = sorted(covdict[ref].values())
            cutoff = max(20, percentile(allcovs, 0.5) * 1e-2)
            print >>sys.stderr, 'Using cutoff %s for %s' % (cutoff, ref)
        
        # Get intervals that are strictly greater than cutoff
        ivs = []
        for pos in sorted(covdict[ref].keys()):
            if covdict[ref][pos] > cutoff:
                if ivs and ivs[-1][1] + 1 == pos:
                    ivs[-1][1] = pos
                else:
                    ivs.append([pos, pos])
        
        # Allow for small indels (indel parameter)
        ivs = merge_interval_list(ivs, indel)
        ret[ref] = [(iv[0],iv[1]) for iv in ivs]
        """
        if debug:
            print >>sys.stderr, 'Initial intervals:'
            for iv in ivs:
                print >>sys.stderr, '%s:%d-%d' % (ref, iv[0], iv[1])
        
        # Trim beginning and end according to cutoff
        ret[ref] = []
        for iv in ivs:
            news = newe = None
            for pos in range(iv[0], iv[1]):
                if pos in covdict[ref] and covdict[ref][pos] > cutoff:
                    news = pos
                    break
            for pos in range(iv[1], iv[0], -1):
                if pos in covdict[ref] and covdict[ref][pos] > cutoff:
                    newe = pos
                    break
            if news is not None and newe is not None:
                ret[ref].append((news,newe))
        """
    return ret

'''
def load_ref_align(aln_fa, refname=None):
    seqiter = pairwise(fastagen(open(aln_fa, 'rU')))
    if refname is None:
        s1, s2 = seqiter.next()
        return ReferenceAlignment(rname=s1[0], raln=s1[1], 
                                  qname=s2[0], qaln=s2[1])
    else:
        for s1, s2 in seqiter:
            if s1[0] == refname:
                return ReferenceAlignment(rname=s1[0], raln=s1[1], 
                                          qname=s2[0], qaln=s2[1])
    raise PipelineStepError("Reference %s not found")
'''

def samtools_depth(bamfile):
    check_dependency('samtools')
    cmd1 = ['samtools', 'depth', bamfile, ]
    return check_output(cmd1)


def post_assembly(alignment=None, ref_fa=None, depth_txt=None, outdir='.',
        min_depth=None, max_ambig=200, 
        keep_tmp=False, debug=False,
    ):
    if depth_txt is None:
        out = samtools_depth(alignment)
        covlines = [l.split('\t') for l in out.strip('\n').split('\n')]
    else:
        covlines = [l.strip('\n').split('\t') for l in open(depth_txt, 'rU')]
    
    # Create covered intervals output
    # cov_ivs = get_covered_intervals(covlines, min_depth)
    cov_ivs = get_intervals(covlines, ref_fa, max_ambig, min_depth, debug)
    print >>sys.stderr, 'Covered intervals:'
    with open(os.path.join(outdir, 'covered.intervals'), 'w') as outh:
        for ref, ivs in cov_ivs.iteritems():
            for iv in ivs:
                print >>sys.stderr, '%s:%d-%d' % (ref, iv[0], iv[1])
                print >>outh, '%s:%d-%d' % (ref, iv[0], iv[1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix consensus sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
