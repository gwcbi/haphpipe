#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import json
import argparse
from subprocess import check_output
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq

from ..utils.sysutils import existing_file, args_params
from ..utils.sequtils import wrap, parse_seq_id, region_to_tuple
from ..utils.gtfparse import gtf_parser, GTFRow
from ..utils.blastalign import load_slot_json

from ..utils.sysutils import PipelineStepError

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--align_json', type=existing_file, required=True,
                        help='''JSON file describing alignment (output of pairwise_align
                                stage)''')
    group1.add_argument('--outfile',
                        help='''Output file. Default is stdout''')
    group2 = parser.add_argument_group('Extract pairwise options')                                
    group2.add_argument('--outfmt', type=str, default='nuc_fa',
                        choices=['nuc_fa', 'aln_fa', 'amp_gtf', 'tsv',
                                  'prot_fa', ],
                        help='Format for output')
    group2.add_argument('--refreg',
                        help='''Reference region''')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=extract_pairwise)


def extract_pairwise(align_json=None, outfile=None,
        outfmt=None, refreg=None,
        debug=False,
    ):
    outh = sys.stdout if outfile is None else open(outfile, 'w')
    
    if outfmt == 'nuc_fa' or outfmt == 'prot_fa':
        jaln = load_slot_json(align_json, 'padded_alignments')
        if refreg is None:
            for newname, alignment in jaln.iteritems():
                nucstr = ''.join(t[2] for t in alignment if t[3] != -1)
                nucstr = nucstr.replace('*', 'N')
                print >>outh, '>%s' % newname                
                if outfmt == 'nuc_fa':
                    print >>outh, wrap(nucstr)
                else:
                    s = Seq(nucstr[:(len(nucstr)/3)*3])
                    print >>outh, wrap(str(s.translate()))
        else:
            refmap = {parse_seq_id(k)['ref']:k for k in jaln.keys()}
            chrom, ref_s, ref_e = region_to_tuple(refreg)
            ref_s = ref_s - 1
            alignment = jaln[refmap[chrom]]
            
            # Get alignment start
            for aln_s in xrange(len(alignment)):
                if alignment[aln_s][0] == ref_s:
                    break
                while alignment[aln_s][3] == -1:
                    aln_s += 1
            
            # Get alignment end
            for aln_e in xrange(len(alignment)-1, -1, -1):
                if alignment[aln_e][0] == ref_e:
                    break
            while alignment[aln_e][3] == -1:
                aln_e += -1

            nucstr = ''.join(t[2] for t in alignment[aln_s:aln_e] if t[3] != -1)
            nucstr = nucstr.replace('*', 'N')
            print >>outh, '>%s (%s)' % (refmap[chrom], refreg)
            if outfmt == 'nuc_fa':            
                print >>outh, wrap(nucstr)
            else:
                s = Seq(nucstr[:(len(nucstr)/3)*3])
                print >>outh, wrap(str(s.translate()))                

    elif outfmt == 'aln_fa':
        jaln = load_slot_json(align_json, 'padded_alignments')
        for newname, alignment in jaln.iteritems():
            aid = parse_seq_id(newname)
            rstr = ''.join(t[1] for t in alignment).replace('*', 'N')
            qstr = ''.join(t[2] for t in alignment).replace('*', 'N')
            print >>outh, '>ref|%s' % aid['ref']
            print >>outh, wrap(rstr)
            print >>outh, '>sid|%s' % aid['sid']
            print >>outh, wrap(qstr)

    elif outfmt == 'amp_gtf':
        jgtf = load_slot_json(align_json, 'padded_gtf')
        print >>outh, '\n'.join(_ for _ in jgtf)

    elif outfmt == 'tsv':
        jaln = load_slot_json(align_json, 'padded_alignments')
        for newname, alignment in jaln.iteritems():
            print >>outh, '# %s' % newname
            for l in alignment:
                print >>outh, '\t'.join(str(_) for _ in l)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Extract from consensus based on pairwise alignment')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
