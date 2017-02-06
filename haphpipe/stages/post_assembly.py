#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse

from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.sequtils import fastagen
from ..utils.alignobj import ReferenceAlignment
from ..utils.helpers import pairwise

from ..utils.sysutils import PipelineStepError

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--aln_fa', type=existing_file,
                        help='Alignment file')
    parser.set_defaults(func=post_assembly)

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

"""
def samtools_depth(self, ):
    check_dependency('samtools')
    cmd1 = ['samtools', 'depth', bamfile,]
    Popen(cmd1, shell=True, stdout=PIPE, stderr=PIPE)
    o,e = p1.communicate()
"""

def post_assembly(aln_fa=None):
    ra = load_ref_align(aln_fa)
    print ra


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix consensus sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
