#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse

from Bio import SeqIO

from utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from utils.sysutils import create_tempdir, remove_tempdir
from utils.sequtils import clean_seqnames, wrap
from utils.alignutils import assemble_to_ref

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--contigs_fa', type=existing_file, required=True,
                        help='Fasta file with assembled contigs')
    group1.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Fasta file with reference genome to scaffold against')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('Scaffold options')
    group2.add_argument('--seqname',
                        help='Name to append to scaffold sequence.')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Additional options')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=assemble_scaffold)


def assemble_scaffold(contigs_fa=None, ref_fa=None, outdir='.',
                      seqname='scaffold',
                      keep_tmp=False, debug=False):
    """ Scaffold contigs using reference sequence
    """
    # Check dependencies
    check_dependency('nucmer')
    check_dependency('delta-filter')
    check_dependency('show-tiling')
    
    # Outputs
    out_scaffold = os.path.join(outdir, 'scaffold.fa')
    out_imputed = os.path.join(outdir, 'imputed.fa')
    out_aln = os.path.join(outdir, 'aligned.fa')
    
    # Temporary directory
    tempdir = create_tempdir('assemble_scaffold')

    # Create fasta file with sequence IDs only (remove decription)
    tmp_contigs_fa = os.path.join(tempdir, 'query.fna')
    with open(tmp_contigs_fa, 'w') as outh:
        for n,s in clean_seqnames(open(contigs_fa, 'rU')):
            print >>outh, '>%s\t%s' % (n, wrap(s))
    
    with open(os.path.join(outdir, 'padded.out'), 'w') as pad_fh:
        scaffolds = assemble_to_ref(ref_fa, tmp_contigs_fa, tempdir, pad_fh=pad_fh)
        
    # Output scaffolds as FASTA
    with open(out_scaffold, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            s = scaffolds[ref].scaffold()
            print >>outh, '>%s\n%s' % (n, wrap(s))

    # Output imputed as FASTA
    with open(out_imputed, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            s = scaffolds[ref].imputed()
            print >>outh, '>%s\n%s' % (n, wrap(s))
    
    # Output alignments for other pipeline stages
    with open(out_aln, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)        
            print >>outh, '>REF|%s\n%s' % (n, scaffolds[ref].raln())
            print >>outh, '>%s\n%s' % (n, scaffolds[ref].qaln())
    
    if not keep_tmp:
        remove_tempdir(tempdir, 'assemble_scaffold')
    
    return out_scaffold

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scaffold contigs using reference sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
