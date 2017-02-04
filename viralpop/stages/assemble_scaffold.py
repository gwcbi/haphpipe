#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse
from collections import defaultdict

from Bio import SeqIO

from utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from utils.sysutils import create_tempdir, remove_tempdir
from utils.sequtils import clean_seqnames, wrap
from utils.alignutils import TilingRow, NucmerAlignment, EmptyAlignment
from utils.alignutils import align_nucmer, show_aligns
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


def new_assemble_scaffold(contigs_fa=None, ref_fa=None, outdir='.',
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
    out_aln = os.path.join(outdir, 'scaffold.tsv')
    
    # Temporary directory
    tempdir = create_tempdir('assemble_scaffold')

    # Create fasta file with sequence IDs only (remove decription)
    tmp_contigs_fa = os.path.join(tempdir, 'query.fna')
    with open(tmp_contigs_fa, 'w') as outh:
        for n,s in clean_seqnames(open(contigs_fa, 'rU')):
            print >>outh, '>%s\t%s' % (n, wrap(s))
    
    scaffolds = assemble_to_ref(ref_fa, tmp_contigs_fa, tempdir)
    
    # Output scaffolds as FASTA
    with open(out_scaffold, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            s = scaffolds[ref].get_scaffold()
            print >>outh, '>%s\n%s' % (n, wrap(s))
    
    # Output alignments for other pipeline stages
    with open(out_aln, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            print >>outh, '>%s' % n
            print >>outh, '\n'.join('\t'.join(row) for row in scaffolds[ref].as_matrix())
    
    if not keep_tmp:
        remove_tempdir(tempdir, 'assemble_scaffold')
    
    return out_scaffold

def old_assemble_scaffold(contigs_fa=None, ref_fa=None, outdir='.',
                      seqname='scaffold',
                      keep_tmp=False, debug=False):
    """ Scaffold contigs using reference sequence
    """
    # Check dependencies
    check_dependency('nucmer')
    check_dependency('delta-filter')
    check_dependency('show-tiling')
    
    # Outputs
    out_scaffold = os.path.join(outdir, 'oscaffold.fa')
    out_aln = os.path.join(outdir, 'oscaffold.tsv')
    
    # Temporary directory
    tempdir = create_tempdir('assemble_scaffold')
    
    # Create fasta file with sequence IDs only (remove decription)
    tmp_contigs_fa = os.path.join(tempdir, 'query.fna')
    with open(tmp_contigs_fa, 'w') as outh:
        for n,s in clean_seqnames(open(contigs_fa, 'rU')):
            print >>outh, '>%s\t%s' % (n, wrap(s))
    
    # Align contigs
    fil, til = align_nucmer(tmp_contigs_fa, ref_fa, tempdir)
    print >>sys.stderr, '\n'.join(l.strip() for l in open(til, 'rU'))
    
    # Parse tiling
    tr_list = [TilingRow(l) for l in open(til, 'rU')]
    
    # Only references that have alignments
    ref_dict = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    refs = list(set([tr.ref for tr in tr_list]))    
    scaffolds = {}
    
    for ref in refs:
        ranked = [tr for tr in tr_list if tr.ref == ref]
        # Rank by percent ID
        ranked.sort(key=lambda x:x.pid, reverse=True)
        # Then rank by aligned contig length
        ranked.sort(key=lambda x:x.qry_alen, reverse=True)
    
        # Extract nucmer alignments
        nuc_alns = []
        for tr in ranked:
            out1 = show_aligns(tr.ref, tr.qry, fil)
            nuc_alns.append(NucmerAlignment(out1.strip('\n').split('\n')))
        
        for na, tr in zip(nuc_alns, ranked):
            assert tr.ref_alen == na.aln_positions[-1].ref_pos - na.aln_positions[0].ref_pos + 1
            assert tr.qry_alen == na.aln_positions[-1].qry_pos - na.aln_positions[0].qry_pos + 1        

        cur_scaf = EmptyAlignment(str(ref_dict[ref].seq).lower())
        assert cur_scaf.get_aligned_str()[0] == str(ref_dict[ref].seq).lower()
        assert cur_scaf.get_aligned_str()[1] == '?' * len(ref_dict[ref])
        
        for na in reversed(nuc_alns):
            # print '%d-%d %d' % (na.aln_positions[0].ref_pos, na.aln_positions[-1].ref_pos, na.aln_positions[-1].ref_pos - na.aln_positions[0].ref_pos)
            cur_scaf = cur_scaf.combine_alignments(na)
            assert cur_scaf.get_aligned_str()[0].replace('.','') == str(ref_dict[ref].seq).lower()
            # print cur_scaf.get_aligned_str()[1]
        
        scaffolds[ref] = cur_scaf
    
    # Output scaffolds as FASTA
    with open(out_scaffold, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            s = scaffolds[ref].get_scaffold()
            print >>outh, '>%s\n%s' % (n, wrap(s))
    
    # Output alignments for other pipeline stages
    with open(out_aln, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            print >>outh, '>%s' % n
            print >>outh, '\n'.join('\t'.join(row) for row in scaffolds[ref].as_matrix())
    
    if not keep_tmp:
        remove_tempdir(tempdir, 'assemble_scaffold')
    
    return out_scaffold

def assemble_scaffold(contigs_fa=None, ref_fa=None, outdir='.',
                      seqname='scaffold',
                      keep_tmp=False, debug=False):
    new_assemble_scaffold(contigs_fa=contigs_fa, ref_fa=ref_fa, outdir=outdir, seqname=seqname, keep_tmp=keep_tmp, debug=debug)
    old_assemble_scaffold(contigs_fa=contigs_fa, ref_fa=ref_fa, outdir=outdir, seqname=seqname, keep_tmp=keep_tmp, debug=debug)    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scaffold contigs using reference sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
