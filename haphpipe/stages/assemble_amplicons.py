#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse

from Bio import SeqIO

from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.sysutils import create_tempdir, remove_tempdir
from ..utils.sequtils import clean_seqnames, wrap
from ..utils import alignobj
from ..utils.alignutils import align_nucmer, show_aligns, parse_show_aligns
from ..utils.gtfparse import gtf_parser

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--contigs_fa', type=existing_file, required=True,
                        help='Fasta file with assembled contigs')
    group1.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Fasta file with reference genome to scaffold against')
    group1.add_argument('--ref_gtf', type=existing_file, required=True,
                        help='GTF format file containing amplicon regions')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('Scaffold options')
    group2.add_argument('--seqname',
                        help='Name to append to scaffold sequence.')
    group2.add_argument('--padding', type=int,
                        help='Name to append to scaffold sequence.')
        
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Additional options')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=assemble_amplicons)

def assemble_amplicons(contigs_fa=None, ref_fa=None, ref_gtf=None, outdir='.',
                      seqname='newseq', padding=50,
                      keep_tmp=False, debug=False):
    """ Scaffold contigs using reference sequence
    """
    # Check dependencies
    check_dependency('nucmer')
    check_dependency('delta-filter')
    check_dependency('show-tiling')

    # Outputs
    out_padded = os.path.join(outdir, 'padded.out')
    if os.path.exists(out_padded): os.unlink(out_padded)
    out_assembly = os.path.join(outdir, 'assembly.fa')
    out_summary = os.path.join(outdir, 'summary.txt')
    
    # Temporary directory
    tempdir = create_tempdir('assemble_amplicons')

    # Create fasta file with sequence IDs only (remove decription)
    tmp_contigs_fa = os.path.join(tempdir, 'query.fna')
    with open(tmp_contigs_fa, 'w') as outh:
        for n,s in clean_seqnames(open(contigs_fa, 'rU')):
            print >>outh, '>%s\n%s' % (n, wrap(s))
    
    # Load reference sequence(s)
    refseqs = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    
    # For each amplicon, extract the sequence from the reference and scaffold using nucmer    
    amplicon_alignments = []
    amps = [gtf_line for gtf_line in gtf_parser(ref_gtf) if gtf_line.feature == 'amplicon']
    for gl in amps:
        print >>sys.stderr, 'Amplicon ref|%s|reg|%s' % (gl.chrom, gl.attrs['name'])
        # Extract reference amplicon
        amp_s = max(0, (gl.start - 1) - padding)
        amp_e = min(len(refseqs[gl.chrom]), gl.end + padding)
        ampseq = refseqs[gl.chrom].seq[amp_s:amp_e]
        amplicon_fa = os.path.join(tempdir, 'subject.fa')
        with open(amplicon_fa, 'w') as outh:
            print >>outh, '>ref|%s|reg|%s' % (gl.chrom, gl.attrs['name'])
            print >>outh, wrap(str(ampseq))
        
        # Align with nucmer 
        fil,til = align_nucmer(amplicon_fa, tmp_contigs_fa, tempdir)
        
        # Parse tiling and show alignments
        trows = [alignobj.TilingRow(l) for l in open(til, 'rU')]
        if not trows:
            amplicon_alignments.append((gl.chrom, gl.attrs['name'], None))
        else:
            # Initialize alignment
            amp_seq = SeqIO.read(amplicon_fa, 'fasta')
            combined = alignobj.EmptyReferenceAlignment(str(amp_seq.seq).lower())
            for tr in trows:
                out = show_aligns(tr.ref, tr.qry, fil)
                for nucaln in parse_show_aligns(out):
                    combined = combined.merge_alignments(nucaln)
                    with open(out_padded, 'a') as outh:
                        print >>outh, '%s\n%s\n%s' % (tr, combined.raln(), combined.qaln())
            amplicon_alignments.append((gl.chrom, gl.attrs['name'], combined))
        
        # Cleanup
        for f in [fil, til, amplicon_fa]:
            if os.path.isfile(f):
                os.unlink(f)

    with open(out_assembly, 'w') as outseq, open(out_summary, 'w') as outsum:
        for chrom, name, combined in amplicon_alignments:
            # chrom, start, end, name, combined = t[0], int(t[1]), int(t[2]), t[3], t[4]
            amp_id = 'sid|%s|ref|%s|reg|%s' % (seqname, chrom, name)
            if combined is None:
                print >>outsum, '%s\tFAIL\t%d' % (amp_id, 0)
                print >>sys.stderr, '%s\tFAIL\t%d\t%s' % (amp_id, 0,"👎🏼")
            else:
                scaf, s, e = combined.scaffold2()
                print >>outsum, '%s\tPASS\t%d' % (amp_id, len(scaf))
                print >>sys.stderr, '%s\tPASS\t%d\t%s' % (amp_id, len(scaf), "👍🏼")
                print >>outseq, '>%s' % (amp_id)
                print >>outseq, '%s' % wrap(scaf)
    
    if not keep_tmp:
        remove_tempdir(tempdir, 'assemble_scaffold')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Scaffold contigs using reference sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))       
