#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import str
import os
import argparse

from Bio import SeqIO

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils
from haphpipe.utils import alignutils
from haphpipe.utils import gtfparse


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


def stageparser(parser):
    """ Add stage-specific options to argparse parser

    Args:
        parser (argparse.ArgumentParser): ArgumentParser object

    Returns:
        None

    """
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--contigs_fa', type=sysutils.existing_file,
                        required=True,
                        help='Fasta file with assembled contigs')
    group1.add_argument('--ref_fa', type=sysutils.existing_file,
                        required=True,
                        help='''Fasta file with reference genome to scaffold
                                against''')
    group1.add_argument('--ref_gtf', type=sysutils.existing_file,
                        required=True,
                        help='GTF format file containing amplicon regions')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('Scaffold options')
    group2.add_argument('--sample_id', default='sampleXX',
                        help='Sample ID.')
    group2.add_argument('--padding', type=int, default=50,
                        help='Bases to include outside reference annotation.')
    group2.add_argument('--min_contig_len', type=int, default=200,
                        help='Minimum contig length for tiling path')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Additional options')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')

    parser.set_defaults(func=assemble_amplicons)


def assemble_amplicons(
        contigs_fa=None, ref_fa=None, ref_gtf=None, outdir='.',
        sample_id='sampleXX', padding=50, min_contig_len=200,
        keep_tmp=False, quiet=False, logfile=None, debug=False
    ):
    """ Pipeline step to assemble contigs using reference and amplicon regions

    Args:
        contigs_fa (str): Path to fasta file with assembled contigs
        ref_fa (str): Path to reference fasta file
        ref_gtf (str): Path to reference GTF file with amplicons
        outdir (str): Path to output directory
        sample_id (str): Name to append to scaffold sequence
        padding (int): Bases to include outside reference annotation
        min_contig_len (int): Minimum contig length for tiling path
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out_assembly (str): Path to assembled amplicons (FASTA)
        out_summary (str): Path to assembly summary
        out_padded (str): Path to padded output file

    """
    # Check dependencies
    sysutils.check_dependency('nucmer')
    sysutils.check_dependency('delta-filter')
    sysutils.check_dependency('show-tiling')

    # Outputs
    out_assembly = os.path.join(outdir, 'amplicon_assembly.fna')
    out_summary = os.path.join(outdir, 'amplicon_summary.txt')
    out_padded = os.path.join(outdir, 'amplicon_padded.out')
    ###--Uzma--Remove/delete file path. Allows assemble_amplicon to be run more than once?
    if os.path.exists(out_padded): os.unlink(out_padded)

    # Temporary directory
    tempdir = sysutils.create_tempdir(
        'assemble_amplicons', None, quiet, logfile
    )

    # Create fasta file with sequence IDs only (remove decription)
     ###--Uzma--Yields file with tuples of acc# and sequence (descriptions removed in header)?
    tmp_contigs_fa = sequtils.clean_seqnames_file(contigs_fa, tempdir)

    # Load reference sequence(s)
    refseqs = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    
    # For each amplicon, extract the sequence from the reference and scaffold using nucmer    
    amplicon_alignments = []
    amps = [gl for gl in gtfparse.gtf_parser(ref_gtf) if gl.feature == 'amplicon']

    ###--Uzma--What are <>.feature, <>.chrom, and <>.attrs?
    for gl in amps:
        msg = 'Amplicon ref|%s|reg|%s\n' % (gl.chrom, gl.attrs['name'])
        sysutils.log_message(msg, quiet, logfile)
        # Extract reference amplicon
        amp_s = max(0, (gl.start - 1) - padding)
        amp_e = min(len(refseqs[gl.chrom]), gl.end + padding)
        ampseq = refseqs[gl.chrom].seq[amp_s:amp_e]
        amplicon_fa = os.path.join(tempdir, 'subject.fa')
        with open(amplicon_fa, 'w') as outh:
            print('>ref|%s|reg|%s' % (gl.chrom, gl.attrs['name']), file=outh)
            print(sequtils.wrap(str(ampseq)), file=outh)
        
        # Align with nucmer 
        ###--Uzma--Why is fil,til being set equal to the function?
        fil,til = alignutils.align_nucmer(
            tmp_contigs_fa, amplicon_fa, tempdir,
            min_contig_len=min_contig_len,
            quiet=quiet, logfile=logfile, debug=debug
        )

        # Skip everything else if debugging
        if debug: continue

        # Parse tiling and show alignments
        trows = [alignutils.TilingRow(l) for l in open(til, 'rU')]
        if not trows:
            amplicon_alignments.append((gl.chrom, gl.attrs['name'], None))
        else:
            # Initialize alignment
            amp_seq = SeqIO.read(amplicon_fa, 'fasta')
            combined = alignutils.EmptyReferenceAlignment(str(amp_seq.seq).lower())
            for tr in trows:
                out = alignutils.show_aligns(tr.ref, tr.qry, fil)
                for nucaln in alignutils.parse_show_aligns(out):
                    combined = combined.merge_alignments(nucaln)
                    with open(out_padded, 'a') as outh:
                        print('%s\n%s\n%s' % (tr, combined.raln(), combined.qaln()), file=outh)
            amplicon_alignments.append((gl.chrom, gl.attrs['name'], combined))
        
        # Cleanup
        for f in [fil, til, amplicon_fa]:
            if os.path.isfile(f):
                os.unlink(f)

    # Write to output files
    with open(out_assembly, 'w') as outseq, open(out_summary, 'w') as outsum:
        for ref_id, reg, combined in amplicon_alignments:
            amp_id = sequtils.make_seq_id(sid=sample_id, ref=ref_id, reg=reg)
            if combined is None:
                msg1 = '%s\tFAIL\t%d' % (amp_id, 0)
                msg2 = u'%s\tFAIL\t%d\t%s\n' % (amp_id, 0, u"👎🏼")
                if logfile is not None:
                    print(u'%s\tFAIL\t%d\t%s' % (amp_id, 0, u"👎🏼"),
                          file=logfile)
            else:
                scaf, s, e = combined.scaffold2()
                msg1 = '%s\tPASS\t%d' % (amp_id, len(scaf))
                msg2 = u'%s\tPASS\t%d\t%s\n' % (amp_id, len(scaf), u"👍🏼")
                print('>%s' % (amp_id), file=outseq)
                print('%s' % sequtils.wrap(scaf), file=outseq)

            print(msg1, file=outsum)
            sysutils.log_message(msg2, quiet, logfile)

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'assemble_amplicons', quiet, logfile)

    return out_assembly, out_summary, out_padded


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Assemble contigs to amplicon regions.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
