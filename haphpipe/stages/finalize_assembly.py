#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import absolute_import

from builtins import str
import os
import argparse
import shutil

from Bio import SeqIO

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils
from haphpipe.stages import align_reads
from haphpipe.stages import call_variants


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--fq1', type=sysutils.existing_file,
                        help='Fastq file with read 1')
    group1.add_argument('--fq2', type=sysutils.existing_file,
                        help='Fastq file with read 1')
    group1.add_argument('--fqU', type=sysutils.existing_file,
                        help='Fastq file with unpaired reads')
    group1.add_argument('--ref_fa', type=sysutils.existing_file, required=True,
                        help='Consensus fasta file')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('Fix consensus options')
    group2.add_argument('--bt2_preset', default='very-sensitive',
                        choices=['very-fast', 'fast', 'sensitive', 'very-sensitive',
                                 'very-fast-local', 'fast-local', 'sensitive-local',
                                 'very-sensitive-local',],
                        help='Bowtie2 preset to use')
    group2.add_argument('--sample_id', default='sampleXX',
                        help='Sample ID')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=finalize_assembly)


def finalize_assembly(fq1=None, fq2=None, fqU=None, ref_fa=None, outdir='.',
        bt2_preset='very-sensitive', sample_id='sampleXX',
        ncpu=1,
        keep_tmp=False, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to finalize consensus
    """
    # Outputs
    out_ref = os.path.join(outdir, 'final.fna')
    out_aligned = os.path.join(outdir, 'final.bam')
    out_bt2 = os.path.join(outdir, 'final_bt2.out')
    out_vcf = os.path.join(outdir, 'final.vcf.gz')

    # Temporary directory
    tempdir = sysutils.create_tempdir(
        'finalize_assembly', None, quiet, logfile
    )

    # Copy reference and rename sequences
    with open(out_ref, 'w') as outh:
        for s in SeqIO.parse(ref_fa, 'fasta'):
            if sample_id == 'sampleXX':
                dline = s.description
            else:
                dline = sequtils.update_seq_id(s.id, sample_id)
                dline += ' SM:%s' % sample_id
            print('>%s' % dline, file=outh)
            print(sequtils.wrap(str(s.seq).upper()), file=outh)

    # Align to reference
    tmp_aligned, tmp_bt2 = align_reads.align_reads(
        fq1=fq1, fq2=fq2, fqU=fqU, ref_fa=out_ref, outdir=tempdir,
        bt2_preset=bt2_preset, sample_id=sample_id,
        ncpu=ncpu,
        keep_tmp=keep_tmp, quiet=quiet, logfile=logfile, debug=debug,
    )

    # Call variants
    tmp_vcf = call_variants.call_variants(
        aln_bam=tmp_aligned, ref_fa=out_ref, outdir=tempdir,
        emit_all=False,
        ncpu=ncpu,
        keep_tmp=keep_tmp, quiet=quiet, logfile=logfile, debug=debug,
    )

    shutil.copy(tmp_aligned, out_aligned)
    shutil.copy(tmp_bt2, out_bt2)
    shutil.copy(tmp_vcf, out_vcf)

    # Index BAM and VCF
    cmds = [
        ['tabix', out_vcf],
        ['samtools', 'index', out_aligned],
    ]
    sysutils.command_runner(
        cmds, 'finalize_assembly', quiet, logfile, debug,
    )

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'finalize_assembly', quiet, logfile)

    return out_ref, out_aligned, out_vcf, out_bt2


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='''Finalize consensus sequence.''',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))

if __name__ == '__main__':
    console()
