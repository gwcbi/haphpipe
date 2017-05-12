#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse

from ..utils.helpers import guess_encoding
from ..utils.sysutils import PipelineStepError, command_runner
from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.sysutils import create_tempdir, remove_tempdir
# from ..utils.sequtils import wrap


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--fq1', type=existing_file,
                        help='Fastq file with read 1')
    group1.add_argument('--fq2', type=existing_file,
                        help='Fastq file with read 2')
    group1.add_argument('--fqU', type=existing_file,
                        help='Fastq file with unpaired reads')
    group1.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Reference fasta file.')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('Alignment options')
    group2.add_argument('--bt2_preset', 
                        choices=['very-fast', 'fast', 'sensitive', 'very-sensitive',
                                 'very-fast-local', 'fast-local', 'sensitive-local',
                                 'very-sensitive-local',],
                        help='Bowtie2 preset to use')
    group2.add_argument('--rgid',
                        help='Read group ID')
    group2.add_argument('--no_realign', action='store_true',
                        help='Do not realign indels')
    group2.add_argument('--no_markdup', action='store_true',
                        help='Do not remove duplicates')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--encoding',
                        help='Quality score encoding')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=align_reads)

def align_reads(fq1=None, fq2=None, fqU=None, ref_fa=None, outdir='.',
        bt2_preset='sensitive-local', rgid='test', no_realign=None, no_markdup=None,
        ncpu=1, encoding=None, keep_tmp=False, debug=False,
    ):
    """ Pipeline step to align reads
    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired" # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single" # Single end
    elif fq1 is not None and fq2 is not None and fqU is not None:
        input_reads = "both"
    else:
        msg = "Incorrect combination of reads: fq1=%s fq2=%s fqU=%s" % (fq1, fq2, fqU)
        raise PipelineStepError(msg)
    
    if encoding is None:
        encoding = guess_encoding(fq1) if fq1 is not None else guess_encoding(fqU)
    
    # Check dependencies
    check_dependency('bowtie2')
    check_dependency('samtools')
    check_dependency('picard')
    check_dependency('gatk')
    
    # Outputs
    out_aligned = os.path.join(outdir, 'aligned.bam')
    
    # Temporary directory
    tempdir = create_tempdir('align_reads')
    
    # Copy and index initial reference
    curref = os.path.join(tempdir, 'initial.fasta')
    cmd1 = ['cp', ref_fa, curref]
    cmd2 = ['samtools', 'faidx', curref]
    cmd3 = ['picard', 'CreateSequenceDictionary', 
            'R=%s' % curref, 'O=%s' % os.path.join(tempdir, 'initial.dict')]
    cmd4 = ['bowtie2-build', curref, os.path.join(tempdir, 'initial')]
    command_runner([cmd1,cmd2,cmd3,cmd4], 'align_reads:index_ref', debug)
    
    # Align with bowtie2
    out_bt2 = os.path.join(outdir, 'bowtie2.out')
    cmd5 = [
        'bowtie2',
        '-p', '%d' % ncpu,
        '--phred33' if encoding=="Phred+33" else '--phred64',
        '--no-unal',
        '--rg-id', rgid,
        '--rg', 'SM:%s' % rgid,
        '--rg', 'LB:1',
        '--rg', 'PU:1',
        '--rg', 'PL:illumina',
        '--%s' % bt2_preset,
        '-x', '%s' % os.path.join(tempdir, 'initial'),
    ]
    if input_reads in ['paired', 'both', ]:
        cmd5 += ['-1', fq1, '-2', fq2,]
    elif input_reads in ['single', 'both', ]:
        cmd5 += ['-U', fqU, ]
    cmd5 += ['-S', os.path.join(tempdir, 'tmp.sam'), ]
    cmd5 += ['2>&1', '|', 'tee', out_bt2, ]
    cmd6a = ['samtools', 'view', '-bS', os.path.join(tempdir, 'tmp.sam'), '>', os.path.join(tempdir, 'unsorted.bam'),]
    cmd6b = ['samtools', 'sort', os.path.join(tempdir, 'unsorted.bam'), os.path.join(tempdir, 'sorted'),]
    cmd6c = ['samtools', 'index',  os.path.join(tempdir, 'sorted.bam'),]
    cmd7 = ['rm', '-f', os.path.join(tempdir, 'tmp.sam'), os.path.join(tempdir, 'unsorted.bam'), ]
    command_runner([cmd5,cmd6a,cmd6b,cmd6c,cmd7], 'align_reads:align', debug)
    
    cur_bam = os.path.join(tempdir, 'sorted.bam')
    
    if no_markdup:
        print >>sys.stderr, '[--- Skipping MarkDuplicates ---]'
    else:
        # MarkDuplicates
        cmd8 = [
            'picard', 'MarkDuplicates',
            'REMOVE_DUPLICATES=true',
            'CREATE_INDEX=true',
            'M=%s' % os.path.join(tempdir, 'rmdup.metrics.txt'),
            'I=%s' % cur_bam,
            'O=%s' % os.path.join(tempdir, 'rmdup.bam'),
        ]
        command_runner([cmd8,], 'align_reads:markdups', debug)
        cur_bam = os.path.join(tempdir, 'rmdup.bam')
    
    if no_realign:
        print >>sys.stderr, '[--- Skipping realignment ---]'
    else:
        # RealignerTargetCreator
        cmd9 = [
            'gatk', '-T', 'RealignerTargetCreator',
            '-I', cur_bam,
            '-R', curref,
            '-o', os.path.join(tempdir, 'tmp.intervals'),
        ]
        # IndelRealigner
        cmd10 = [
            'gatk', '-T', 'IndelRealigner',
            '-maxReads', '1000000',
            '-dt', 'NONE',
            '-I', cur_bam,
            '-R', curref,
            '-targetIntervals', os.path.join(tempdir, 'tmp.intervals'),
            '-o', os.path.join(tempdir, 'realign.bam')
        ]
        command_runner([cmd9,cmd10,], 'align_reads:realign', debug)
        cur_bam = os.path.join(tempdir, 'realign.bam')
    
    # Check that cur_bam was created
    if not os.path.exists(cur_bam):
        msg = "BAM does not exist: %s" % cur_bam
        raise PipelineStepError(msg)
    
    cmd11a = ['rm', '-f', out_aligned, ]
    cmd11b = ['mv', cur_bam, out_aligned, ]
    cmd11c = ['samtools', 'index', out_aligned, ]
    command_runner([cmd11a, cmd11b, cmd11c, ], 'align_reads:copy', debug)

    if not keep_tmp:
        remove_tempdir(tempdir, 'assemble_scaffold')
    
    return out_aligned

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Refine assembly')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
