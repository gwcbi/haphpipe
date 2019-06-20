#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import os
import argparse

from haphpipe.utils import helpers
from haphpipe.utils import sysutils


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
    group1.add_argument('--fq1', type=sysutils.existing_file,
                        help='Fastq file with read 1')
    group1.add_argument('--fq2', type=sysutils.existing_file,
                        help='Fastq file with read 2')
    group1.add_argument('--fqU', type=sysutils.existing_file,
                        help='Fastq file with unpaired reads')
    group1.add_argument('--ref_fa', type=sysutils.existing_file, required=True,
                        help='Reference fasta file.')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('Alignment options')
    group2.add_argument('--bt2_preset', default='sensitive-local',
                        choices=['very-fast', 'fast', 'sensitive',
                                 'very-sensitive', 'very-fast-local',
                                 'fast-local', 'sensitive-local',
                                 'very-sensitive-local',],
                        help='Bowtie2 preset')
    group2.add_argument('--sample_id', default='sampleXX',
                        help='Sample ID. Used as read group ID in BAM')
    group2.add_argument('--no_realign', action='store_true',
                        help='Do not realign indels')
    group2.add_argument('--remove_duplicates', action='store_true',
                        help='''Remove duplicates from final alignment.
                                Otherwise duplicates are marked but not
                                removed.''')
    group2.add_argument('--encoding',
                        choices=['Phred+33', 'Phred+64'],
                        help='Quality score encoding')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int, default=1,
                        help='Number of CPUs to use')
    group3.add_argument('--xmx', type=int,
                        default=sysutils.get_java_heap_size(),
                        help='Maximum heap size for Java VM, in GB.')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=align_reads)


def align_reads(
        fq1=None, fq2=None, fqU=None, ref_fa=None, outdir='.',
        bt2_preset='sensitive-local', sample_id='sampleXX',
        no_realign=False, remove_duplicates=False, encoding=None,
        ncpu=1, xmx=sysutils.get_java_heap_size(),
        keep_tmp=False, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to align reads

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        fqU (str): Path to fastq file with unpaired reads
        ref_fa (str): Path to reference fasta file
        outdir (str): Path to output directory
        bt2_preset (str): Bowtie2 preset to use for alignment
        sample_id (str): Read group ID
        no_realign (bool): Do not realign indels
        remove_duplicates (bool): Remove duplicates from final alignment
        encoding (str): Quality score encoding
        ncpu (int): Number of CPUs to use
        xmx (int): Maximum heap size for JVM in GB
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out_aligned (str): Path to aligned BAM file
        out_bt2 (str): Path to bowtie2 report

    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired" # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single" # Single end
    elif fq1 is not None and fq2 is not None and fqU is not None:
        input_reads = "both"
    else:
        msg = "Incorrect combination of reads: "
        msg += "fq1=%s fq2=%s fqU=%s" % (fq1, fq2, fqU)
        raise sysutils.PipelineStepError(msg)
    
    if encoding is None:
        if input_reads == 'single':
            encoding = helpers.guess_encoding(fqU)
        else:
            encoding = helpers.guess_encoding(fq1)

    # Check dependencies
    sysutils.check_dependency('bowtie2')
    sysutils.check_dependency('samtools')
    sysutils.check_dependency('picard')

    # Identify correct command for GATK
    GATK_BIN = sysutils.determine_dependency_path(['gatk', 'gatk3'])

    # Set JVM heap argument (for GATK)
    JAVA_HEAP = '_JAVA_OPTIONS="-Xmx%dg"' % xmx

    # Outputs
    out_aligned = os.path.join(outdir, 'aligned.bam')
    out_bt2 = os.path.join(outdir, 'aligned.bt2.out')
    
    # Temporary directory
    tempdir = sysutils.create_tempdir('align_reads', None, quiet, logfile)
    
    # Copy and index initial reference
    curref = os.path.join(tempdir, 'initial.fasta')
    cmd1 = ['cp', ref_fa, curref]
    cmd2 = ['samtools', 'faidx', curref]
    cmd3 = ['picard', 'CreateSequenceDictionary', 
            'R=%s' % curref, 'O=%s' % os.path.join(tempdir, 'initial.dict')]
    cmd4 = ['bowtie2-build', curref, os.path.join(tempdir, 'initial')]
    sysutils.command_runner(
        [cmd1, cmd2, cmd3, cmd4], 'align_reads:index', quiet, logfile, debug
    )
    
    # Align with bowtie2
    cmd5 = [
        'bowtie2',
        '-p', '%d' % ncpu,
        '--phred33' if encoding=="Phred+33" else '--phred64',
        '--no-unal',
        '--rg-id', sample_id,
        '--rg', 'SM:%s' % sample_id,
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
    cmd6a = [
        'samtools', 'view', '-u', os.path.join(tempdir, 'tmp.sam'), '|',
        'samtools', 'sort', '>', os.path.join(tempdir, 'sorted.bam'),
    ]
    # cmd6b = ['samtools', 'sort', os.path.join(tempdir, 'unsorted.bam'),
    # os.path.join(tempdir, 'sorted'),]
    cmd6c = ['samtools', 'index',  os.path.join(tempdir, 'sorted.bam'),]
    # cmd7 = ['rm', '-f', os.path.join(tempdir, 'tmp.sam'), ]
    #  os.path.join(tempdir, 'unsorted.bam'), ]
    sysutils.command_runner(
        [cmd5, cmd6a, cmd6c,], 'align_reads:align', quiet, logfile, debug
    )
    
    cur_bam = os.path.join(tempdir, 'sorted.bam')
    
    if remove_duplicates:
        sysutils.log_message('[--- Removing duplicates ---]', quiet, logfile)
    else:
        sysutils.log_message('[--- Marking duplicates ---]', quiet, logfile)

    # MarkDuplicates
    cmd8 = [
        'picard', 'MarkDuplicates',
        'CREATE_INDEX=true',
        'USE_JDK_DEFLATER=true',
        'USE_JDK_INFLATER=true',
        'M=%s' % os.path.join(tempdir, 'rmdup.metrics.txt'),
        'I=%s' % cur_bam,
        'O=%s' % os.path.join(tempdir, 'rmdup.bam'),
    ]
    if remove_duplicates:
        cmd8 += ['REMOVE_DUPLICATES=true', ]
    sysutils.command_runner(
        [cmd8,], 'align_reads:markdups', quiet, logfile, debug
    )
    cur_bam = os.path.join(tempdir, 'rmdup.bam')
    
    if no_realign:
        print('[--- Skipping realignment ---]', file=sys.stderr)
    else:
        # RealignerTargetCreator
        cmd9 = [
            JAVA_HEAP, GATK_BIN, '-T', 'RealignerTargetCreator',
            '-I', cur_bam,
            '-R', curref,
            '-o', os.path.join(tempdir, 'tmp.intervals'),
        ]
        # IndelRealigner
        cmd10 = [
            JAVA_HEAP, GATK_BIN, '-T', 'IndelRealigner',
            '--use_jdk_deflater', '--use_jdk_inflater',
            '-maxReads', '1000000',
            '-dt', 'NONE',
            '-I', cur_bam,
            '-R', curref,
            '-targetIntervals', os.path.join(tempdir, 'tmp.intervals'),
            '-o', os.path.join(tempdir, 'realign.bam')
        ]
        sysutils.command_runner(
            [cmd9, cmd10, ], 'align_reads:realign', quiet, logfile, debug
        )
        cur_bam = os.path.join(tempdir, 'realign.bam')
    
    # Check that cur_bam was created
    if not os.path.exists(cur_bam):
        msg = "BAM does not exist: %s" % cur_bam
        raise sysutils.PipelineStepError(msg)
    
    cmd11a = ['rm', '-f', out_aligned, ]
    cmd11b = ['mv', cur_bam, out_aligned, ]
    cmd11c = ['samtools', 'index', out_aligned, ]
    sysutils.command_runner(
        [cmd11a, cmd11b, cmd11c, ], 'align_reads:copy', quiet, logfile, debug
    )

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'align_reads', quiet, logfile)
    
    return out_aligned, out_bt2


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='''Align reads to reference.''',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()

