#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import argparse

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
    group1.add_argument('--aln_bam', type=sysutils.existing_file,
                        required=True,
                        help='Alignment file.')    
    group1.add_argument('--ref_fa',
                        type=sysutils.existing_file,
                        required=True,
                        help='Reference fasta file.')
    group1.add_argument('--outdir', type=sysutils.existing_dir,
                        default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('Variant calling options')    
    group2.add_argument('--emit_all', action='store_true',
                        help='Output calls for all sites.')
    group2.add_argument('--min_base_qual', type=int, default=15,
                        help='''Minimum base quality required to consider a
                                base for calling.''')
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
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
    parser.set_defaults(func=call_variants)


def call_variants(
        aln_bam=None, ref_fa=None, outdir='.',
        emit_all=False, min_base_qual=15,
        ncpu=1, xmx=sysutils.get_java_heap_size(),
        keep_tmp=False, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to call variants

    Args:
        aln_bam (str): Path to alignment file (BAM)
        ref_fa (str): Path to reference fasta file
        outdir (str): Path to output directory
        emit_all (bool): Output calls for all sites
        min_base_qual (int): Minimum base quality for calling
        ncpu (int): Number of CPUs to use
        xmx (int): Maximum heap size for JVM in GB
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out_vcf (str): Path to output VCF

    """
    # Check dependencies
    sysutils.check_dependency('samtools')
    sysutils.check_dependency('picard')

    # Identify correct command for GATK
    GATK_BIN = sysutils.determine_dependency_path(['gatk', 'gatk3'])

    # Set JVM heap argument (for GATK)
    JAVA_HEAP = '_JAVA_OPTIONS="-Xmx%dg"' % xmx

    # Outputs
    out_vcf = os.path.join(outdir, 'variants.vcf.gz')    
    
    # Temporary directory
    tempdir = sysutils.create_tempdir('call_variants', None, quiet, logfile)
    
    # Copy and index initial reference
    curref = os.path.join(tempdir, 'initial.fasta')
    cmd1 = ['cp', ref_fa, curref]
    cmd2 = ['samtools', 'faidx', curref]
    cmd3 = ['picard', 'CreateSequenceDictionary', 
            'R=%s' % curref, 'O=%s' % os.path.join(tempdir, 'initial.dict')]
        
    # UnifiedGenotyper
    cmd4 = [JAVA_HEAP, GATK_BIN, '-T', 'UnifiedGenotyper',
        '--use_jdk_deflater', '--use_jdk_inflater',
        '--num_threads', '%d' % ncpu,
        '-gt_mode', 'DISCOVERY',
        '-glm', 'BOTH',
        '--baq', 'OFF',
        '--useOriginalQualities',
        '-dt', 'NONE',
        '--min_base_quality_score', '%d' % min_base_qual,
        '-ploidy', '4',
        '-I', aln_bam,
        '-R', curref,
        '-o', out_vcf,
    ]
    if emit_all:
        cmd4 += ['-out_mode', 'EMIT_ALL_SITES']

    sysutils.command_runner(
        [cmd1, cmd2, cmd3, cmd4,], 'call_variants:GATK', quiet, logfile, debug
    )

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'call_variants:GATK', quiet, logfile)

    return out_vcf


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Call variants.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
