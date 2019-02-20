#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import argparse
import random

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
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')

    group2 = parser.add_argument_group('Sample options')
    group2a = group2.add_mutually_exclusive_group(required=True)
    group2a.add_argument('--nreads', type=int,
                        help='''Number of reads to sample. If greater than the
                                number of reads in file, all reads will be
                                sampled.''')
    group2a.add_argument('--frac', type=float,
                        help='''Fraction of reads to sample, between 0 and 1.
                                Each read has [frac] probability of being
                                sampled, so number of sampled reads is not
                                precisely [frac * num_reads].''')
    group2.add_argument('--seed', type=int,
                        help='''Seed for random number generator.''')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=sample_reads)


def sample_reads(
        fq1=None, fq2=None, fqU=None, outdir='.',
        nreads=None, frac=None, seed=None,
        quiet=False, logfile=None, debug=False,
    ):
    """

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        fqU (str): Path to fastq file with unpaired reads
        outdir (str): Path to output directory
        nreads (int): Number of reads to sample
        frac (float): Fraction of reads to sample
        seed (int): Seed for random number generator
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out1 (str): Path to sampled fastq file with read 1
        out2 (str): Path to sampled fastq file with read 2
        outU (str): Path to sampled fastq file with unpaired reads
    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired"  # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single"  # Single end
    elif fq1 is not None and fq2 is not None and fqU is not None:
        input_reads = "both"
    else:
        msg = "Incorrect combination of reads: "
        msg += "fq1=%s fq2=%s fqU=%s" % (fq1, fq2, fqU)
        raise sysutils.PipelineStepError(msg)

    # Check dependencies
    sysutils.check_dependency('seqtk')

    # Set seed
    seed = seed if seed is not None else random.randrange(1,1000)
    sysutils.log_message(
        '[--- sample_reads ---] Random seed = %d\n' % seed, quiet, logfile
    )

    # Set nreads/frac
    if frac is not None:
        if frac <= 0 or frac > 1:
            raise sysutils.PipelineStepError('--frac must be > 0 and <= 1.')
        frac_arg = '%f' % frac
    else:
        frac_arg = '%d' % nreads

    cmds  = None
    if input_reads == 'single':
        out1 = out2 = None
        outU = os.path.join(outdir, 'sample_U.fastq')
        cmds = [
            [ 'seqtk', 'sample', '-s%d' % seed, fqU, frac_arg, '>', outU, ],
        ]
    elif input_reads == 'paired':
        out1 = os.path.join(outdir, 'sample_1.fastq')
        out2 = os.path.join(outdir, 'sample_2.fastq')
        outU = None
        cmds = [
            ['seqtk', 'sample', '-s%d' % seed, fq1, frac_arg, '>', out1, ],
            ['seqtk', 'sample', '-s%d' % seed, fq2, frac_arg, '>', out2, ],
        ]
    elif input_reads == 'both':
        out1 = os.path.join(outdir, 'sample_1.fastq')
        out2 = os.path.join(outdir, 'sample_2.fastq')
        outU = os.path.join(outdir, 'sample_U.fastq')
        cmds = [
            ['seqtk', 'sample', '-s%d' % seed, fq1, frac_arg, '>', out1, ],
            ['seqtk', 'sample', '-s%d' % seed, fq2, frac_arg, '>', out2, ],
            ['seqtk', 'sample', '-s%d' % seed, fqU, frac_arg, '>', outU, ],
        ]

    sysutils.command_runner(cmds, 'sample_reads', quiet, logfile, debug)
    return out1, out2, outU


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Subsample reads using seqtk.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()