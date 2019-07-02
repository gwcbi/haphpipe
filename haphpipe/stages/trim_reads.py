#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import os
import argparse

from haphpipe.utils import helpers
from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import PipelineStepError
from haphpipe.utils.sysutils import MissingRequiredArgument


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


## Constants
# Set default trimming parameters
TRIMMERS = [
    "LEADING:3",
    "TRAILING:3",
    "SLIDINGWINDOW:4:15",
    "MINLEN:36",
]


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
                        help='Fastq file with read 1')
    group1.add_argument('--fqU', type=sysutils.existing_file,
                        help='Fastq file with unpaired reads')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('Trimmomatic options')
    group2.add_argument('--adapter_file',
                        help='Adapter file')
    group2.add_argument('--trimmers', action='append', default=TRIMMERS,
                        help='Trim commands for trimmomatic')
    group2.add_argument('--encoding',
                        choices=['Phred+33', 'Phred+64'],
                        help='Quality score encoding')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int, default=1,
                        help='Number of CPU to use')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')

    parser.set_defaults(func=trim_reads)


def trim_reads(
        fq1=None, fq2=None, fqU=None, outdir=".",
        adapter_file=None, trimmers=TRIMMERS, encoding=None,
        ncpu=1, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to trim reads

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        fqU (str): Path to fastq file with unpaired reads
        outdir (str): Path to output directory
        adapter_file (str): Path to adapter file (fasta)
        trimmers (`list` of `str`): Trim commands for trimmomatic
        encoding (str): Quality score encoding
        ncpu (int): Number of CPUs to use
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out1 (str): Path to trimmed fastq file with read 1
        out2 (str): Path to trimmed fastq file with read 2
        outU (str): Path to trimmed fastq file with unpaired reads
        out_summary (str): Path to summary file
    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired" # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single" # Single end
    else:
        msg = "incorrect input reads; requires either "
        msg += "(--fq1 and --fq2) OR (--fqU)"
        raise MissingRequiredArgument(msg)
        
    """ There are two different ways to call Trimmomatic. If using modules on
        C1, the path to the jar file is stored in the "$Trimmomatic"
        environment variable. Otherwise, if using conda, the "trimmomatic"
        script is in PATH.
    """
    # Check dependencies
    try:
        sysutils.check_dependency('trimmomatic')
        cmd1 = ['trimmomatic']
    except PipelineStepError as e:
        if 'Trimmomatic' in os.environ:
            cmd1 = ['java', '-jar', '$Trimmomatic']
        else:
            raise e

    # Get encoding
    if encoding is None:
        if input_reads == 'single':
            encoding = helpers.guess_encoding(fqU)
        else:
            encoding = helpers.guess_encoding(fq1)

    # Outputs for both single and paired
    out_summary = os.path.join(outdir, 'trimmomatic_summary.out')
    outU = os.path.join(outdir, 'trimmed_U.fastq')

    if input_reads is 'single':
        # Outputs
        out1 = out2 = None
        # Trimmomatic command
        cmd1 += [
            'SE',
            '-threads', '%d' % ncpu,
            '-phred33' if encoding == "Phred+33" else '-phred64',
            '-summary', out_summary,
            fqU,
            outU,
        ]
        # Specify trimming steps
        if adapter_file is not None:
            adapter_file = adapter_file.replace('PE', 'SE')
            cmd1.append("ILLUMINACLIP:%s:2:30:10" % adapter_file)
        cmd1 += trimmers
        
        # Run command
        sysutils.command_runner(
            [cmd1,], 'trim_reads', quiet, logfile, debug
        )
        return out1, out2, outU
    elif input_reads is 'paired':
        # Outputs
        out1 = os.path.join(outdir, 'trimmed_1.fastq')
        out2 = os.path.join(outdir, 'trimmed_2.fastq')
        tmp1U = os.path.join(outdir, 'tmp1U.fq')
        tmp2U = os.path.join(outdir, 'tmp2U.fq')
        # Trimmomatic command
        cmd1 += [
            'PE',
            '-threads', '%d' % ncpu,
            '-phred33' if encoding == "Phred+33" else '-phred64',
            '-summary', out_summary,
            fq1, fq2,
            out1, tmp1U,
            out2, tmp2U,
        ]
        # Specify trimming steps
        if adapter_file is not None:
            cmd1.append("ILLUMINACLIP:%s:2:30:10" % adapter_file)
        cmd1 += trimmers
        
        # Concat files command
        cmd2 = ['cat', tmp1U, tmp2U, '>>',  outU, ]
        cmd3 = ['rm', '-f', tmp1U, tmp2U, ]
        
        # Run commands
        sysutils.command_runner(
            [cmd1, cmd2, cmd3, ], 'trim_reads', quiet, logfile, debug
        )
        return out1, out2, outU, out_summary


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Trim reads using Trimmomatic.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    try:
        args.func(**sysutils.args_params(args))
    except MissingRequiredArgument as e:
        parser.print_usage()
        print('error: %s' % e, file=sys.stderr)


if __name__ == '__main__':
    console()
