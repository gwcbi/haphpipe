#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from subprocess import check_output

from utils.sysutils import PipelineStepError, check_dependency
from utils.sysutils import existing_file, existing_dir, command_runner, args_params

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--fq1', type=existing_file,
                        help='Fastq file with read 1')
    group1.add_argument('--fq2', type=existing_file,
                        help='Fastq file with read 1')
    group1.add_argument('--fqU', type=existing_file,
                        help='Fastq file with unpaired reads')              
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('BLESS settings')
    group2.add_argument('--kmerlength', type=int,
                        help='''Length of k-mers.''')
    group2.add_argument('--extend', type=int,
                        help='''Read extension amount.''')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=ec_reads)


def ec_reads(fq1=None, fq2=None, fqU=None, outdir='.',
             kmerlength=31, extend=None,
             ncpu=1, debug=False,
    ):
    """ Error correct reads using BLESS
    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired" # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single" # Single end
    else:
        msg = "Incorrect combination of reads: fq1=%s fq2=%s fqU=%s" % (fq1, fq2, fqU)
        raise PipelineStepError(msg)
    
    # Check dependencies
    check_dependency('bless')
    
    # NOTE: bless should be called with the full path
    cmd1 = [check_output('which bless', shell=True).strip(), ]
    if input_reads is "single":
        # Outputs
        out1 = out2 = None
        outU = os.path.join(outdir, 'bless.corrected.fastq')
        cmd1 += ['-read', fqU, ]
    elif input_reads is "paired":
        # Outputs
        out1 = os.path.join(outdir, 'bless.1.corrected.fastq')
        out2 = os.path.join(outdir, 'bless.2.corrected.fastq')
        outU = None
        cmd1 += ['-read1', fq1, '-read2', fq2, ]
    # Include options
    cmd1 += [
        '-prefix', os.path.join(outdir, 'bless'),
        '-kmerlength', '%d' % kmerlength,
    ]
    if extend is not None:
        cmd1 += ['-extend', '%d' % extend, ]
    
    command_runner([cmd1,], 'ec_reads', debug)
    return out1, out2, outU


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Error correct reads using BLESS')
    ec_reads_parser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))