#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse

from ..utils.helpers import guess_encoding
from ..utils.sysutils import PipelineStepError, check_dependency
from ..utils.sysutils import existing_file, existing_dir, command_runner, args_params

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--fq1', type=existing_file, required=True,
                        help='Fastq file with read 1')
    group1.add_argument('--fq2', type=existing_file, required=True,
                        help='Fastq file with read 1')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('FLAsh settings')
    group2.add_argument('--min_overlap', type=int,
                        help='''The minimum required overlap length between two
                                reads to provide a confident overlap.''')
    group2.add_argument('--max_overlap', type=int,
                        help='''Maximum overlap length expected in approximately
                                90 pct of read pairs.''')
    group2.add_argument('--allow_outies', action='store_true',
                        help='''Also try combining read pairs in the "outie"
                                orientation''')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--encoding',
                        help='Quality score encoding')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=join_reads)

def join_reads(fq1=None, fq2=None, outdir=".",
        min_overlap=None, max_overlap=None, allow_outies=None,
        ncpu=1, encoding=None, debug=False,
    ):
    """ Join reads using FLASH
    """
    # Check inputs
    if fq1 is not None and fq2 is not None:
        pass # Both are present
    else:
        msg = "Incorrect combination of reads: fq1=%s fq2=%s" % (fq1, fq2)
        raise PipelineStepError(msg)
    
    # Check for executable
    check_dependency('flash')
    
    if encoding is None:
        encoding = guess_encoding(fq1)
    
    # Outputs
    outU = os.path.join(outdir, 'flash.extendedFrags.fastq')
    out1 = os.path.join(outdir, 'flash.notCombined_1.fastq')
    out2 = os.path.join(outdir, 'flash.notCombined_2.fastq')
    
    # Flash command
    cmd1 = [
        'flash',
        '-t', '%d' % ncpu,
        '-d', '%s' % outdir,
        '-o', 'flash',
    ]
    if encoding != "Phred+33":
        cmd1 += ['-p', '64']
    if min_overlap is not None:
        cmd1 += ['-m', '%d' % min_overlap]
    if max_overlap is not None:
        cmd1 += ['-M', '%d' % max_overlap]
    if allow_outies is True:
        cmd1 += ['-O']        
    cmd1 += [fq1, fq2]
    
    command_runner([cmd1,], 'join_reads', debug)
    return out1, out2, outU

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Join reads using FLASH')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
