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

# Set default trimming parameters
TRIMMERS = [
    "LEADING:3",
    "TRAILING:3",
    "SLIDINGWINDOW:4:15",
    "MINLEN:36",
]

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
    
    group2 = parser.add_argument_group('Trimmomatic settings')
    group2.add_argument('--adapter_file',
                        help='Adapter file')
    group2.add_argument('--trimmers', action='append',
                        help='Trim commands for trimmomatic')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--encoding',
                        help='Quality score encoding')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=trim_reads)

def trim_reads(fq1=None, fq2=None, fqU=None, outdir=".",
               adapter_file=None, trimmers=TRIMMERS,
               ncpu=1, encoding=None, debug=False,
    ):
    """ Trim reads using Trimmomatic
        
    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired" # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single" # Single end
    else:
        msg = "Incorrect combination of reads: fq1=%s fq2=%s fqU=%s" % (fq1, fq2, fqU)
        raise PipelineStepError(msg)
        
    """ There are two different ways to call Trimmomatic. If using modules on C1, the
        path to the jar file is stored in the "$Trimmomatic" environment variable.
        Otherwise, if using conda, the "trimmomatic" script is in PATH
    """
    # Check dependencies
    try:
        check_dependency('trimmomatic')
        cmd1 = ['trimmomatic']
    except PipelineStepError as e:
        if 'Trimmomatic' in os.environ:
            cmd1 = ['java', '-jar', '$Trimmomatic']
        else:
            raise e
    
    if encoding is None:
        encoding = guess_encoding(fq1) if input_reads == 'paired' else guess_encoding(fqU)
    
    if input_reads is 'single':
        # Assume single-end reads
        out1 = out2 = None
        outU = os.path.join(outdir, 'trimmed_U.fastq')
        
        # Trimmomatic command
        cmd1 += [
            'SE',
            '-threads', '%d' % ncpu,
            '-phred33' if encoding == "Phred+33" else '-phred64',
            '-trimlog', '%s' % os.path.join(outdir, 'trimmomatic.log'),
            fq1,
            outU,
        ]
        # Specify trimming steps
        if adapter_file is not None:
            adapter_file = adapter_file.replace('PE', 'SE')
            cmd1.append("ILLUMINACLIP:%s:2:30:10" % adapter_file)
        cmd1 += trimmers
        
        # Run command
        command_runner([cmd1,], 'trim_reads', debug)
        return out1, out2, outU
    elif input_reads is 'paired':
        # Outputs
        out1 = os.path.join(outdir, 'trimmed_1.fastq')
        out2 = os.path.join(outdir, 'trimmed_2.fastq')
        outU = os.path.join(outdir, 'trimmed_U.fastq')
        tmp1U = os.path.join(outdir, 'tmp1U.fq')
        tmp2U = os.path.join(outdir, 'tmp2U.fq')
        # Trimmomatic command
        cmd1 += [
            'PE',
            '-threads', '%d' % ncpu,
            '-phred33' if encoding == "Phred+33" else '-phred64',
            '-trimlog', '%s' % os.path.join(outdir, 'trimmomatic.log'),
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
        cmd3 = ['rm', '-f', tmp1U, ]
        cmd4 = ['rm', '-f', tmp2U, ]
        
        # Run commands
        command_runner([cmd1,], 'trim_reads', debug)
        return out1, out2, outU

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Trim reads using Trimmomatic')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
