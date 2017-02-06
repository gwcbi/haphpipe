#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse

from utils.sysutils import PipelineStepError, check_dependency
from utils.sysutils import existing_file, existing_dir, command_runner, args_params
from utils.sysutils import create_tempdir, remove_tempdir
from utils.sequtils import assembly_stats

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
    
    group2 = parser.add_argument_group('Assembly options')
    parser.add_argument('--min_contig_length', type=int,
                        help='Minimum assembled contig length to report')    
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--max_memory', type=int,
                        help='Maximum memory to use (in GB)')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Additional options')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')                        
    parser.set_defaults(func=assemble_denovo)

def assemble_denovo(fq1=None, fq2=None, fqU=None, outdir='.',
        min_contig_length=None,
        ncpu=1, max_memory=50, keep_tmp=False, debug=False,
    ):
    """ Assemble reads using Trinity
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
    
    # Check dependencies
    check_dependency('Trinity')
    
    # Outputs
    out1 = os.path.join(outdir, 'contigs.fa')

    # Temporary directory
    tempdir = create_tempdir('tmp_assemble_trinity')        
    
    # Trinity command
    cmd1 = ['Trinity',
        '--CPU', '%d' % ncpu,
        '--max_memory', '%dG' % max_memory,
        '--seqType', 'fq',
        '--output', tempdir,
    ]
    if input_reads in ['paired', 'both', ]:
        cmd1 += ['--left', fq1, '--right', fq2, ]
    elif input_reads in ['single', 'both', ]:
        cmd1 += ['--single', fqU, ]
    
    # Copy command
    cmd2 = ['cp',
        os.path.join(tempdir, 'Trinity.fasta'),
        out1,
    ]

    command_runner([cmd1, cmd2, ], 'assemble_denovo', debug)

    if not keep_tmp:
        remove_tempdir(tempdir, 'assemble_scaffold')       
    
    if os.path.isfile(out1):
        with open(os.path.join(outdir, 'summary.txt'), 'w') as outh:
            assembly_stats(open(out1, 'rU'), outh)
    
    return out1

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assemble reads using Trinity')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))