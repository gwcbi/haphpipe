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



__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--aln_bam', type=existing_file, required=True,
                        help='Alignment file.')    
    group1.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Reference fasta file.')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('Variant calling options')    
    group2.add_argument('--emit_all', action='store_true',
                        help='Output calls for all sites.')
    
    group2.add_argument('--min_dp', type=int,
                        help='Minimum depth to call position')
    group2.add_argument('--rgid',
                        help='Read group ID')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=call_variants)

def call_variants(aln_bam=None, ref_fa=None, outdir='.',
        emit_all=False,
        ncpu=1, keep_tmp=False, debug=False,
    ):
    """ Pipeline step to refine assembly
    """
    # Check dependencies
    check_dependency('samtools')
    check_dependency('picard')
    check_dependency('gatk')
    
    # Outputs
    # out_refined = os.path.join(outdir, 'refined.fa')
    out_vcf = os.path.join(outdir, 'variants.vcf.gz')    
    
    # Temporary directory
    tempdir = create_tempdir('call_variants')
    
    # Copy and index initial reference
    curref = os.path.join(tempdir, 'initial.fasta')
    cmd1 = ['cp', ref_fa, curref]
    cmd2 = ['samtools', 'faidx', curref]
    cmd3 = ['picard', 'CreateSequenceDictionary', 
            'R=%s' % curref, 'O=%s' % os.path.join(tempdir, 'initial.dict')]
        
    # UnifiedGenotyper
    cmd11 = [
        'gatk', '-T', 'UnifiedGenotyper',
        '--num_threads', '%d' % ncpu,
        '-gt_mode', 'DISCOVERY',
        '-glm', 'BOTH',
        '--baq', 'OFF',
        '--useOriginalQualities',
        '-dt', 'NONE',
        '--min_base_quality_score', '15',
        '-ploidy', '4',
        '-I', aln_bam,
        '-R', curref,
        '-o', out_vcf,
    ]
    if emit_all:
        cmd11 += ['-out_mode', 'EMIT_ALL_SITES']
    command_runner([cmd1, cmd2, cmd3, cmd11,], 'call_variants:GATK', debug)

    return out_vcf

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Refine assembly')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
