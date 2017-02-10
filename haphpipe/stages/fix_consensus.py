#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse
import shutil

from Bio import SeqIO

from ..utils.helpers import guess_encoding
from ..utils.sysutils import PipelineStepError, command_runner
from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.sysutils import create_tempdir, remove_tempdir
from ..utils.sequtils import wrap, extract_amplicons
from vcf_to_fasta import vcf_to_fasta
from ..utils.alignutils import assemble_to_ref

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
    group1.add_argument('--assembly_fa', type=existing_file, required=True,
                        help='Fasta file with assembly')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('Fix consensus options')
    group2.add_argument('--bt2_preset', 
                        choices=['very-fast', 'fast', 'sensitive', 'very-sensitive',
                                 'very-fast-local', 'fast-local', 'sensitive-local',
                                 'very-sensitive-local',],
                        help='Bowtie2 preset to use')
    group2.add_argument('--rgid',
                        help='Read group ID')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--encoding',
                        help='Quality score encoding')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=fix_consensus)

def fix_consensus(fq1=None, fq2=None, fqU=None, assembly_fa=None, outdir='.',
        bt2_preset='very-sensitive', min_dp=1, rgid='test',
        ncpu=1, encoding=None, keep_tmp=False, debug=False,
    ):
    """ Pipeline step to refine assembly
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
    check_dependency('bowtie2')
    check_dependency('samtools')
    check_dependency('picard')
    check_dependency('gatk')
    
    # Check encoding
    if encoding is None:
        encoding = guess_encoding(fq1) if input_reads == 'paired' else guess_encoding(fqU)
    
    # Temporary directory
    tempdir = create_tempdir('fix_consensus')
     
    # Copy assembly_fa to consensus
    shutil.copyfile(assembly_fa, os.path.join(outdir, 'consensus.fasta'))
    
    # Copy and index consensus reference
    curref = os.path.join(tempdir, 'initial.fasta')
    cmd1 = ['cp', os.path.join(outdir, 'consensus.fasta'), curref]
    cmd2 = ['samtools', 'faidx', curref]
    cmd3 = ['picard', 'CreateSequenceDictionary', 
            'R=%s' % curref, 'O=%s' % os.path.join(tempdir, 'initial.dict')]
    cmd4 = ['bowtie2-build', curref, os.path.join(tempdir, 'initial')]
    command_runner([cmd1,cmd2,cmd3,cmd4], 'fix_consensus:index_ref', debug)
    
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
    cmd5 += ['-S', os.path.join(tempdir, 'unsorted.sam'), ]
    cmd5 += ['2>&1', '|', 'tee', out_bt2, ]
    cmd6 = ['samtools', 'view', '-bS', os.path.join(tempdir, 'unsorted.sam'), '>', os.path.join(tempdir, 'unsorted.bam'),]
    cmd7 = ['samtools', 'sort', os.path.join(tempdir, 'unsorted.bam'), os.path.join(tempdir, 'withdups'),]
    cmd8 = ['samtools', 'index',  os.path.join(tempdir, 'withdups.bam'),]
    cmd9 = ['rm', '-f', os.path.join(tempdir, 'unsorted.sam'), os.path.join(tempdir, 'unsorted.bam'), ]
    command_runner([cmd5,cmd6,cmd7,cmd8,cmd9], 'fix_consensus:align', debug)
    
    # MarkDuplicates
    # For now, just mark duplicates but not remove them
    cmd10a = [
        'picard', 'MarkDuplicates',
        'REMOVE_DUPLICATES=false',
        'CREATE_INDEX=true',
        'M=%s' % os.path.join(tempdir, 'rmdup.metrics.txt'),
        'I=%s' % os.path.join(tempdir, 'withdups.bam'),
        'O=%s' % os.path.join(tempdir, 'all.bam'),
    ]
    cmd10b = ['mv', 
        os.path.join(tempdir, 'all.bai'), 
        os.path.join(tempdir, 'all.bam.bai'),
    ]
    cmd10c = ['rm', '-f', 'withdups.bam*', ]
    command_runner([cmd10a, cmd10b, cmd10c, ], 'fix_consensus:mark_dups', debug)
    
    # Filtering
    # Duplicates were marked with 0x400
    if input_reads == 'paired' or input_reads == 'both':
        # Include (-f):PAIRED,PROPER_PAIR 
        # Exclude (-F): UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY    
        cmdF1 = ['samtools', 'view',
            '-b',
            '-f', '3',
            '-F', '3852',
            os.path.join(tempdir, 'all.bam'),
            '>',
            os.path.join(tempdir, 'f1.bam'),
        ]
    elif input_reads == 'single':
        # Exclude: PAIRED,UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
        cmdF1 = ['samtools', 'view',
            '-b',
            '-F', '3853',
            os.path.join(tempdir, 'all.bam'),
            '>',
            os.path.join(tempdir, 'f1.bam'),
        ]
    cmdF2 = ['samtools', 'index', os.path.join(tempdir, 'f1.bam'), ]
    command_runner([cmdF1, cmdF2, ], 'fix_consensus:filtering', debug)    
    
    # RealignerTargetCreator
    cmd11 = [
        'gatk', '-T', 'RealignerTargetCreator',
        '-I', os.path.join(tempdir, 'f1.bam'),
        '-R', curref,
        '-o', os.path.join(tempdir, 'tmp.intervals'),
    ]
    # IndelRealigner
    cmd12a = [
        'gatk', '-T', 'IndelRealigner',
        '-maxReads', '1000000',
        '-dt', 'NONE',
        '-I', os.path.join(tempdir, 'f1.bam'),
        '-R', curref,
        '-targetIntervals', os.path.join(tempdir, 'tmp.intervals'),
        '-o', os.path.join(tempdir, 'final.bam')
    ]
    cmd12b = ['rm', '-f',  os.path.join(tempdir, 'final.bai'),]
    cmd12c = ['samtools', 'index',  os.path.join(tempdir, 'final.bam'),]
    command_runner([cmd11, cmd12a, cmd12b, cmd12c, ], 'fix_consensus:realign', debug)

    # UnifiedGenotyper
    cmd13 = [
        'gatk', '-T', 'UnifiedGenotyper',
        '--num_threads', '%d' % ncpu,
        '-glm', 'BOTH',
        '--baq', 'OFF',
        '--useOriginalQualities',
        '-dt', 'NONE',
        '-A', 'AlleleBalance',
        '--min_base_quality_score', '15',
        '-ploidy', '4',
        '-I', os.path.join(tempdir, 'final.bam'),
        '-R', curref,
        '-o', os.path.join(outdir, 'variants.ug.vcf.gz'),
    ]
    command_runner([cmd13,], 'fix_consensus:Unified Genotyper', debug)
    
    cmd14 = [
        'freebayes',
         '--min-alternate-fraction', '0.01'
         '--pooled-continuous',
         '--standard-filters',
         '--ploidy 1',
         '--haplotype-length', '0',
         '-f', curref,
         os.path.join(tempdir, 'final.bam'),
    ]
    cmd14 += ['|', 'bcftools', 'view', '-Oz',]
    cmd14 += ['|', 'bcftools', 'filter', '-m', "'+'", '-Oz', '-e' '"QA > 4000"', '-s', "'HQ'",]
    cmd14 += ['|', 'bcftools', 'filter', '-m', "'+'", '-Oz', '-e' '"AO/DP > 0.50"', '-s', "'GT50'",]
    cmd14 += ['|', 'bcftools', 'filter', '-m', "'+'", '-Oz', '-e' '"AO/DP <= 0.20"', '-s', "'LT20'",]
    cmd14 += ['>', os.path.join(outdir, 'variants.fb.vcf.gz'), ]
    # command_runner([cmd14,], 'fix_consensus:Freebayes', debug)      

    command_runner([
        ['cp', os.path.join(tempdir, 'all.bam'), outdir,],
        ['cp', os.path.join(tempdir, 'all.bam.bai'), outdir,],        
        ['cp', os.path.join(tempdir, 'final.bam'), outdir,],
        ['cp', os.path.join(tempdir, 'final.bam.bai'), outdir],
    ], 'fix_consensus:cleanup', debug)
    
    if not keep_tmp:
        remove_tempdir(tempdir, 'fix_consensus')
    
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix consensus sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))

    



