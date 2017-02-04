#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse

from Bio import SeqIO

from utils.sysutils import PipelineStepError, command_runner
from utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from utils.sysutils import create_tempdir, remove_tempdir
from utils.sequtils import wrap
from vcf_to_fasta import vcf_to_fasta

from utils.sequtils import wrap, extract_amplicons
from utils.alignutils import TilingRow, align_nucmer, show_aligns, overlapper

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
    group1.add_argument('--con_fa', type=existing_file, required=True,
                        help='Fasta file with consensus sequence')
    group1.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Fasta file with reference sequence')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('Fix consensus options')
    # group2.add_argument('--bt2_preset', 
    #                     choices=['very-fast', 'fast', 'sensitive', 'very-sensitive',
    #                              'very-fast-local', 'fast-local', 'sensitive-local',
    #                              'very-sensitive-local',],
    #                     help='Bowtie2 preset to use')
    # group2.add_argument('--min_dp', type=int,
    #                     help='Minimum depth to call position')
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

def fix_consensus(fq1=None, fq2=None, fqU=None, con_fa=None, ref_fa=None, outdir='.',
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
    
    # Outputs
    out_con = os.path.join(outdir, 'consensus.fasta')
    out_bam = os.path.join(outdir, 'aligned.bam')
    out_vcf = os.path.join(outdir, 'variants.vcf')    
    
    # Temporary directory
    tempdir = create_tempdir('fix_consensus')
    
    # Align consensus to reference
    check_dependency('nucmer')
    check_dependency('delta-filter')
    check_dependency('show-tiling')
    
    print >>sys.stderr, 'Align consensus to reference'
    asm_dict = {s.id:s for s in SeqIO.parse(con_fa, 'fasta')}
    ref_dict = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    
    for chrom in sorted(asm_dict.keys()):
        # Extract amplicons, allow for longer strings of "n"    
        amps = extract_amplicons(chrom, str(asm_dict[chrom].seq), 200)
        tmp_amplicons_fa = os.path.join(tempdir, '%s.amplicons.fasta' % chrom)
        with open(tmp_amplicons_fa, 'w') as outh:
            for i, amp in enumerate(amps):
                print >>outh, '>%s.%d' % (amp[0].split()[0], i)
                print >>outh, '%s' % amp[1]

        # Align amplicons
        fil, til = align_nucmer(tmp_amplicons_fa, ref_fa, tempdir)
        
        # Parse tiling
        tr_list = [TilingRow(l) for l in open(til, 'rU')]

        # Get references
        refs = list(set([tr.ref for tr in tr_list]))
        if not len(refs) == 1:
            raise PipelineStepError("Amplicons align to multiple references")

        ref = refs[0]
        ranked = sorted(tr_list, key=lambda x:x.qry_alen, reverse=True)
        for tr in tr_list:
            tr.ref_align(fil)
        
    if not keep_tmp:
        remove_tempdir(tempdir, 'fix_consensus')

    return out_con, out_bam, out_vcf    
"""    
    
    # Copy and index initial reference
    curref = os.path.join(tempdir, 'initial.fasta')
    cmd1 = ['cp', ref_fa, curref]
    cmd2 = ['samtools', 'faidx', curref]
    cmd3 = ['picard', 'CreateSequenceDictionary', 
            'R=%s' % curref, 'O=%s' % os.path.join(tempdir, 'initial.dict')]
    cmd4 = ['bowtie2-build', curref, os.path.join(tempdir, 'initial')]
    command_runner([cmd1,cmd2,cmd3,cmd4], 'refine_assembly:index_ref', debug)
    
    # Align with bowtie2
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
    if fq1 is not None:
        cmd5 += ['-1', fq1, '-2', fq2,]
    if fqU is not None:
        cmd5 += ['-U', fqU, ]
    cmd5 += ['|', 'samtools', 'view', '-bS', '-', '>', os.path.join(tempdir, 'unsorted.bam'),]
    cmd6 = ['samtools', 'sort', os.path.join(tempdir, 'unsorted.bam'), os.path.join(tempdir, 'aligned'),]
    cmd7 = ['samtools', 'index',  os.path.join(tempdir, 'aligned.bam'),]
    command_runner([cmd5,cmd6,cmd7,], 'refine_assembly:align', debug)
    
    # MarkDuplicates
    cmd8 = [
        'picard', 'MarkDuplicates',
        'REMOVE_DUPLICATES=true',
        'CREATE_INDEX=true',
        'M=%s' % os.path.join(tempdir, 'rmdup.metrics.txt'),
        'I=%s' % os.path.join(tempdir, 'aligned.bam'),
        'O=%s' % os.path.join(tempdir, 'rmdup.bam'),
    ]
    # RealignerTargetCreator
    cmd9 = [
        'gatk', '-T', 'RealignerTargetCreator',
        '-I', os.path.join(tempdir, 'rmdup.bam'),
        '-R', curref,
        '-o', os.path.join(tempdir, 'tmp.intervals'),
    ]
    # IndelRealigner
    cmd10 = [
        'gatk', '-T', 'IndelRealigner',
        '-maxReads', '1000000',
        '-dt', 'NONE',
        '-I', os.path.join(tempdir, 'rmdup.bam'),
        '-R', curref,
        '-targetIntervals', os.path.join(tempdir, 'tmp.intervals'),
        '-o', os.path.join(tempdir, 'realign.bam')
    ]
    # UnifiedGenotyper
    cmd11 = [
        'gatk', '-T', 'UnifiedGenotyper',
        '--num_threads', '%d' % ncpu,
        '-out_mode', 'EMIT_ALL_SITES',
        '-glm', 'BOTH',
        '--baq', 'OFF',
        '--useOriginalQualities',
        '-dt', 'NONE',
        # '-stand_call_conf', '0',
        # '-stand_emit_conf', '0',
        '-A', 'AlleleBalance',
        '--min_base_quality_score', '15',
        '-ploidy', '4',
        '-I', os.path.join(tempdir, 'realign.bam'),
        '-R', curref,
        '-o', os.path.join(tempdir, 'tmp.vcf.gz'),
    ]
    command_runner([cmd8,cmd9,cmd10,cmd11,], 'refine_assembly:consensus', debug)
    
    # Call consensus from VCF
    if os.path.exists(os.path.join(tempdir, 'tmp.vcf.gz')):
        with open(out1, 'w') as outh:
            vcf_to_fasta(os.path.join(tempdir, 'tmp.vcf.gz'), sampidx=0, min_dp=min_dp, outfile=outh)
    else:
        print >>sys.stderr, 'VCF file %s not found' % os.path.join(tempdir, 'tmp.vcf.gz')

    if not keep_tmp:
        remove_tempdir(tempdir, 'assemble_scaffold')
    
    return out1
"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix consensus sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))

    



