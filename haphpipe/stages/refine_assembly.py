#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import absolute_import
from builtins import input

import os
import shutil
import re
import argparse
import random

from Bio import SeqIO
from Bio import pairwise2

from haphpipe.utils import sysutils
from haphpipe.stages import align_reads
from haphpipe.stages import call_variants
from haphpipe.stages import vcf_to_consensus
from haphpipe.stages import sample_reads


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
    group1.add_argument('--ref_fa', type=sysutils.existing_file, required=True,
                        help='Assembly to refine')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('Refinement options')
    group2.add_argument('--max_step', type=int, default=1,
                        help='Maximum number of refinement steps')
    group2.add_argument('--subsample', type=int,
                        help='Use a subsample of reads for refinement.')
    group2.add_argument('--seed', type=int,
                        help='''Seed for random number generator (ignored if
                                not subsampling).''')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int, default=1,
                        help='Number of CPUs to use')
    group3.add_argument('--xmx', type=int,
                        help='Maximum heap size for Java VM, in gigabytes.')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=refine_assembly)


def refine_assembly(**kwargs):
    """

    Args:
        **kwargs:

    Returns:

    """
    #     fq1=None, fq2=None, fqU=None, ref_fa=None, outdir='.',
    #     max_step=1, subsample=None,
    #     ncpu=1, xmx=None, keep_tmp=False,
    #     quiet=False, logfile=None, debug=False,
    #):
    max_step = kwargs.pop('max_step')
    if max_step == 1:
        kwargs['iteration'] = None
        return refine_assembly_step(**kwargs)
        #     fq1=fq1, fq2=fq2, fqU=fqU, ref_fa=ref_fa, outdir=outdir,
        #     iteration=None, subsample=subsample,
        #    ncpu=ncpu, xmx=xmx, keep_tmp=keep_tmp,
        #    quiet=quiet, logfile=logfile, debug=debug
        #)
    else:
        kwargs['max_step'] = max_step
        return progressive_refine_assembly(**kwargs)
        #    fq1=fq1, fq2=fq2, fqU=fqU, ref_fa=ref_fa, outdir=outdir,
        #    max_step=max_step, subsample=subsample,
        #    ncpu=ncpu, keep_tmp=keep_tmp,
        #    quiet=quiet, logfile=logfile, debug=debug
        #)


def refine_assembly_step(
        fq1=None, fq2=None, fqU=None, ref_fa=None, outdir='.',
        iteration=None, subsample=None, seed=None,
        ncpu=1, keep_tmp=False,
        quiet=False, logfile=None, debug=False,
    ):
    # Temporary directory
    tempdir = sysutils.create_tempdir('refine_assembly', None, quiet, logfile)

    if subsample is None:
        sq1, sq2, sqU = fq1, fq2, fqU
    else:
        seed = seed if seed is not None else random.randrange(1, 1000)
        sq1, sq2, sqU = sample_reads.sample_reads(
            fq1=fq1, fq2=fq2, fqU=fqU, outdir=tempdir,
            nreads=subsample, seed=seed,
            quiet=False, logfile=logfile, debug=debug
        )


    tmp_aligned, tmp_bt2 = align_reads.align_reads(
        fq1=sq1, fq2=sq2, fqU=sqU, ref_fa=ref_fa, outdir=tempdir,
        ncpu=ncpu, keep_tmp=keep_tmp,
        quiet=quiet, logfile=logfile, debug=debug,
    )

    tmp_vcf = call_variants.call_variants(
        aln_bam=tmp_aligned, ref_fa=ref_fa, outdir=tempdir,
        emit_all=True,
        ncpu=ncpu, keep_tmp=keep_tmp,
        quiet=quiet, logfile=logfile, debug=debug,
    )

    tmp_fasta = vcf_to_consensus.vcf_to_consensus(
        vcf=tmp_vcf, outdir=tempdir, sampidx=0,
        keep_tmp=keep_tmp,
        quiet=quiet, logfile=logfile
    )

    # Copy command
    if iteration is None:
        out_refined = os.path.join(outdir, 'refined.fasta')
        out_bt2 = os.path.join(outdir, 'refined_bt2.out')
    else:
        out_refined = os.path.join(outdir, 'refined.%02d.fasta' % iteration)
        out_bt2 = os.path.join(outdir, 'refined_bt2.%02d.out' % iteration)

    shutil.copy(tmp_fasta, out_refined)
    shutil.copy(tmp_bt2, out_bt2)

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'refine_assembly', quiet, logfile)

    return out_refined, out_bt2


def progressive_refine_assembly(
        fq1=None, fq2=None, fqU=None, ref_fa=None, outdir='.',
        max_step=None, subsample=None, seed=None,
        ncpu=1, keep_tmp=False,
        quiet=False, logfile=None, debug=False,
    ):

    # Outputs
    out_refined = os.path.join(outdir, 'refined.fasta')
    out_bt2 = os.path.join(outdir, 'refined_bt2.out')
    out_summary = os.path.join(outdir, 'refined_summary.out')

    # Initialize
    seq_ids = []
    cur_asm = ref_fa
    cur_alnrate = None
    cur_seqs = {}

    for s in SeqIO.parse(cur_asm, 'fasta'):
        seq_ids.append(s.id)
        cur_seqs[s.id] = s

    # Message log for summary
    summary = [
        ['iteration', 'alnrate', 'diffs'] + ['diff:%s' % s for s in seq_ids]
    ]

    # Seed random number generator
    random.seed(seed)

    for i in range(1, max_step+1):
        # Generate a refined assembly
        tmp_refined, tmp_bt2 = refine_assembly_step(
            fq1=fq1, fq2=fq2, fqU=fqU, ref_fa=cur_asm, outdir=outdir,
            iteration=i, subsample=subsample,
            ncpu=ncpu, keep_tmp=keep_tmp,
            quiet=True, logfile=logfile, debug=debug
        )

        # Check whether alignments are different
        new_seqs = {s.id: s for s in SeqIO.parse(tmp_refined, 'fasta')}
        diffs = []
        for sid in seq_ids:
            alns = pairwise2.align.globalms(
                new_seqs[sid].seq, cur_seqs[sid].seq, 2, -1, -3, -1
            )
            d = min(sum(nc != cc for nc, cc in zip(t[0], t[1])) for t in alns)
            diffs.append(d)

        # Check new alignment rate
        with open(tmp_bt2, 'rU') as fh:
            m = re.search('(\d+\.\d+)\% overall alignment rate', fh.read())
            new_alnrate = float(m.group(1))

        # Create messages for log
        row = [str(i), '%.02f' % new_alnrate, '%d' % sum(diffs), ]
        row += list(map(str, diffs))
        summary.append(row)

        # Create messages for console
        sysutils.log_message('\nRefinement result:\n', quiet, logfile)
        sysutils.log_message('\tDifferences:\n', quiet, logfile)
        for s,d in zip(seq_ids, diffs):
            sysutils.log_message('\t\t%s\t%d\n' % (s,d), quiet, logfile)
        if sum(diffs) > 0:
            msg = '\t%d differences found with previous\n' % sum(diffs)
        else:
            msg = '\tNo differences with previous\n'
        sysutils.log_message(msg, quiet, logfile)

        if cur_alnrate is None:
            msg = '\tAlignment rate: %0.2f\n' % new_alnrate
        elif new_alnrate > cur_alnrate:
            msg = '\tAlignment rate has improved: '
            msg += '%.02f > %.02f\n' % (new_alnrate, cur_alnrate)
        else:
            msg = '\tAlignment rate has not improved: '
            msg += '%.02f <= %.02f\n' % (new_alnrate, cur_alnrate)
        sysutils.log_message(msg, quiet, logfile)

        # Decide whether to keep going
        keep_going = True
        if sum(diffs) == 0:
            keep_going = False
            msg = 'Stopping because no differences found\n'
            sysutils.log_message(msg, quiet, logfile)

        # We should also quit if alignment rate does not improve
        # However, subsampling reads can lead to changes in alignment rate
        # that can be ignore. When subsampling is implemented the first
        # boolean value should check whether subsampling is enabled
        if subsample is None: # not subsampling
            if cur_alnrate is not None and new_alnrate <= cur_alnrate:
                keep_going = False
                msg = 'Stopping because alignment rate did not improve\n'
                sysutils.log_message(msg, quiet, logfile)

        cur_asm = tmp_refined
        cur_alnrate = new_alnrate
        cur_seqs = new_seqs

        if not keep_going:
            break

    # Final outputs
    shutil.copy(cur_asm, out_refined)
    shutil.copy(tmp_bt2, out_bt2)

    with open(out_summary, 'w') as outh:
        print('\n'.join('\t'.join(r) for r in summary), file=outh)

    return out_refined, out_bt2, out_summary


'''
def refine_assembly(fq1=None, fq2=None, fqU=None, assembly_fa=None, outdir='.',
        bt2_preset='sensitive', min_dp=1, rgid='test',
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
    
    if encoding is None:
        encoding = guess_encoding(fq1) if fq1 is not None else guess_encoding(fqU)
    
    # Check dependencies
    check_dependency('bowtie2')
    check_dependency('samtools')
    check_dependency('picard')
    check_dependency('gatk')
    
    # Outputs
    out_refined = os.path.join(outdir, 'refined.fa')
    out_vcf = os.path.join(outdir, 'variants.vcf.gz')    
    
    # Temporary directory
    tempdir = create_tempdir('refine_assembly')
    
    # Copy and index initial reference
    curref = os.path.join(tempdir, 'initial.fasta')
    cmd1 = ['cp', assembly_fa, curref]
    cmd2 = ['samtools', 'faidx', curref]
    cmd3 = ['picard', 'CreateSequenceDictionary', 
            'R=%s' % curref, 'O=%s' % os.path.join(tempdir, 'initial.dict')]
    cmd4 = ['bowtie2-build', curref, os.path.join(tempdir, 'initial')]
    command_runner([cmd1,cmd2,cmd3,cmd4], 'refine_assembly:index_ref', debug)
    
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
    cmd5 += ['-S', os.path.join(tempdir, 'tmp.sam'), ]
    cmd5 += ['2>&1', '|', 'tee', out_bt2, ]
    cmd6a = ['samtools', 'view', '-bS', os.path.join(tempdir, 'tmp.sam'), '>', os.path.join(tempdir, 'unsorted.bam'),]
    cmd6b = ['samtools', 'sort', os.path.join(tempdir, 'unsorted.bam'), os.path.join(tempdir, 'aligned'),]
    cmd6c = ['samtools', 'index',  os.path.join(tempdir, 'aligned.bam'),]
    cmd7 = ['rm', '-f', os.path.join(tempdir, 'tmp.sam'), os.path.join(tempdir, 'unsorted.bam'), ]
    command_runner([cmd5,cmd6a,cmd6b,cmd6c,cmd7], 'refine_assembly:align', debug)
    
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
        '-o', out_vcf,
    ]
    command_runner([cmd8,cmd9,cmd10,cmd11,], 'refine_assembly:consensus', debug)
    
    # Call consensus from VCF
    if os.path.exists(out_vcf):
        with open(out_refined, 'w') as outh:
            vcf_to_fasta(out_vcf, sampidx=0, min_dp=min_dp, outfile=outh)
    else:
        print('VCF file %s not found' % out_vcf, file=sys.stderr)

    if not keep_tmp:
        remove_tempdir(tempdir, 'assemble_scaffold')
    
    return out_refined
'''


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='''Three step assembly refinement: align reads, call
                       variants, and update reference''',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
