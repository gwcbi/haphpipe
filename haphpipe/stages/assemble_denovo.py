#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import shutil

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils
from haphpipe.stages import sample_reads
from haphpipe.utils.sysutils import MissingRequiredArgument


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
    
    group2 = parser.add_argument_group('Assembly options')
    try:
        sysutils.check_dependency('Trinity')
        is_trinity = True
    except sysutils.PipelineStepError:
        is_trinity = False
    try:
        sysutils.check_dependency('spades.py')
        is_spades = True
    except sysutils.PipelineStepError:
        is_spades = False

    if is_trinity and is_spades:
        group2.add_argument('--assembler', default='spades',
                            choices=['spades', 'trinity', ],
                            help='''Assembler to use.''')
    elif is_trinity:
        group2.set_defaults(assembler="trinity")
    elif is_spades:
        group2.set_defaults(assembler="spades")

    if is_spades:
        group2.add_argument('--no_error_correction', action='store_true',
                            help='Do not perform error correction [spades only]')
    if is_trinity:
        group2.add_argument('--min_contig_length', type=int, default=200,
                            help='''Minimum assembled contig length to report
                                    [Trinity only]''')
    group2.add_argument('--subsample', type=int,
                        help='Use a subsample of reads for assembly.')
    group2.add_argument('--seed', type=int,
                        help='''Seed for random number generator (ignored if
                                not subsampling).''')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int, default=1,
                        help='Number of CPU to use')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Keep temporary directory')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=assemble_denovo)


def assemble_denovo(**kwargs):
    """ Wrapper function for denovo assembly

    Arguments are passed to specific assembly functions depending on the
    'assembler' argument

    Args:
        **kwargs: Arguments passed to specific assembly functions

    Returns:

    """
    if kwargs['assembler'] == 'spades':
        return assemble_denovo_spades(**kwargs)
    elif kwargs['assembler'] == 'trinity':
        return assemble_denovo_trinity(**kwargs)
    else:
        raise sysutils.PipelineStepError('INVALID ASSEMBLER.')

def assemble_denovo_spades(
        fq1=None, fq2=None, fqU=None, outdir='.',
        no_error_correction=False, subsample=None, seed=None,
        ncpu=1, keep_tmp=False, quiet=False, logfile=None, debug=False,
        **kwargs
    ):
    """ Pipeline step to assemble reads using spades (denovo)

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        fqU (str): Path to fastq file with unpaired reads
        outdir (str): Path to output directory
        no_error_correction (bool): do not perform error correction
        subsample (int): use a subsample of reads for assembly
        seed (int): Seed for random number generator
        ncpu (int): Number of CPUs to use
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run
        **kwargs: Not used.

    Returns:
        out_fa (str): Path to assembled contigs file (fasta format)
        out_summary (str): Path to assembly summary

    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired"  # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single"  # Single end
    elif fq1 is not None and fq2 is not None and fqU is not None:
        input_reads = "both"
    else:
        msg = "incorrect input reads; requires either "
        msg += "(--fq1 AND --fq2) OR (--fqU) OR (--fq1 AND --fq2 AND --fqU)"
        raise MissingRequiredArgument(msg)

    # Check dependencies
    sysutils.check_dependency('spades.py')

    # Outputs
    out_fa = os.path.join(outdir, 'denovo_contigs.fna')
    out_summary = os.path.join(outdir, 'denovo_summary.txt')

    # Temporary directory
    tempdir = sysutils.create_tempdir('assemble_spades', None, quiet, logfile)

    # Subsample
    if subsample is not None:
        full1, full2, fullU = fq1, fq2, fqU
        fq1, fq2, fqU = sample_reads.sample_reads(
                fq1=full1, fq2=full2, fqU=fullU, outdir=tempdir,
                nreads=subsample, seed=seed,
                quiet=quiet, logfile=logfile, debug=debug
            )

    # spades command
    cmd1 = [
        'spades.py',
        '-o', tempdir,
        '-t', '%d' % ncpu,
    ]
    if input_reads in ['paired', 'both', ]:
        cmd1 += ['-1', os.path.abspath(fq1),
                 '-2', os.path.abspath(fq2),
                 ]
    if input_reads in ['single', 'both', ]:
        cmd1 += ['-s', os.path.abspath(fqU), ]
    if no_error_correction:
        cmd1 += ['--only-assembler', ]

    sysutils.command_runner(
        [cmd1, ], 'assemble_spades', quiet, logfile, debug
    )
    shutil.copy(os.path.join(tempdir, 'contigs.fasta'), out_fa)

    if os.path.isfile(out_fa):
        with open(out_summary, 'w') as outh:
            sequtils.assembly_stats(open(out_fa, 'rU'), outh)

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'assemble_spades', quiet, logfile)

    return out_fa, out_summary


def assemble_denovo_trinity(
        fq1=None, fq2=None, fqU=None, outdir='.',
        min_contig_length=200, subsample=None, seed=None,
        ncpu=1, keep_tmp=False, quiet=False, logfile=None, debug=False,
        **kwargs
    ):
    """ Pipeline step to assemble reads using Trinity (denovo)

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        fqU (str): Path to fastq file with unpaired reads
        outdir (str): Path to output directory
        min_contig_length (int): minimum assembled contig length to report
        subsample (int): use a subsample of reads for assembly
        seed (int): Seed for random number generator
        ncpu (int): Number of CPUs to use
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run
        **kwargs: Not used.

    Returns:
        out1 (str): Path to assembled contigs file (fasta format)

    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired" # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single" # Single end
    elif fq1 is not None and fq2 is not None and fqU is not None:
        input_reads = "both"
    else:
        msg = "incorrect input reads; requires either "
        msg += "(--fq1 AND --fq2) OR (--fqU) OR (--fq1 AND --fq2 AND --fqU)"
        raise MissingRequiredArgument(msg)
    
    # Check dependencies
    sysutils.check_dependency('Trinity')
    
    # Outputs
    out1 = os.path.join(outdir, 'contigs.fa')

    # Temporary directory
    tempdir = sysutils.create_tempdir('assemble_trinity', None, quiet, logfile)
    
    # Trinity command
    cmd1 = [
        'Trinity',
        '--min_contig_length', '%d' % min_contig_length,
        '--CPU', '%d' % ncpu,
        #'--max_memory', '%dG' % max_memory,
        '--seqType', 'fq',
        '--output', tempdir,
    ]
    if input_reads in ['paired', 'both', ]:
        cmd1 += ['--left', os.path.abspath(fq1),
                 '--right', os.path.abspath(fq2),
                 ]
    elif input_reads in ['single', 'both', ]:
        cmd1 += ['--single', os.path.abspath(fqU), ]
    
    # Copy command
    cmd2 = ['cp',
        os.path.join(tempdir, 'Trinity.fasta'),
        out1,
    ]

    sysutils.command_runner(
        [cmd1, cmd2, ], 'assemble_trinity', quiet, logfile, debug
    )

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'assemble_trinity', quiet, logfile)
    
    if os.path.isfile(out1):
        with open(os.path.join(outdir, 'assembly_summary.txt'), 'w') as outh:
            sequtils.assembly_stats(open(out1, 'rU'), outh)
    
    return out1


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Assemble reads denovo.',
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
