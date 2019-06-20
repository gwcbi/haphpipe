# -*- coding: utf-8 -*-
from __future__ import print_function

import os
import argparse

import yaml

from haphpipe.utils import sysutils


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
    group1.add_argument('--outdir', type=sysutils.existing_dir,
                        help='Output directory')
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
    parser.set_defaults(func=ec_reads)


def ec_reads(
        fq1=None, fq2=None, fqU=None, outdir='.',
        ncpu=1, keep_tmp=False, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to error-correct reads using spades

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        fqU (str): Path to fastq file with unpaired reads
        outdir (str): Path to output directory
        ncpu (int): Number of CPUs to use
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out1 (str): Path to corrected fastq file with read 1
        out2 (str): Path to corrected fastq file with read 2
        outU (str): Path to corrected fastq file with unpaired reads

    """
    # Check inputs
    if fq1 is not None and fq2 is not None and fqU is None:
        input_reads = "paired"  # Paired end
    elif fq1 is None and fq2 is None and fqU is not None:
        input_reads = "single"  # Single end
    elif fq1 is not None and fq2 is not None and fqU is not None:
        input_reads = "both"
    else:
        msg = "Incorrect combination of reads: "
        msg += "fq1=%s fq2=%s fqU=%s" % (fq1, fq2, fqU)
        raise sysutils.PipelineStepError(msg)

    # Check dependencies
    sysutils.check_dependency('spades.py')

    # Outputs
    out1 = os.path.join(outdir, 'corrected_1.fastq')
    out2 = os.path.join(outdir, 'corrected_2.fastq')
    outU = os.path.join(outdir, 'corrected_U.fastq')

    # Temporary directory
    tempdir = sysutils.create_tempdir('ec_reads', None, quiet, logfile)

    # spades command
    cmd1 = [
        'spades.py',
        '-o', tempdir,
        '-t', '%d' % ncpu,
        '--only-error-correction',
    ]
    if input_reads in ['paired', 'both', ]:
        cmd1 += ['-1', os.path.abspath(fq1),
                 '-2', os.path.abspath(fq2),
                 ]
    if input_reads in ['single', 'both', ]:
        cmd1 += ['-s', os.path.abspath(fqU), ]

    sysutils.command_runner([cmd1,], 'ec_reads', quiet, logfile, debug)

    # Copy files
    yaml_file = os.path.join(tempdir, 'corrected/corrected.yaml')
    if not os.path.exists(yaml_file):
        sysutils.PipelineStepError("YAML file %s not found" % yaml_file)

    with open(yaml_file, 'rU') as fh:
        d = yaml.load(fh, Loader=yaml.FullLoader)[0]
    cmds = []
    if 'left reads' in d:
        cmds.append(
            ['gunzip', '-c', ] + sorted(d['left reads']) + ['>', out1]
        )
    if 'right reads' in d:
        cmds.append(
            ['gunzip', '-c', ] + sorted(d['right reads']) + ['>', out2]
        )
    if 'single reads' in d:
        cmds.append(
            ['gunzip', '-c', ] + sorted(d['single reads']) + ['>', outU]
        )

    sysutils.command_runner(cmds, 'ec_reads', quiet, logfile, debug)

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'ec_reads', quiet, logfile)

    return out1, out2, outU


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Error correct reads using spades',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
