#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse

from haphpipe.utils import helpers
from haphpipe.utils import sysutils


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

def stageparser(parser):
    """ Add stage-specific options to argparse parser

    Args:
        parser (argparse.ArgumentParser): ArgumentParser object

    Returns:
        None

    """
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--fq1', type=sysutils.existing_file, required=True,
                        help='Fastq file with read 1')
    group1.add_argument('--fq2', type=sysutils.existing_file, required=True,
                        help='Fastq file with read 1')
    group1.add_argument('--outdir', type=sysutils.existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('FLAsh settings')
    group2.add_argument('--min_overlap', type=int, default=10,
                        help='''The minimum required overlap length between two
                                reads to provide a confident overlap.''')
    group2.add_argument('--max_overlap', type=int,
                        help='''Maximum overlap length expected in
                                approximately 90 pct of read pairs, longer
                                overlaps are penalized.''')
    group2.add_argument('--allow_outies', action='store_true',
                        help='''Also try combining read pairs in the "outie"
                                orientation''')
    group2.add_argument('--encoding',
                        help='Quality score encoding')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
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
    parser.set_defaults(func=join_reads)

def join_reads(
        fq1=None, fq2=None, outdir=".",
        min_overlap=None, max_overlap=None, allow_outies=None,
        encoding=None,
        ncpu=1, keep_tmp=False, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to join paired-end reads

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        outdir (str): Path to output directory
        min_overlap (int): The minimum required overlap length
        max_overlap (int): Maximum overlap length
        allow_outies (bool): Try combining "outie" reads
        encoding (str): Quality score encoding
        ncpu (int): Number of CPUs to use
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out1 (str): Path to fastq file with unjoined read 1
        out2 (str): Path to fastq file with unjoined read 2
        outU (str): Path to fastq file with joined reads

    """
    # Check inputs
    if fq1 is not None and fq2 is not None:
        pass # Both are present
    else:
        msg = "Incorrect combination of reads: fq1=%s fq2=%s" % (fq1, fq2)
        raise sysutils.PipelineStepError(msg)
    
    # Check for executable
    sysutils.check_dependency('flash')

    # Get encoding
    if encoding is None:
        encoding = helpers.guess_encoding(fq1)
    
    # Outputs
    outU = os.path.join(outdir, 'joined.fastq')
    out1 = os.path.join(outdir, 'notjoined_1.fastq')
    out2 = os.path.join(outdir, 'notjoined_2.fastq')

    # Temporary directory
    tempdir = sysutils.create_tempdir('join_reads', None, quiet, logfile)

    # Flash command
    cmd1 = [
        'flash',
        '-t', '%d' % ncpu,
        '-d', tempdir,
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

    cmd2 = ['mv', os.path.join(tempdir, 'out.extendedFrags.fastq'), outU, ]
    cmd3 = ['mv', os.path.join(tempdir, 'out.notCombined_1.fastq'), out1, ]
    cmd4 = ['mv', os.path.join(tempdir, 'out.notCombined_2.fastq'), out2, ]
    sysutils.command_runner(
        [cmd1, cmd2, cmd3, cmd4, ], 'join_reads', quiet, logfile, debug
    )

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'join_reads', quiet, logfile)

    return out1, out2, outU


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Join reads using FLASh.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()

