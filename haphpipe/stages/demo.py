#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import os
import argparse
import shutil

from haphpipe.utils import helpers
from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import PipelineStepError
from haphpipe.utils.sysutils import MissingRequiredArgument


__author__ = 'Matthew L. Bendall and Keylie M. Gibson'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall and 2020 Keylie M. Gibson"

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
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    group1.add_argument('--refonly', action='store_true',
                        help='Does not run entire demo, only pulls the reference files')
    parser.set_defaults(func=demo)

def demo(outdir=".", refonly=False):
    try:
        _ = FileNotFoundError()
    except NameError:
        class FileNotFoundError(OSError):
            pass

    # This file, demo.py, is located within "stages", so the package root is
    # up one directory
    _base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    _data = os.path.abspath(os.path.join(_base,'refs'))
    #_data = os.path.abspath(os.path.join(os.path.dirname(_base), 'bin/refs'))
    #print(_data)
    #return
    #_data = os.path.join(_base, 'data')

    # get paths for reference files
    def get_data(path):
        return os.path.join(_data, path)

    data_amp_fasta = get_data('HIV_B.K03455.HXB2.amplicons.fasta')
    data_fasta = get_data('HIV_B.K03455.HXB2.fasta')
    data_gtf = get_data('HIV_B.K03455.HXB2.gtf')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    refs = os.path.join(outdir, 'haphpipe_demo/refs')
    if not os.path.exists(refs):
        os.makedirs(refs)

    # Copy reference files to directory
    shutil.copy(data_amp_fasta, refs)
    shutil.copy(data_fasta, refs)
    shutil.copy(data_gtf, refs)

    #dest = os.path.abspath(outdir)
    #if not os.path.exists(os.path.join(outdir,))

    print(_base, file=sys.stderr)

    if refonly is False:
        print("Setting up demo directories and references in outdirectory %s. Demo samples will now run." % os.path.join(outdir, 'haphpipe_demo'))

        # Check for executable
        sysutils.check_dependency("fastq-dump")

        # Demo command
        cmd1 = [
            'haphpipe_demo', 'haphpipe_demo'
        ]

        sysutils.command_runner(
            [cmd1, ], 'demo'
        )
    else:
        print("Demo was run with --refonly. References are now in outdirectory: %s." %refs)

def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Set up demo directory and run demo.',
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
