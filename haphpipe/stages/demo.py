#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import os
import argparse

from haphpipe.utils import helpers
from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import PipelineStepError
from haphpipe.utils.sysutils import MissingRequiredArgument


__author__ = 'Matthew L. Bendall and Keylie M. Gibson'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall and 2020 Keylie M. Gibson"


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
    parser.set_defaults(func=demo)

def demo(
        outdir="."
    ):
    try:
        _ = FileNotFoundError()
    except NameError:
        class FileNotFoundError(OSError):
            pass

    # This file, demo.py, is located within "stages", so the package root is
    # up one directory
    _base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    _data = os.path.join(_base, 'data')

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    #dest = os.path.abspath(outdir)
    #if not os.path.exists(os.path.join(outdir,))


    print(_base, file=sys.stderr)

    # Check for executable
    sysutils.check_dependency("fastq-dump")

    # Demo command
    # fix once everything is uploaded
    cmd1 = [
        'haphpipe_demo', 'haphpipe_demo'
    ]

    sysutils.command_runner(
        [cmd1, ], 'demo'
    )

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
