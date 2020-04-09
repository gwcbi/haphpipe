#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from haphpipe.utils import sysutils

from haphpipe.stages import demo

__author__ = 'Keylie M. Gibson'
__copyright__ = "Copyright (C) 2020 Keylie M. Gibson"


def main():
    parser = argparse.ArgumentParser(
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    sub = parser.add_subparsers(
        help='''setup demo directory and test data.
             '''
    )
    demo.stageparser(sub.add_parser('demo'))

    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    main()
