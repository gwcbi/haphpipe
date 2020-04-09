#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from haphpipe.utils import sysutils
from haphpipe.stages import multiple_align
from haphpipe.stages import model_test
#from haphpipe.stages import build_tree


__author__ = 'Keylie M. Gibson'
__copyright__ = "Copyright (C) 2020 Keylie M. Gibson"


def console():
    parser = argparse.ArgumentParser(
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    sub = parser.add_subparsers(
        help='''Insert descirption about phylo'''
    )
    multiple_align.stageparser(sub.add_parser('multiple_align'))
    model_test.stageparser(sub.add_parser('model_test'))
    #build_tree.stageparser(sub.add_parser('build_tree'))


    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
