#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from haphpipe.utils import sysutils

from haphpipe.stages import pairwise_align
from haphpipe.stages import post_assembly
from haphpipe.stages import extract_pairwise
from haphpipe.stages import annotate_from_ref
from haphpipe.stages import summary_stats



__author__ = 'Matthew L. Bendall and Keylie M. Gibson'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall; (C) 2020 Keylie M. Gibson"


def main():
    parser = argparse.ArgumentParser(
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    sub = parser.add_subparsers(
        help='''Annotate consensus sequence(s). An assembled consensus sequence
                can be amino acid aligned to a reference (i.e. HXB2) so the
                same coordinate system can be used. Sequences can then be
                extracted based on the shared coordinate systems.
             '''
    )
    pairwise_align.stageparser(sub.add_parser('pairwise_align'))
    extract_pairwise.stageparser(sub.add_parser('extract_pairwise'))
    post_assembly.stageparser(sub.add_parser('post_assembly'))
    annotate_from_ref.stageparser(sub.add_parser('annotate_from_ref'))
    summary_stats.stageparser(sub.add_parser('summary_stats'))

    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    main()
