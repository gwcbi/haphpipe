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

def main():
    parser = argparse.ArgumentParser(
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    sub = parser.add_subparsers(
        help='''Annotate consensus sequence(s). An assembled consensus sequence
                can be amino acid aligned to a reference (i.e. HXB2) so the
                same coordinate system can be used. Sequences can then be
                extracted based on the shared coordinate systems. Summary 
                statistics can also be calculated.
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
