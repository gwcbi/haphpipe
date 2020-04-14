#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from haphpipe.utils import sysutils
from haphpipe.stages import multiple_align
from haphpipe.stages import model_test
from haphpipe.stages import build_tree


__author__ = 'Keylie M. Gibson'
__copyright__ = "Copyright (C) 2020 Keylie M. Gibson"

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

def console():
    parser = argparse.ArgumentParser(
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    sub = parser.add_subparsers(
        help='''Phylogenomic stages. Includes multiple alignment, evolutionary 
                model testing using ModelTest-ng, and building phylogenetic trees 
                with RAxML.
                '''
    )
    multiple_align.stageparser(sub.add_parser('multiple_align'))
    model_test.stageparser(sub.add_parser('model_test'))
    build_tree.stageparser(sub.add_parser('build_tree'))


    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
