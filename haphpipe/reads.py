#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from haphpipe.utils import sysutils
from haphpipe.stages import sample_reads
from haphpipe.stages import trim_reads
from haphpipe.stages import join_reads
from haphpipe.stages import ec_reads


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

def console():
    parser = argparse.ArgumentParser(
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    sub = parser.add_subparsers(
        help='''Manipulate reads. Input is reads in fastq format, output is
                modified reads in fastq format.'''
    )
    sample_reads.stageparser(sub.add_parser('sample_reads'))
    trim_reads.stageparser(sub.add_parser('trim_reads'))
    join_reads.stageparser(sub.add_parser('join_reads'))
    ec_reads.stageparser(sub.add_parser('ec_reads'))

    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
