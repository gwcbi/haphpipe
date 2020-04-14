#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from haphpipe.utils import sysutils

from haphpipe.stages import predict_haplo

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

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(help='Assembly stages')
    predict_haplo.stageparser(sub.add_parser('predict_haplo'))
        
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))

if __name__ == '__main__':
    main()
