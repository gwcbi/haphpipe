#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from .utils.sysutils import args_params

from .stages import predict_haplo

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(help='Assembly stages')
    predict_haplo.stageparser(sub.add_parser('predict_haplo'))
        
    args = parser.parse_args()
    args.func(**args_params(args))

if __name__ == '__main__':
    main()
