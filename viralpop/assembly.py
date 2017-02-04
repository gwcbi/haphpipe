#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from utils.sysutils import args_params

from stages import assemble_scaffold

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(help='Assembly stages')    
    assemble_scaffold.stageparser(sub.add_parser('assemble_scaffold'))
    
    args = parser.parse_args()
    args.func(**args_params(args))
