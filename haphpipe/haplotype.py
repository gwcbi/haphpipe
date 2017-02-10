#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from utils.sysutils import args_params

from stages import predict_haplo

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(help='Assembly stages')
    predict_haplo.stageparser(sub.add_parser('predict_haplo'))
        
    args = parser.parse_args()
    args.func(**args_params(args))

if __name__ == '__main__':
    main()
