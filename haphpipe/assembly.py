#! /usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

from utils.sysutils import args_params

import stages
from stages import trim_reads
from stages import join_reads
from stages import ec_reads
from stages import assemble_denovo
# from stages import assign_contigs
from stages import assemble_scaffold
# from stages import impute_ref
from stages import refine_assembly
from stages import fix_consensus
from stages import vcf_to_fasta

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(help='Assembly stages')
    trim_reads.stageparser(sub.add_parser('trim_reads'))
    join_reads.stageparser(sub.add_parser('join_reads'))    
    ec_reads.stageparser(sub.add_parser('ec_reads'))
    assemble_denovo.stageparser(sub.add_parser('assemble_denovo'))
    # assign_contigs.stageparser(sub.add_parser('assign_contigs'))
    assemble_scaffold.stageparser(sub.add_parser('assemble_scaffold'))
    # impute_ref.stageparser(sub.add_parser('impute_ref'))
    refine_assembly.stageparser(sub.add_parser('refine_assembly'))
    fix_consensus.stageparser(sub.add_parser('fix_consensus'))
    vcf_to_fasta.stageparser(sub.add_parser('vcf_to_fasta'))
    
    args = parser.parse_args()
    args.func(**args_params(args))

if __name__ == '__main__':
    main()
