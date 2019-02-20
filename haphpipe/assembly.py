#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from haphpipe.utils.sysutils import args_params

from haphpipe.stages import assemble_denovo
from haphpipe.stages import assemble_scaffold
from haphpipe.stages import assemble_amplicons
from haphpipe.stages import align_reads
from haphpipe.stages import call_variants
from haphpipe.stages import refine_assembly
from haphpipe.stages import fix_consensus
from haphpipe.stages import pairwise_align
from haphpipe.stages import vcf_to_fasta
from haphpipe.stages import post_assembly
from haphpipe.stages import extract_pairwise
from haphpipe.stages import annotate_from_ref

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"

def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(help='Assembly stages')
    assemble_denovo.stageparser(sub.add_parser('assemble_denovo'))
    # assign_contigs.stageparser(sub.add_parser('assign_contigs'))
    assemble_scaffold.stageparser(sub.add_parser('assemble_scaffold'))
    # impute_ref.stageparser(sub.add_parser('impute_ref'))
    assemble_amplicons.stageparser(sub.add_parser('assemble_amplicons'))
    refine_assembly.stageparser(sub.add_parser('refine_assembly'))
    align_reads.stageparser(sub.add_parser('align_reads'))
    call_variants.stageparser(sub.add_parser('call_variants'))
    pairwise_align.stageparser(sub.add_parser('pairwise_align'))    
    fix_consensus.stageparser(sub.add_parser('fix_consensus'))
    vcf_to_fasta.stageparser(sub.add_parser('vcf_to_fasta'))
    
    post_assembly.stageparser(sub.add_parser('post_assembly'))
    extract_pairwise.stageparser(sub.add_parser('extract_pairwise'))
    annotate_from_ref.stageparser(sub.add_parser('annotate_from_ref'))
    
    args = parser.parse_args()
    args.func(**args_params(args))

if __name__ == '__main__':
    main()
