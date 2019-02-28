#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import argparse

from haphpipe.utils import sysutils

from haphpipe.stages import assemble_denovo
from haphpipe.stages import assemble_amplicons
from haphpipe.stages import assemble_scaffold

from haphpipe.stages import align_reads
from haphpipe.stages import call_variants
from haphpipe.stages import vcf_to_consensus
from haphpipe.stages import refine_assembly
from haphpipe.stages import finalize_assembly


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


def main():
    parser = argparse.ArgumentParser(
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    sub = parser.add_subparsers(
        help='''Assemble consensus sequence(s). Input reads (in fastq format)
                are assembled using either denovo assembly or reference-based
                alignment. Resulting consensus can be further refined or
                annotated.
                '''
    )
    # Denovo
    assemble_denovo.stageparser(sub.add_parser('assemble_denovo'))
    assemble_amplicons.stageparser(sub.add_parser('assemble_amplicons'))
    assemble_scaffold.stageparser(sub.add_parser('assemble_scaffold'))

    # Reference-based
    align_reads.stageparser(sub.add_parser('align_reads'))
    call_variants.stageparser(sub.add_parser('call_variants'))
    vcf_to_consensus.stageparser(sub.add_parser('vcf_to_consensus'))
    refine_assembly.stageparser(sub.add_parser('refine_assembly'))
    finalize_assembly.stageparser(sub.add_parser('finalize_assembly'))

    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    main()
