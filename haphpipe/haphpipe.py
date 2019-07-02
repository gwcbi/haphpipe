#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import
import sys
import argparse

from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import ArgumentDefaultsHelpFormatterSkipNone as HF

from haphpipe._version import VERSION

# Reads stages
from haphpipe.stages import sample_reads
from haphpipe.stages import trim_reads
from haphpipe.stages import join_reads
from haphpipe.stages import ec_reads
# Assemble stages
from haphpipe.stages import assemble_denovo
from haphpipe.stages import assemble_amplicons
from haphpipe.stages import assemble_scaffold
from haphpipe.stages import align_reads
from haphpipe.stages import call_variants
from haphpipe.stages import vcf_to_consensus
from haphpipe.stages import refine_assembly
from haphpipe.stages import finalize_assembly
# Haplotype stages
from haphpipe.stages import predict_haplo
from haphpipe.stages import ph_parser
# Annotate stages
from haphpipe.stages import pairwise_align
from haphpipe.stages import post_assembly
from haphpipe.stages import extract_pairwise
from haphpipe.stages import annotate_from_ref

BASE_USAGE = '''
Program: haphpipe (haplotype and phylodynamics pipeline)
Version: %(version)s

Commands:
 -- Reads
    sample_reads             subsample reads using seqtk
    trim_reads               trim reads using Trimmomatic
    join_reads               join reads using FLASh
    ec_reads                 error correct reads using SPAdes

 -- Assemble
    assemble_denovo          assemble reads denovo
    assemble_amplicons       assemble contigs to amplicon regions
    assemble_scaffold        assemble contigs to genome
    align_reads              align reads to reference
    call_variants            call variants
    vcf_to_consensus         create consensus sequence from VCF
    refine_assembly          iterative refinement: align - variants - consensus
    finalize_assembly        finalize consensus sequence

 -- Haplotype
    predict_haplo            assemble haplotypes with PredictHaplo
    ph_parser                parse output from PredictHaplo.

 -- Annotate
    pairwise_align           align consensus to an annotated reference
    extract_pairwise         extract sequence regions from pairwise alignment
    annotate_from_ref        annotate consensus from reference annotation
'''


class _BaseHelpAction(argparse.Action):
    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help=None
                 ):
        super(_BaseHelpAction, self).__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        print(BASE_USAGE % {'version': VERSION}, file = sys.stdout)
        parser.exit()


def console():
    parser = argparse.ArgumentParser(formatter_class=HF, add_help=False)
    parser.add_argument('-h', '--help', action=_BaseHelpAction)
    sub = parser.add_subparsers()

    # Reads stages
    sample_reads.stageparser(
        sub.add_parser('sample_reads', formatter_class=HF)
    )
    trim_reads.stageparser(
        sub.add_parser('trim_reads', formatter_class=HF)
    )
    join_reads.stageparser(
        sub.add_parser('join_reads', formatter_class=HF)
    )
    ec_reads.stageparser(
        sub.add_parser('ec_reads', formatter_class=HF)
    )

    # Assemble stages
    assemble_denovo.stageparser(
        sub.add_parser('assemble_denovo', formatter_class=HF)
    )
    assemble_amplicons.stageparser(
        sub.add_parser('assemble_amplicons', formatter_class=HF)
    )
    assemble_scaffold.stageparser(
        sub.add_parser('assemble_scaffold', formatter_class=HF)
    )
    align_reads.stageparser(
        sub.add_parser('align_reads', formatter_class=HF)
    )
    call_variants.stageparser(
        sub.add_parser('call_variants', formatter_class=HF)
    )
    vcf_to_consensus.stageparser(
        sub.add_parser('vcf_to_consensus', formatter_class=HF)
    )
    refine_assembly.stageparser(
        sub.add_parser('refine_assembly', formatter_class=HF)
    )
    finalize_assembly.stageparser(
        sub.add_parser('finalize_assembly', formatter_class=HF)
    )

    # Haplotype
    predict_haplo.stageparser(
        sub.add_parser('predict_haplo', formatter_class=HF)
    )
    ph_parser.stageparser(
        sub.add_parser('ph_parser', formatter_class=HF)
    )

    # Annotate stages
    pairwise_align.stageparser(
        sub.add_parser('pairwise_align', formatter_class=HF)
    )
    extract_pairwise.stageparser(
        sub.add_parser('extract_pairwise', formatter_class=HF)
    )
    annotate_from_ref.stageparser(
        sub.add_parser('annotate_from_ref', formatter_class=HF)
    )
    # post_assembly.stageparser(
    #     sub.add_parser('post_assembly', formatter_class=HF)
    # )

    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
