#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
from past.utils import old_div
import argparse
import sys
import os

from haphpipe.utils import sysutils

__author__ = 'Keylie M. Gibson'
__copyright__ = "Copyright (C) 2019 Keylie M. Gibson"


def stageparser(parser):
    """ Add stage-specific options to argparse parser

    Args:
        parser (argparse.ArgumentParser): ArgumentParser object

    Returns:
        None

    """
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--haplotypes_fa', type=sysutils.existing_file,
                        required=True,
                        help='Haplotype file created by PredictHaplo.')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory.')

    group2 = parser.add_argument_group('Parser options')
    group2.add_argument('--prefix',
                        help='Prefix to add to sequence names')
    group2.add_argument('--keep_gaps', action='store_true',
                        help='Do not remove gaps from alignment')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=ph_parser)


def ph_parser(
        haplotypes_fa=None, outdir='.',
        prefix=None, keep_gaps=False,
        quiet=False, logfile=None, debug=False,
    ):
    """

    Args:
        haplotypes_fa (str): Path to haplotype file from PredictHaplo (fasta-ish)
        outdir (str): Path to output directory
        prefix (str): Prefix to add to sequence names
        keep_gaps (bool): Do not remove gaps from alignment
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:

    """
    summary_txt = open(os.path.join(outdir, 'ph_summary.txt'), 'w')
    newseq_fa = open(os.path.join(outdir, 'ph_haplotypes.fna'), 'w')

    num_hap = 0
    freq = []
    fasta = []
    newseq = None
    for l in open(haplotypes_fa, 'r'):
        l = l.strip('\n')
        if l.startswith('>'):
            num_hap += 1
            if newseq is not None:
                fasta.append(newseq)
            newseq = [l.strip(">"), None, ""]
        elif l.startswith(';'):
            parts = l.strip(';').split(':')
            if parts[0] == 'Freq':
                freq.append(float(parts[1]))
                newseq[1] = float(parts[1])
            else:
                pass
        else:
            newseq[2] += l.strip('\n')

    fasta.append(newseq)

    if len(fasta) == num_hap:
        sysutils.log_message("Number of haplotypes is correct.\n", quiet,
                             logfile)
        freq_sqrd = [x ** 2 for x in freq]
        freq_sqrd_sum = sum(freq_sqrd)

        hap_div = ((old_div(7000, (7000 - 1))) * (1 - freq_sqrd_sum))

        print("PH_num_hap %s" % num_hap, file=summary_txt)
        print("PH_hap_diversity %s" % hap_div, file=summary_txt)

        seqlen = len(fasta[0][2])
        equal_len = True
        for seq in fasta:
            sl = len(seq[2])
            if sl != seqlen:
                sysutils.log_message(
                    "Sequence length is different for each haplotype.\n",
                    quiet, logfile
                )
                equal_len = False
            else:
                pass
        if equal_len == True:
            print("PH_seq_len %s" % seqlen, file=summary_txt)

        for sub_list in fasta:
            if prefix is None:
                print('>%s Freq=%s' % (sub_list[0], sub_list[1]),
                      file=newseq_fa)
            else:
                print('>%s_%s Freq=%s' % (prefix, sub_list[0], sub_list[1]),
                      file=newseq_fa)
            if keep_gaps:
                print("%s" % (sub_list[2]), file=newseq_fa)
            else:
                print("%s" % (sub_list[2].replace('-', "")), file=newseq_fa)

        sysutils.log_message(
            "Summary and FASTA file completed for %s.\n" % haplotypes_fa,
            quiet, logfile
        )

    summary_txt.close()
    newseq_fa.close()

def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Call variants.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()

