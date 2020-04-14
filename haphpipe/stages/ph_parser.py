#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
from past.utils import old_div
import argparse
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
    parser.set_defaults(func=ph_parser)


def ph_parser(
        haplotypes_fa=None, outdir='.',
        prefix=None, keep_gaps=False,
        quiet=False, logfile=None,
    ):
    """

    Args:
        haplotypes_fa (str): Path to haplotype file from PredictHaplo (fasta-ish)
        outdir (str): Path to output directory
        prefix (str): Prefix to add to sequence names
        keep_gaps (bool): Do not remove gaps from alignment
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file

    Returns:

    """
    summary_txt = open(os.path.join(outdir, 'ph_summary.txt'), 'w')
    newseq_fa = open(os.path.join(outdir, 'ph_haplotypes.fna'), 'w')

    num_hap = 0
    freq = []
    fasta = []
    newseq = None
    ###--Uzma--I understand what each part is doing, but it's difficult to keep track of why because I'm not sure how the haplotype file 
    ###--is suposed to be formatted
    ph = os.path.basename(haplotypes_fa).split(".")[0]
    for l in open(haplotypes_fa, 'r'):
        l = l.strip('\n')
        if l.startswith('>'):
            num_hap += 1
            if newseq is not None:
                fasta.append(newseq)
            newseq = [ph, l.strip(">"), None, ""] ###--Uzma--What does  '...None, ""]' do?
        elif l.startswith(';'):
            parts = l.strip(';').split(':')
            if parts[0] == 'Freq':
                freq.append(float(parts[1]))
                newseq[2] = float(parts[1])
            else:
                pass
        else:
            newseq[3] += l.strip('\n')

    fasta.append(newseq)

    if len(fasta) == num_hap:
        sysutils.log_message("Number of haplotypes is correct.\n", quiet,
                             logfile)
        freq_sqrd = [x ** 2 for x in freq]
        freq_sqrd_sum = sum(freq_sqrd)

        ###--Uzma--Why 7000?
        hap_div = ((old_div(7000, (7000 - 1))) * (1 - freq_sqrd_sum))

        print("PH_num_hap %s" % num_hap, file=summary_txt)
        print("PH_hap_diversity %s" % hap_div, file=summary_txt)

        seqlen = len(fasta[0][-1])
        equal_len = True
        for seq in fasta:
            sl = len(seq[-1])
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
                print('>sid|%s_%s|reg|%s| Freq=%s' % (sub_list[0], sub_list[1], sub_list[0].split("_")[-1], sub_list[2]),
                      file=newseq_fa)
            else:
                print('>sid|%s_%s_%s|reg|%s| Freq=%s' % (prefix, sub_list[0], sub_list[1], sub_list[0].split("_")[-1], sub_list[2]),
                      file=newseq_fa)
            if keep_gaps:
                print("%s" % (sub_list[-1]), file=newseq_fa)
            else:
                print("%s" % (sub_list[-1].replace('-', "")), file=newseq_fa)

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
        description='Parse output from PredictHaplo.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()

