#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import sys
# import os
# import json
import argparse
# from subprocess import check_output
# from collections import defaultdict

# from Bio import SeqIO
from Bio.Seq import Seq

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils
# from haphpipe.utils.gtfparse import gtf_parser, GTFRow
from haphpipe.utils.blastalign import load_slot_json

from ..utils.sysutils import PipelineStepError

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

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--align_json', type=sysutils.existing_file,
                        required=True,
                        help='''JSON file describing alignment (output of
                                pairwise_align stage)''')
    group1.add_argument('--outfile',
                        help='''Output file. Default is stdout''')
    group2 = parser.add_argument_group('Extract pairwise options')                                
    group2.add_argument('--outfmt', type=str, default='nuc_fa',
                        choices=['nuc_fa', 'aln_fa', 'amp_gtf', 'tsv',
                                  'prot_fa', ],
                        help='Format for output')
    group2.add_argument('--refreg',
                        help='''Reference region. String format is
                                ref:start-stop. For example, the region string
                                to extract pol when aligned to HXB2 is
                                HIV_B.K03455.HXB2:2085-5096''')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=extract_pairwise)


def extract_pairwise(align_json=None, outfile=None,
        outfmt=None, refreg=None,
        debug=False,
    ):
    outh = sys.stdout if outfile is None else open(outfile, 'w')
    
    if outfmt == 'nuc_fa' or outfmt == 'prot_fa':
        jaln = load_slot_json(align_json, 'padded_alignments')
        if refreg is None:
            for newname, alignment in list(jaln.items()):
                nucstr = ''.join(t[2] for t in alignment if t[3] != -1)
                nucstr = nucstr.replace('*', 'N')
                print('>%s' % newname, file=outh)                
                if outfmt == 'nuc_fa':
                    print(sequtils.wrap(nucstr), file=outh)
                else:
                    s = Seq(nucstr[:(old_div(len(nucstr),3))*3])
                    print(sequtils.wrap(str(s.translate())), file=outh)
        else:
            refmap = {sequtils.parse_seq_id(k)['ref']:k for k in list(jaln.keys())}
            chrom, ref_s, ref_e = sequtils.region_to_tuple(refreg)
            ref_s = ref_s - 1
            alignment = jaln[refmap[chrom]]
            
            # Get alignment start
            for aln_s in range(len(alignment)):
                if alignment[aln_s][0] == ref_s:
                    break
                while alignment[aln_s][3] == -1:
                    aln_s += 1
            
            # Get alignment end
            for aln_e in range(len(alignment)-1, -1, -1):
                if alignment[aln_e][0] == ref_e:
                    break
            while alignment[aln_e][3] == -1:
                aln_e += -1

            nucstr = ''.join(t[2] for t in alignment[aln_s:aln_e] if t[3] != -1)
            nucstr = nucstr.replace('*', 'N')
            print('>%s (%s)' % (refmap[chrom], refreg), file=outh)
            if outfmt == 'nuc_fa':            
                print(sequtils.wrap(nucstr), file=outh)
            else:
                s = Seq(nucstr[:(old_div(len(nucstr),3))*3])
                print(sequtils.wrap(str(s.translate())), file=outh)

    elif outfmt == 'aln_fa':
        jaln = load_slot_json(align_json, 'padded_alignments')
        for newname, alignment in list(jaln.items()):
            aid = sequtils.parse_seq_id(newname)
            rstr = ''.join(t[1] for t in alignment).replace('*', 'N')
            qstr = ''.join(t[2] for t in alignment).replace('*', 'N')
            print('>ref|%s|' % aid['ref'], file=outh)
            print(sequtils.wrap(rstr), file=outh)
            print('>sid|%s|' % aid['sid'], file=outh)
            print(sequtils.wrap(qstr), file=outh)

    elif outfmt == 'amp_gtf':
        jgtf = load_slot_json(align_json, 'padded_gtf')
        print('\n'.join(_ for _ in jgtf), file=outh)

    elif outfmt == 'tsv':
        jaln = load_slot_json(align_json, 'padded_alignments')
        for newname, alignment in list(jaln.items()):
            print('# %s' % newname, file=outh)
            for l in alignment:
                print('\t'.join(str(_) for _ in l), file=outh)


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Extract sequence regions from pairwise alignment.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
