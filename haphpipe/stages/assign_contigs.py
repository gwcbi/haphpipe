# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import str
import sys
import os
import argparse
from collections import defaultdict
from subprocess import Popen, PIPE
from Bio import SeqIO

from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.sequtils import wrap
from ..utils.helpers import merge_interval_list


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

'''
fmt6_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
             "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
'''

def blast_tophits(query, reference):
    p = Popen(['blastn',
               '-db', reference,
               '-outfmt', '6',
               '-num_alignments', '1'
               ], 
               stdin=open(query, 'rU'), stdout=PIPE, stderr=PIPE)
    o, e = p.communicate()
    if e.strip('\n'):
        print(e, file=sys.stderr)
    return [l.split('\t') for l in o.strip('\n').split('\n')]


def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--contigs_fa', type=existing_file, required=True,
                        help='File containing contigs, fasta format')
    group1.add_argument('--ref_db', required=True,
                        help='''Prefix for reference blast database''')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')                        
    parser.set_defaults(func=assign_contigs)


def assign_contigs(contigs_fa=None, ref_db=None, outdir='.'):
    """ Pipeline step to blast contigs against subtype database """
    # Check for executable
    check_dependency('blastn')    
    
    # Outputs
    out1 = os.path.join(outdir, 'blast.hits')
    out2 = []
    out3 = os.path.join(outdir, 'subtype.config')
    
    hits = blast_tophits(contigs_fa, ref_db)
    with open(out1, 'w') as outh:
        print('\n'.join('\t'.join(h) for h in hits), file=outh)
    
    contig_dict = {s.id:s for s in SeqIO.parse(contigs_fa, 'fasta')}
    by_subtype = defaultdict(list)
    for h in hits:
        by_subtype[h[1].split('.')[0]].append(contig_dict[h[0]])
    
    subtypes = sorted(by_subtype.keys())
    for st in subtypes:
        st_contigs = os.path.join(outdir, '%s.contigs.fa' % st)
        out2.append((st, st_contigs))
        print('%s%d' % (st.ljust(15), len(by_subtype[st])), file=sys.stderr)
        with open(st_contigs, 'w') as outh:
            for s in by_subtype[st]:
                print('>%s\n%s' % (s.id, wrap(str(s.seq))), file=outh)
    
    with open(out3, 'w') as outh:
        for t in out2:    
            print('%s\t%s' % t, file=outh)
    
    return out1, out2, out3


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assign contigs to subtypes')
    assign_contigs_parser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
