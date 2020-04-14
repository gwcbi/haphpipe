#! /usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import argparse
import sys
from Bio import SeqIO
from bisect import bisect
import random

# Copyright (C) 2019 Matthew L. Bendall
#
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


parser = argparse.ArgumentParser()
parser.add_argument("-input", help="Fasta file from PredictHaplo. This is after PHparser.py is ran.")
parser.add_argument("-out", help="Fasta output. Used to put the PredictHaplo sequences in FASTA format")
parser.add_argument("-num", type=int, help="Desired number of sequences to generate")
args = parser.parse_args()


seqs = []
freqs = []
fasta = []
cumfreq = []


equal_len = True

for s in SeqIO.parse(args.input, 'fasta'):
  seqs.append(s)
  
freqs = [float(s.description.split()[1].split('=')[1]) for s in seqs]
freqs = [old_div(f,sum(freqs)) for f in freqs] 
  
if len(seqs) != len(freqs):
    equal_len = False
    print("Length of sequences is not the same as length of frequencies. Something is wrong.")
else:
    pass

if len(seqs) == len(freqs):
  counter = 0
  for f in freqs:
    counter += f
    cumfreq.append(counter)

if len(freqs) == len(cumfreq):
  for n in range(args.num):
    s = random.random()
    index = bisect(cumfreq, s)
    fasta.append(seqs[index])

if len(fasta) == args.num:
  fafile = open(args.out, 'w')
  for i,seq in enumerate(fasta):
    print('>%s.%05d %s' % (seq.description.split()[0], i, seq.description.split()[1]), file=fafile)
    print('%s' % str(seq.seq), file=fafile)

  print("FASTA file comprised of %s sequences is completed for %s." %(args.num, args.input), file=sys.stderr)
