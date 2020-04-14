#! /usr/bin/env python

from __future__ import print_function
from __future__ import division
from past.utils import old_div
import argparse
import sys
import re

# Copyright (C) 2019 Matthew L. Bendall

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
parser.add_argument("-input", help="input .fa file from QuasiRecomb")
parser.add_argument("-sum", help="first output. Used to put the QuasiRecomb Summary")
parser.add_argument("-fa", help="second output. Used to put the QuasiRecomb sequences in FASTA format without gaps")
parser.add_argument("-pre", help="prefix to add to fasta names")
args = parser.parse_args()


num_hap = 0
freq = []
fasta = []
newseq = None
for l in open(args.input, 'rU'):
  l = l.strip('\n')
  if l.startswith('>'):
    num_hap += 1
    parts = l.split('_')
    if newseq is not None:
        fasta.append(newseq)
    newseq = [parts[0].strip('>'), None, ""]
    
    freq.append(float(parts[1]))
    newseq[1] = float(parts[1])
  else:
    newseq[2] += l.strip('\n')



fasta.append(newseq)


if len(fasta) == num_hap:
  print("Number of haplotypes is correct.")
  
  freq_sqrd = [x**2 for x in freq]
  freq_sqrd_sum = sum(freq_sqrd)

  hap_div = ((old_div(7000, (7000 - 1))) * (1 - freq_sqrd_sum))


  sumfile = open(args.sum, 'w')
  print("QR_num_hap %s" %num_hap, file=sumfile)
  print("QR_hap_diversity %s" %hap_div, file=sumfile)
  
  seqlen = len(fasta[0][2])
  equal_len = True
  for seq in fasta:
    sl = len(seq[2])
    if sl != seqlen:
      print("Sequence length is different for each haplotype.")
      equal_len = False
    else:
      pass
  if equal_len == True:
    print("QR_seq_len %s" %seqlen, file=sumfile)

  fafile = open(args.fa, 'w')
  for sub_list in fasta:
    print('>%s_%s Freq=%s' %(args.pre, sub_list[0], sub_list[1]), file=fafile)
    print("%s" %(sub_list[2].replace('-', "")), file=fafile)


  print("Summary and FASTA file completed.", file=sys.stderr)

elif len(fasta) != num_hap:
  print("Number of haplotypes is off. This needs checking before computing haplotype diversity statistic.")
  sumfile = open(args.sum, 'w')
  
  seqlen = len(fasta[0][2])
  equal_len = True
  for seq in fasta:
    sl = len(seq[2])
    if sl != seqlen:
      print("Sequence length is different for each haplotype.")
      equal_len = False
    else:
      pass
  if equal_len == True:
    print("PH_seq_len %s" %seqlen, file=sumfile)
  fafile = open(args.fa, 'w')
  for sub_list in fasta:
    print(">%s Freq=%s" %(sub_list[0], sub_list[1]), file=fafile)
    print("%s" %(sub_list[2]), file=fafile)
    
  print("FASTA file completed. Summary file had trouble with number of haplotypes. Stats could not be calculated.", file=sys.stderr)
