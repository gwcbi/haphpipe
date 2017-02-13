#! /usr/bin/env python

import argparse
import sys
import re

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
  print "Number of haplotypes is correct."
  
  freq_sqrd = [x**2 for x in freq]
  freq_sqrd_sum = sum(freq_sqrd)

  hap_div = ((7000 / (7000 - 1)) * (1 - freq_sqrd_sum))


  sumfile = open(args.sum, 'w')
  print >>sumfile, "QR_num_hap %s" %num_hap
  print >>sumfile, "QR_hap_diversity %s" %hap_div
  
  seqlen = len(fasta[0][2])
  equal_len = True
  for seq in fasta:
    sl = len(seq[2])
    if sl != seqlen:
      print "Sequence length is different for each haplotype."
      equal_len = False
    else:
      pass
  if equal_len == True:
    print >>sumfile, "QR_seq_len %s" %seqlen

  fafile = open(args.fa, 'w')
  for sub_list in fasta:
    print >>fafile, '>%s_%s Freq=%s' %(args.pre, sub_list[0], sub_list[1])
    print >>fafile, "%s" %(sub_list[2].replace('-', ""))


  print >>sys.stderr, "Summary and FASTA file completed."

elif len(fasta) != num_hap:
  print "Number of haplotypes is off. This needs checking before computing haplotype diversity statistic."
  sumfile = open(args.sum, 'w')
  
  seqlen = len(fasta[0][2])
  equal_len = True
  for seq in fasta:
    sl = len(seq[2])
    if sl != seqlen:
      print "Sequence length is different for each haplotype."
      equal_len = False
    else:
      pass
  if equal_len == True:
    print >>sumfile, "PH_seq_len %s" %seqlen
  fafile = open(args.fa, 'w')
  for sub_list in fasta:
    print >>fafile, ">%s Freq=%s" %(sub_list[0], sub_list[1])
    print >>fafile, "%s" %(sub_list[2])
    
  print >>sys.stderr, "FASTA file completed. Summary file had trouble with number of haplotypes. Stats could not be calculated."
