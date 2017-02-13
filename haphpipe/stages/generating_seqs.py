#! /usr/bin/env python

import argparse
import sys
from Bio import SeqIO
from bisect import bisect
import random

parser = argparse.ArgumentParser()
parser.add_argument("-input", help="Fasta file from PredictHaplo. This is after PHparser.py is ran.")
parser.add_argument("-out", help="Fasta output. Used to put the PredictHaplo sequences in FASTA format")
parser.add_argument("-num", help="Desired number of sequences to generate")
args = parser.parse_args()


seqs = []
freqs = []
fasta = []
cumfreq = []


equal_len = True

for s in SeqIO.parse(args.input, 'fasta'):
  seqs.append(s)
  
freqs = [float(s.description.split()[1].split('=')[1]) for s in seqs]
freqs = [f/sum(freqs) for f in freqs] 
  
if len(seqs) != len(freqs):
    equal_len = False
    print "Length of sequences is not the same as length of frequencies. Something is wrong."
else:
    pass

if len(seqs) == len(freqs):
  counter = 0
  for f in freqs:
    counter += f
    cumfreq.append(counter)

if len(freqs) == len(cumfreq):
  for n in range(int(args.num)):
    s = random.random()
    index = bisect(cumfreq, s)
    fasta.append(seqs[index])

SeqIO.write(fasta, args.out, 'fasta')
print >>sys.stderr, "FASTA file comprised of %s sequences is completed for %s." %(args.num, args.input)
