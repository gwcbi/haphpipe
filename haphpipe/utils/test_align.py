# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import re

from haphpipe.utils.alignutils import *

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


def issequential(l):
    tmp = sorted(l)
    for i,v in enumerate(tmp[:-1]):
        if not v+1 == tmp[i+1]:
            return False
    return True


print("Alignment")
print('%s%s' % ('full ref:'.ljust(40), 'acgtacgtacgtacgt'))
print()
print('%s%s' % ('aligned ref:'.ljust(40), 'ac-gtac'))
print('%s%s' % ('aligned qry:'.ljust(40), 'accg--c'))

ealn = EmptyReferenceAlignment('acgtacgtacgtacgt')
aln = ReferenceAlignment('ac-gtac', 'accg--c')
assert aln.rseq() == 'acgtac'
assert aln.qseq() == 'accgc'
m = ealn.merge_alignments(aln)
assert m.rseq() == 'acgtacgtacgtacgt'
assert m.qseq() == 'accgc??????????'
assert m.raln() == 'ac.gtacgtacgtacgt'
assert m.qaln() == 'accg..c??????????'

print()
print('%s%s' % ('ref:'.ljust(40), m.rseq()))
print('%s%s' % ('qry:'.ljust(40), m.qseq()))
print()
print('%s%s' % ('refaln:'.ljust(40), m.raln()))
print('%s%s' % ('qryaln:'.ljust(40), m.qaln()))

# Offset
# 'acgtacgtacgtacgt'
#    'tacgtacgt'
print("Offset")
print('%s%s' % ('full ref:'.ljust(40), 'acgtacgtacgtacgt'))
print()
print('%s%s' % ('aligned ref:'.ljust(40), '   ta--cgtacgt'))
print('%s%s' % ('aligned qry:'.ljust(40), '   tacccgt--gt'))

ealn = EmptyReferenceAlignment('acgtacgtacgtacgt')
aln = ReferenceAlignment('ta--cgtacgt', 'tacccgt--gt')
assert aln.rseq() == 'tacgtacgt'
assert aln.qseq() == 'tacccgtgt'
aln.adjust_ref_start(3)
m = ealn.merge_alignments(aln)
assert m.rseq() == 'acgtacgtacgtacgt', 'ERROR: ref seq should not change'
assert m.qseq() == '???tacccgtgt????'
assert m.raln() == 'acgta..cgtacgtacgt'
assert m.qaln() == '???tacccgt..gt????'

print()
print('%s%s' % ('ref:'.ljust(40), m.rseq()))
print('%s%s' % ('qry:'.ljust(40), m.qseq()))
print()
print('%s%s' % ('refaln:'.ljust(40), m.raln()))
print('%s%s' % ('qryaln:'.ljust(40), m.qaln()))



ref_fa_str = '''>HIV_B.K03455.HXB2_LAI_IIIB_BRU
TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG
GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC
TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG
AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGG
GGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTA
ACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTC
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG
GCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAG
AGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCA
'''

qry_fa_str = '''>HIV_B.SAMPLE
TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG
AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGG
GGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTA
ACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTC
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGG
'''
with open('ref.fa','w') as outh:
    print(ref_fa_str, file=outh)
with open('qry.fa','w') as outh:
    print(qry_fa_str, file=outh)

from haphpipe.utils.alignutils import align_nucmer, show_aligns

fil, til = align_nucmer('qry.fa', 'ref.fa',  '.')
out = show_aligns('HIV_B.SAMPLE', 'HIV_B.K03455.HXB2_LAI_IIIB_BRU', fil)

ealn = EmptyReferenceAlignment(''.join(ref_fa_str.split('\n')[1:-1]).lower())
if ealn.rseq() == ''.join(ref_fa_str.split('\n')[1:-1]).lower():
    print("Loaded ref successful")

aln = NucmerReferenceAlignment(out.split('\n'))
if aln.qseq() == ''.join(qry_fa_str.split('\n')[1:-1]).lower():
    print("Loaded alignment successful")

m = ealn.merge_alignments(aln)
teststr = '?' * 200
teststr += ''.join(qry_fa_str.split('\n')[1:-1]).lower()
teststr += '?' * 262
if m.qseq() == teststr:
    print('merging successful')

ref_fa_str = '''>HIV_B.K03455.HXB2_LAI_IIIB_BRU
TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG
GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC
TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG
AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGG
GGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTA
ACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTC
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG
GCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAG
AGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCA
'''

qry_fa_str = '''>HIV_B.SAMPLE
TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG
AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGG
GGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTA
ACTAGGGAACCCACTGCTTAAGCCTCAATAAAGGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTC
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGG
'''
with open('ref.fa','w') as outh:
    print(ref_fa_str, file=outh)
with open('qry.fa','w') as outh:
    print(qry_fa_str, file=outh)

from haphpipe.utils.alignutils import align_nucmer, show_aligns

fil, til = align_nucmer('qry.fa', 'ref.fa',  '.')
out = show_aligns( 'HIV_B.SAMPLE', 'HIV_B.K03455.HXB2_LAI_IIIB_BRU', fil)

ealn = EmptyReferenceAlignment(''.join(ref_fa_str.split('\n')[1:-1]).lower())
if ealn.rseq() == ''.join(ref_fa_str.split('\n')[1:-1]).lower():
    print("Loaded ref successful")

aln = NucmerReferenceAlignment(out.split('\n'))

if aln.qseq() == ''.join(qry_fa_str.split('\n')[1:-1]).lower():
    print("Loaded alignment successful")
else:
    print("Loaded alignment fail")

m = ealn.merge_alignments(aln)
teststr = '?' * 200
teststr += ''.join(qry_fa_str.split('\n')[1:-1]).lower()
teststr += '?' * 262
if m.qseq() == teststr:
    print('merging successful')
else:
    print('merging fail')
    
print(m.raln())
print(m.qaln())

#########
ref_fa_str = '''>HIV_B.K03455.HXB2_LAI_IIIB_BRU
TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGG
GCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGC
TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG
AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGG
GGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTA
ACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTC
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCG
GCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAG
AGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCA
'''

qry_fa_str = '''>HIV_B.SAMPLE1
TTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAG
AGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGG
GGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTA
ACTAGGGAACCCACTGCTTAAGCCTCAATAAAGGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTC
AGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGG
>HIV_B.SAMPLE2
GAAAATTTCTAGCAGTGGCGCCCGAACAGCGACCTGAAAGCGAAAGGGCAAACCAGAGGAGCTCTCTCGACGCAGGACTCG
GCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAG
AGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTA
'''

# GAAAATCTCTAGCAGTGG
# GAAAATTTCTAGCAGTGG


with open('ref.fa','w') as outh:
    print(ref_fa_str, file=outh)
with open('qry.fa','w') as outh:
    print(qry_fa_str, file=outh)

from haphpipe.utils.alignutils import align_nucmer, show_aligns

fil, til = align_nucmer('qry.fa', 'ref.fa',  '.')
print('\n'.join(l.strip() for l in open(til,'rU')))

ealn = EmptyReferenceAlignment(''.join(ref_fa_str.split('\n')[1:-1]).lower())
if ealn.rseq() == ''.join(ref_fa_str.split('\n')[1:-1]).lower():
    print("Loaded ref successful")

qry1_str = ''.join(qry_fa_str.split('>')[1].split('\n')[1:-1]).lower()
qry2_str = ''.join(qry_fa_str.split('>')[2].split('\n')[1:-1]).lower()

out = show_aligns('HIV_B.SAMPLE1', 'HIV_B.K03455.HXB2_LAI_IIIB_BRU', fil)
aln1 = NucmerReferenceAlignment(out.split('\n'))
if aln1.qseq() == qry1_str:
    print("Loaded alignment 1 successful")
else:
    print("Loaded alignment 1 fail")


out = show_aligns('HIV_B.SAMPLE2', 'HIV_B.K03455.HXB2_LAI_IIIB_BRU', fil)
aln2 = NucmerReferenceAlignment(out.split('\n'))
if aln2.qseq() == qry2_str:
    print("Loaded alignment 2 successful")
else:
    print("Loaded alignment 2 fail")

m = ealn.merge_alignments(aln1)
teststr = '?' * 200
teststr += qry1_str
teststr += '?' * 262
if m.qseq() == teststr:
    print('merging 1 successful')
else:
    print('merging 1 fail')

m = ealn.merge_alignments(aln2)
teststr = '?' * 620
teststr += qry2_str
teststr += '?' * 6
if m.qseq() == teststr:
    print('merging 2 successful')
else:
    print('merging 2 fail')


mboth = ealn.merge_alignments(aln2)
mboth = mboth.merge_alignments(aln1)

test2first = '?' * 200
test2first += qry1_str
test2first += qry2_str[18:]
test2first += '?' * 6

if mboth.qseq() == test2first:
    print('merging 2 first successful')
else:
    print('merging 2 first fail')

print('%s%s' % ('qry:'.ljust(40), mboth.qseq()))
print('%s%s' % ('expected:'.ljust(40), test2first))
print()
print('%s%s' % ('ref aln:'.ljust(40), mboth.raln()))
print('%s%s' % ('qry aln:'.ljust(40), mboth.qaln()))

mboth = ealn.merge_alignments(aln1)
mboth = mboth.merge_alignments(aln2)

test1first = '?' * 200
test1first += qry1_str[:-18]
test1first += qry2_str
test1first += '?' * 6

if mboth.qseq() == test1first:
    print('merging 1 first successful')
else:
    print('merging 1 first fail')

print('%s%s' % ('qry:'.ljust(40), mboth.qseq()))
print('%s%s' % ('expected:'.ljust(40), test1first))
print()
print('%s%s' % ('ref aln:'.ljust(40), mboth.raln()))
print('%s%s' % ('qry aln:'.ljust(40), mboth.qaln()))


#for p,c1,c2 in zip(range(len(test1first)), test1first,test2first):
#  if c1 != c2:
#    print '%d %s %s' % (p, c1, c2)
'''

if aln.qseq() == ''.join(qry_fa_str.split('\n')[1:-1]).lower():
    print "Loaded alignment successful"
else:
    print "Loaded alignment fail"

m = ealn.merge_alignments(aln)
teststr = '?' * 200
teststr += ''.join(qry_fa_str.split('\n')[1:-1]).lower()
teststr += '?' * 262
if m.qseq() == teststr:
    print 'merging successful'
else:
    print 'merging fail'
    
print m.raln()
print m.qaln()





'''