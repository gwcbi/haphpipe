#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import gzip
import re

from ..utils.sysutils import PipelineStepError
from ..utils.sysutils import existing_file, args_params
from ..utils.sequtils import wrap, get_ambig

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

VCFCOLS = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
idx = {v:i for i,v in enumerate(VCFCOLS)}

def parse_vcf_sample(row, sampidx=0):
    """ Parse row in VCF file
    """
    if type(row) is str:
        row = row.strip('\n').split('\t')
    chrom = row[idx['CHROM']]
    start = int(row[idx['POS']])
    stop = start + len(row[idx['REF']]) - 1
    RA = [row[idx['REF']], ]
    AA = [alt for alt in row[idx['ALT']].split(',') if alt != '.']
    # alleles = RA + AA # [row[idx['REF']]] + [alt for alt in row[idx['ALT']].split(',') if alt != '.']
    info = dict(kv.split('=') if '=' in kv else (kv, True) 
                    for kv in row[idx['INFO']].split(';'))
    # dp =  int(info['DP']) if 'DP' in info else 0
    fmt = row[idx['FORMAT']].split(':')
    svals = dict(zip(fmt, row[sampidx + 9].split(':')))
    return chrom, start, stop, RA, AA, info, svals

def call_gt(RA, AA, svals, min_dp=1, major=0.5, minor=0.2):
    alleles = RA + AA
    if len(alleles) == 1:
        samp_dp = svals['DP'] if 'DP' in svals else 0
        return None if samp_dp < min_dp else alleles
    else:
        samp_ad = map(int, svals['AD'].split(','))
        samp_dp = sum(samp_ad)
        if samp_dp < min_dp:
            return None
        assert len(samp_ad) == len(alleles), alleles
        samp_ad = sorted(zip(alleles, samp_ad), key=lambda x:x[1], reverse=True)
        if samp_ad[0][1] > samp_dp * major:
            return [samp_ad[0][0]]
        else:
            return [b for b,c in samp_ad if c >= samp_dp * minor]

def stageparser(parser):
    parser.add_argument('--vcf', type=existing_file, required=True,
                        help='VCF file (created with all sites).')
    parser.add_argument('--sampidx', type=int,
                        help='Index for sample if multi-sample VCF')
    parser.add_argument('--min_dp', type=int,
                        help='Minimum depth to call site')
    parser.add_argument('--major', type=float,
                        help='Allele fraction to make unambiguous call')
    parser.add_argument('--minor', type=float,
                        help='Allele fraction to make ambiguous call')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), 
                        default=sys.stdout, help='Output file')
    parser.set_defaults(func=vcf_to_fasta)


def vcf_to_fasta(vcf=None, sampidx=0, min_dp=1, major=0.5, minor=0.2, outfile=None):
    """ Call consensus sequence from VCF file
    """
    if vcf is None:
        raise PipelineStepError('VCF file is required')

    # Open filehandle for VCF
    if os.path.splitext(vcf)[1] == '.gz':
        fh = gzip.open(vcf, 'rb')
    else:
        fh = open(vcf, 'rU')
    
    chroms = []
    samples = []
    lines = (l.strip('\n') for l in fh)
    
    # Parse headers
    for l in lines:
        if l.startswith('##'):
            m = re.match('##contig=<ID=(\S+),length=(\d+)>', l)
            if m:
                chroms.append((m.group(1), int(m.group(2))))
        else:
            assert l.startswith('#')
            cols = l.strip('#').split('\t')
            samples = cols[9:]
            break
    
    if len(samples) <= sampidx:
        msg = 'Sample index %d does not exist. Samples: %s' % (sampidx, str(samples))
        raise PipelineStepError(msg)
    
    chroms = dict(chroms)
    newseqs = dict((c, ['.'] * chroms[c]) for c in chroms.keys())
    imputed = dict((c, [''] * chroms[c]) for c in chroms.keys())
    for l in lines:
        # chrom, start, stop, gt = parse_vcf_row(l, sampidx, min_dp, major, minor)    
        chrom, start, stop, RA, AA, info, svals = parse_vcf_sample(l, sampidx)
        gt = call_gt(RA, AA, svals, min_dp, major, minor)
        
        if gt is None:
            imputed[chrom][start-1] = RA[0].lower()
        else:
            if len(gt) == 1:
                newseqs[chrom][start-1] = gt[0]
                imputed[chrom][start-1] = gt[0]
            else:
                if all(len(_) == 1 for _ in gt):
                    newseqs[chrom][start-1] = get_ambig(gt)
                    imputed[chrom][start-1] = get_ambig(gt)
                else:
                    newseqs[chrom][start-1] = ''.join(gt[0])
                    imputed[chrom][start-1] = ''.join(gt[0])  
    # newseqs = imputed
    for chrom in sorted(newseqs.keys()):
        print >>outfile, '>%s SM:%s' % (chrom, samples[sampidx])
        ns = ''.join(newseqs[chrom])
        # print >>outfile, ns
        print >>outfile, wrap(ns.replace('.', 'n'))
    return outfile.name

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create consensus sequence from VCF')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
