#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from subprocess import check_output
from collections import defaultdict

from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.gtfparse import gtf_parser, GTFRow

from ..utils.sysutils import PipelineStepError


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--consensus_fa', type=existing_file, required=True,
                        help='''Fasta file containing consensus sequence. This must be the
                                reference used to align reads. The sequence may be padded
                                with "N"s to match reference.''')
    group1.add_argument('--annotations_gtf', type=existing_file, required=True,
                        help='''GTF file with annotated regions''')
    group1.add_argument('--align_bam', type=existing_file, required=True,
                        help='''Aligned reads (BAM)''')                                
    group1.add_argument('--variants_vcf', type=existing_file, required=True,
                        help='''Variants (VCF)''')
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=post_assembly)


def reads_per_region(reg, alnbam):
    cmd = 'samtools view -f 66 -F 268 %s \'%s:%s-%s\' | wc -l' % (alnbam, reg[0], reg[1], reg[2])
    o = check_output(cmd, shell=True)
    return int(o.strip())

def snps_per_region(reg, vcfgz):
    cmd = [
        """bcftools view -i 'TYPE="snp" & DP>100' %s '%s:%s-%s'""" % (vcfgz, reg[0], reg[1], reg[2]),
        'bcftools stats',
        """grep '^SN.*number of SNPs'""",
        """perl -ne '$_=~/(\d+)$/; print "$1\n"'""",
    ]
    o = check_output(' | '.join(cmd), shell=True)
    return int(o.strip())    

def reads_per_amplicon(conbed, alnbam):
    ret = defaultdict(dict)
    bedlines = (l.strip('\n').split('\t') for l in open(conbed, 'rU'))
    for t in bedlines:
        ret[t[0]][t[3]] = reads_per_region((t[0], t[1], t[2]), alnbam)
    return ret

def post_assembly(consensus_fa=None, annotations_gtf=None, align_bam=None, variants_vcf=None,
        debug=False,
    ):
    check_dependency('bcftools')
    
    summary = [
        ['chrom', 'feature', 'start', 'end', 'length',
         'nfrags', 'nsnps',
         'call_reg', 'call_len',
         'match', 'mismatch', 'rgap', 'qgap', 'uncalled'],
    ]
    for gr in gtf_parser(annotations_gtf):
        gr.attrs['nfrags'] = reads_per_region((gr.chrom, gr.start, gr.end), align_bam)
        gr.attrs['nsnps'] = snps_per_region((gr.chrom, gr.start, gr.end), variants_vcf)        
        print gr
        summary.append([
            gr.chrom, gr.attrs['name'], gr.feature, gr.start, gr.end,
            gr.end - gr.start,
            gr.attrs['nfrags'], gr.attrs['nsnps'],      
            gr.attrs['call_reg'] if 'call_reg' in gr.attrs else '.',
            gr.attrs['call_len'] if 'call_len' in gr.attrs else '.',
            gr.attrs['match'] if 'match' in gr.attrs else '.',
            gr.attrs['mismatch'] if 'mismatch' in gr.attrs else '.',
            gr.attrs['rgap'] if 'rgap' in gr.attrs else '.',
            gr.attrs['qgap'] if 'qgap' in gr.attrs else '.',
            gr.attrs['uncalled'] if 'uncalled' in gr.attrs else '.',
            ])
    
    for s in summary:
        print '\t'.join(str(f) for f in s)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix consensus sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))

