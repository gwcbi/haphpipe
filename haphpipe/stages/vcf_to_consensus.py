#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

from builtins import map
from builtins import str
from builtins import zip
from past.builtins import basestring

import os
import argparse
import gzip
import re

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils
from haphpipe.utils.helpers import cast_str


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


VCFCOLS = ['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
idx = {v:i for i,v in enumerate(VCFCOLS)}

def parse_vcf_sample(row, sampidx=0):
    """ Parse row in VCF file
    """
    if isinstance(row, basestring):
        row = row.strip('\n').split('\t')
    chrom = row[idx['CHROM']]
    start = int(row[idx['POS']])
    stop = start + len(row[idx['REF']]) - 1
    RA = [row[idx['REF']], ]
    AA = [alt for alt in row[idx['ALT']].split(',') if alt != '.']
    info = dict(kv.split('=') if '=' in kv else (kv, True) 
                    for kv in row[idx['INFO']].split(';'))
    fmt = row[idx['FORMAT']].split(':')
    svals = dict(list(zip(fmt, row[sampidx + 9].split(':'))))
    if 'AD' in svals: svals['AD'] = cast_str(svals['AD'])
    if 'DP' in svals: svals['DP'] = cast_str(svals['DP'])
    return chrom, start, stop, RA, AA, info, svals

def call_gt(RA, AA, svals, min_dp=1, major=0.5, minor=0.2):
    alleles = RA + AA
    if len(alleles) == 1:
        samp_dp = svals['DP'] if 'DP' in svals else 0
        return None if samp_dp < min_dp else alleles
    else:
        samp_ad = list(map(int, svals['AD'].split(',')))
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
    """ Add stage-specific options to argparse parser

    Args:
        parser (argparse.ArgumentParser): ArgumentParser object

    Returns:
        None

    """
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--vcf', type=sysutils.existing_file, required=True,
                        help='VCF file (created with all sites).')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    group1.add_argument('--sampidx', type=int, default=0,
                        help='Index for sample if multi-sample VCF')

    group2 = parser.add_argument_group('Variant options')
    group2.add_argument('--min_dp', type=int, default=1,
                        help='Minimum depth to call site')
    group2.add_argument('--major', type=float, default=0.5,
                        help='Allele fraction to make unambiguous call')
    group2.add_argument('--minor', type=float, default=0.2,
                        help='Allele fraction to make ambiguous call')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    # group3.add_argument('--debug', action='store_true',
    #                    help='Print commands but do not run')

    parser.set_defaults(func=vcf_to_consensus)


def vcf_to_consensus(
        vcf=None, outdir='.',
        sampidx=0, min_dp=1, major=0.5, minor=0.2,
        keep_tmp=False, quiet=False, logfile=None,
    ):
    """ Pipeline step to create consensus sequence from VCF

    Args:
        vcf (str): Path to variant calls (VCF)
        outdir (str): Path to output directory
        sampidx (int): Index for sample if multi-sample VCF
        min_dp (int): Minimum depth to call site
        major (float): Allele fraction to make unambiguous call
        minor (float): Allele fraction to make ambiguous call
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:

    """
    # Check inputs
    if vcf is None:
        raise sysutils.PipelineStepError('VCF file is required')

    # Outputs
    out_fasta = os.path.join(outdir, 'consensus.fna')

    sysutils.log_message('[--- vcf_to_consensus ---]\n', quiet, logfile)
    sysutils.log_message('VCF:          %s\n' % vcf, quiet, logfile)

    # Parse VCF
    chroms = []
    samples = []

    if os.path.splitext(vcf)[1] == '.gz':
        lines = (l.decode('utf-8').strip('\n') for l in gzip.open(vcf, 'rb'))
    else:
        lines = (l.strip('\n') for l in open(vcf, 'r'))

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
        raise sysutils.PipelineStepError(msg)

    chrom_ordered = [_[0] for _ in chroms]
    chroms = dict(chroms)
    newseqs = dict((c, ['.'] * chroms[c]) for c in list(chroms.keys()))
    imputed = dict((c, [''] * chroms[c]) for c in list(chroms.keys()))
    for l in lines:
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
                    newseqs[chrom][start-1] = sequtils.get_ambig(gt)
                    imputed[chrom][start-1] = sequtils.get_ambig(gt)
                else:
                    newseqs[chrom][start-1] = ''.join(gt[0])
                    imputed[chrom][start-1] = ''.join(gt[0])  
    # newseqs = imputed
    sysutils.log_message('Output FASTA: %s\n' % out_fasta, quiet, logfile)
    with open(out_fasta, 'w') as outh:
        for chrom in chrom_ordered:
            new_seqid = sequtils.update_seq_id(chrom, samples[sampidx])
            print('>%s SM:%s' % (new_seqid, samples[sampidx]), file=outh)
            ns = ''.join(newseqs[chrom])
            print(sequtils.wrap(ns.replace('.', 'n')), file=outh)

    return out_fasta


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Create consensus sequence from VCF',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()

