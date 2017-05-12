#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import json
import argparse
from collections import defaultdict, Counter
from itertools import chain

from Bio import SeqIO

from ..utils.helpers import overlaps
from ..utils.blastalign import *

from ..utils.sysutils import PipelineStepError
from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.sysutils import create_tempdir, remove_tempdir
from ..utils.sequtils import wrap, parse_seq_id
from ..utils.gtfparse import gtf_parser, GTFRow

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--amplicons_fa', type=existing_file, required=True,
                        help='Assembled amplicons (fasta)')
    group1.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Reference sequence (fasta)')
    group1.add_argument('--ref_gtf', type=existing_file, required=True,
                        help='''GTF format file containing amplicon regions. Primary and
                                alternate coding regions should be provided in the
                                attribute field (for amino acid alignment).''')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    # group2 = parser.add_argument_group('Pairwise alignment options')
    # group2.add_argument('--seqname',
    #                     help='Name to append to sequences.')
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=pairwise_align)

    
def pairwise_align(amplicons_fa=None, ref_fa=None, ref_gtf=None, outdir='.',
        keep_tmp=False, debug=False,
    ):
    # Check dependencies
    check_dependency('blastx')
    
    # Temporary directory
    tempdir = create_tempdir('pairwise_align')
    print >>sys.stderr, 'Temporary directory: %s' % tempdir
    
    refseqs = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    
    # Load amplicons from GTF file
    amps = [gtf_line for gtf_line in gtf_parser(ref_gtf) if gtf_line.feature == 'amplicon']
    ampdict = {(gl.chrom, gl.attrs['name']):gl for gl in amps}
    
    out_json = {
        'aa_alignments': {},
        'nuc_alignments': {},
        'padded_alignments': {},
        'padded_gtf': [],
    }
    all_nuc_aln = defaultdict(list) # {(sid, ref): [(reg, list(alignment)), ...], ...}
    
    for amprec in SeqIO.parse(amplicons_fa, 'fasta'):
        # Get amplicon reference and region from sequence ID
        aid = parse_seq_id(amprec.id)
        # Find the GTF line used to orient this amplicon
        gl = ampdict[(aid['ref'], aid['reg'])]
        # chrom, start, end, reg = bedline[0], int(bedline[1]), int(bedline[2]), bedline[3]
        # Start and stop for primary coding region
        pri_s = int(gl.attrs['primary_cds'].split('-')[0]) - 1
        pri_e = int(gl.attrs['primary_cds'].split('-')[1])  
        # Start and stop for additional coding regions
        altcds = []        
        if 'alt_cds' in gl.attrs:            
            for x in gl.attrs['alt_cds'].split(','):
                altcds.append(((int(x.split('-')[0]) - 1), int(x.split('-')[1])))
        
        # Align using amino acids
        alnobj, nuc_aln = alignAA(
            refseqs[aid['ref']],
            amprec,
            (pri_s, pri_e),
            altcds,
            tempdir
        )
        # prialn is a BlastxAlignment object with amplicon aligned to primary cds
        # merged is a nucleotide alignment over the full amplicon, with unaligned regions
        # aligned using alternate cds or nucleotide alignments

        all_nuc_aln[(aid['sid'], aid['ref'])].append((aid['reg'], nuc_aln))
        jid = 'sid|%s|ref|%s|reg|%s' % (aid['sid'], aid['ref'], aid['reg'])
        out_json['aa_alignments'][jid] = alnobj.aa_align
        out_json['nuc_alignments'][jid] = nuc_aln
    
    # Full sequence with padding
    for sid, ref in all_nuc_aln.keys():
        # New name and new alignment
        newname = 'sid|%s|ref|%s' % (sid, ref)
        tmp = []
        # Sort all segments by the start position
        segments = sorted(all_nuc_aln[(sid, ref)], key=lambda x:x[1][0][0])
        rpos = qpos = 0
        for sname, seg in segments:
            gr = GTFRow()
            gr.chrom, gr.source, gr.feature = (newname, 'haphpipe', 'amplicon')
            gr.score, gr.strand, gr.frame = ('.', '+', '.')
            gr.attrs['name'] = sname
                        
            # Pad up to first position of segment
            if rpos < seg[0][0]:
                for p in xrange(rpos, seg[0][0]):
                    tmp.append((p, str(refseqs[ref].seq[p]), '*', qpos))
                    qpos += 1
            gr.start = qpos + 1
            for t in seg:
                if t[3] == -1:
                    tmp.append(t)
                else:
                    tmp.append((t[0], t[1], t[2], qpos))
                    qpos += 1
            # Add annotation line
            gr.end = qpos
            # Include statistics in attributes
            gr.attrs.update(get_seg_stats(seg))
            # Include called regions
            gr.attrs['call_reg'] = '%d-%d' % (gr.start, gr.end)
            gr.attrs['call_len'] = (gr.end - gr.start + 1)
            # Append to json object
            out_json['padded_gtf'].append(str(gr))
            rpos = seg[-1][0] + 1
        
        # Add padding for end of sequence
        if rpos < len(refseqs[ref].seq):
            for p in xrange(rpos, len(refseqs[ref].seq)):
                tmp.append((p, str(refseqs[ref].seq[p]), '*', qpos))
                qpos += 1
        
        # Validate the alignment
        vseq = ''.join(t[2] for t in tmp if t[3] != -1)
        if validate_alignment(tmp, refseqs[ref].seq, vseq):
            print >>sys.stderr, '%s alignment validation passed' % newname
            out_json['padded_alignments'][newname] = tmp
    
    for s in out_json['padded_gtf']:
        print >>sys.stdout, s
    
    with open(os.path.join(outdir, 'alignments.json'), 'w') as outh:
        print >>outh, json.dumps(out_json)
    
    if not keep_tmp:
        remove_tempdir(tempdir, 'pairwise_align')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fix consensus sequence')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
