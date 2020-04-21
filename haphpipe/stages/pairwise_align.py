#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import str
from builtins import range
import sys
import os
import argparse
import json
from collections import defaultdict

from Bio import SeqIO

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils
from haphpipe.utils import gtfparse
from haphpipe.utils.gtfparse import GTFRow
from haphpipe.utils import blastalign as baln


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


def matching_refseq(refseqs, k):
    if k in refseqs:
        return refseqs[k]
    else:
        if len(refseqs) > 1:
            msg = 'No match for "%s" in reference sequences (%s)'
            msg = msg % (k, ','.join(refseqs.keys()))
            raise sysutils.PipelineStepError(msg)
        return next(iter(refseqs.values()))

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--amplicons_fa',
                        type=sysutils.existing_file, required=True,
                        help='Assembled amplicons (fasta)')
    group1.add_argument('--ref_fa',
                        type=sysutils.existing_file, required=True,
                        help='Reference fasta file')
    group1.add_argument('--ref_gtf',
                        type=sysutils.existing_file, required=True,
                        help='''GTF format file containing amplicon regions.
                                Primary and alternate coding regions should be
                                provided in the attribute field (for amino
                                acid alignment).''')
    group1.add_argument('--outdir',
                        type=sysutils.existing_dir, default='.',
                        help='Output directory')
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=pairwise_align)

    
def pairwise_align(
        amplicons_fa=None, ref_fa=None, ref_gtf=None, outdir='.',
        keep_tmp=False, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to align amplicons to reference

    Args:
        amplicons_fa (str): Path to fasta file with amplicon sequences
        ref_fa (str): Path to reference fasta file
        ref_gtf (str): Path to reference GTF file with amplicons
        outdir (str): Path to output directory
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out_aln (str): Path to alignment in JSON format

    """
    # Check dependencies
    sysutils.check_dependency('blastx')

    # Outputs
    out_aln = os.path.join(outdir, 'alignments.json')

    # Temporary directory
    tempdir = sysutils.create_tempdir('pairwise_align', None, quiet, logfile)

    # Load reference sequence(s)
    refseqs = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    
    # Load amplicons from GTF file
    amps = [gl for gl in gtfparse.gtf_parser(ref_gtf) if
            gl.feature == 'amplicon']
    ampdict = {(gl.chrom, gl.attrs['name']):gl for gl in amps}
    
    out_json = {
        'aa_alignments': {},
        'nuc_alignments': {},
        'padded_alignments': {},
        'padded_gtf': [],
    }
    # {(sid, ref): [(reg, list(alignment)), ...], ...}
    ###--Uzma--is ^ above what defauldict is supposed to be formatted like?
    all_nuc_aln = defaultdict(list)

    for amprec in SeqIO.parse(amplicons_fa, 'fasta'):
        # Get amplicon reference and region from sequence ID
        aid = sequtils.parse_seq_id(amprec.id)
        # Find the GTF line used to orient this amplicon
        try:
            ###--Uzma--as in previous stage, what does aid[...] do?
            gl = ampdict[(aid['ref'], aid['reg'])]
        except KeyError:
            poss_gl = [t for t in ampdict.keys() if t[1] == aid['reg']]
            gl = ampdict[poss_gl[0]]

        # Start and stop for primary coding region
        pri_s = int(gl.attrs['primary_cds'].split('-')[0]) - 1
        pri_e = int(gl.attrs['primary_cds'].split('-')[1])  
        # Start and stop for additional coding regions
        altcds = []        
        if 'alt_cds' in gl.attrs:            
            for x in gl.attrs['alt_cds'].split(','):
                altcds.append(((int(x.split('-')[0]) - 1), int(x.split('-')[1])))
        
        # Align using amino acids
        refseq = matching_refseq(refseqs, aid['ref'])
        alnobj, nuc_aln = baln.alignAA(
            refseq,
            amprec,
            (pri_s, pri_e),
            altcds,
            tempdir,
            quiet
        )
        # prialn is a BlastxAlignment object with amplicon aligned to primary cds
        # merged is a nucleotide alignment over the full amplicon, with unaligned regions
        # aligned using alternate cds or nucleotide alignments

        all_nuc_aln[(aid['sid'], aid['ref'])].append((aid['reg'], nuc_aln))
        jid = 'sid|%s|ref|%s|reg|%s|' % (aid['sid'], aid['ref'], aid['reg'])
        out_json['aa_alignments'][jid] = alnobj.aa_align
        out_json['nuc_alignments'][jid] = nuc_aln
    
    # Full sequence with padding
    for sid, ref in list(all_nuc_aln.keys()):
        _refseq = matching_refseq(refseqs, ref)
        # New name and new alignment
        newname = 'sid|%s|ref|%s|' % (sid, _refseq.id)
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
                for p in range(rpos, seg[0][0]):
                    tmp.append((p, str(_refseq.seq[p]), '*', qpos))
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
            gr.attrs.update(baln.get_seg_stats(seg))
            # Include called regions
            gr.attrs['call_reg'] = '%d-%d' % (gr.start, gr.end)
            gr.attrs['call_len'] = (gr.end - gr.start + 1)
            # Append to json object
            out_json['padded_gtf'].append(str(gr))
            rpos = seg[-1][0] + 1
        
        # Add padding for end of sequence
        if rpos < len(_refseq.seq):
            for p in range(rpos, len(_refseq.seq)):
                tmp.append((p, str(_refseq.seq[p]), '*', qpos))
                qpos += 1
        
        # Validate the alignment
        vseq = ''.join(t[2] for t in tmp if t[3] != -1)
        if baln.validate_alignment(tmp, _refseq.seq, vseq):
            if not quiet:
                print('%s alignment validation passed' % newname,
                      file=sys.stderr)
            out_json['padded_alignments'][newname] = tmp
    
    for s in out_json['padded_gtf']:
        if not quiet:
            print(s, file=sys.stdout)
    
    with open(out_aln, 'w') as outh:
        print(json.dumps(out_json), file=outh)
    
    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'pairwise_align', quiet, logfile)

    return out_aln


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Align consensus to an annotated reference.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
