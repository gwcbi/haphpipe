# -*- coding: utf-8 -*-

from __future__ import print_function
from builtins import str
import sys
import os
import argparse

from Bio import SeqIO

from ..utils.sysutils import check_dependency, existing_file, existing_dir, args_params
from ..utils.sysutils import create_tempdir, remove_tempdir
from ..utils.sequtils import wrap, extract_amplicons
from ..utils.alignutils import  assemble_to_ref
from ..utils.alignutils import TilingRow, NucmerAlignment, EmptyAlignment
from ..utils.alignutils import align_nucmer, show_aligns


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--assembly_fa', type=existing_file, required=True,
                        help='Assembly containing ambiguous regions.')
    parser.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Reference file to impute from.')
    parser.add_argument('--maxgap', type=int,
                        help='Maximum number of sequential ambiguous positions before splitting amplicon')
    parser.add_argument('--out_aln',
                        help='Output a text file with alignment in haphpipe-specific format')
    parser.add_argument('--keep_tmp', action='store_true',
                        help='Additional options')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), 
                        default=sys.stdout, help='Output file')                        
    parser.set_defaults(func=impute_ref)    


def new_impute_ref(assembly_fa=None, ref_fa=None, outfile=None,
               maxgap=30, out_aln=None,
               keep_tmp=False,
    ):
    # Check dependencies
    check_dependency('nucmer')
    check_dependency('delta-filter')
    check_dependency('show-tiling')

    # Temporary directory
    tempdir = create_tempdir('impute_ref')
    
    scaffolds = {}
    
    for asm_chrom in SeqIO.parse(assembly_fa, 'fasta'):
        amp_g = extract_amplicons(asm_chrom.id, str(asm_chrom.seq), maxgap=maxgap)
        with open(os.path.join(tempdir, 'qry.fa'), 'w') as outh:
            for i, amp in enumerate(amp_g):
                print('>%s.%d' % (amp[0].split()[0], i), file=outh)
                print('%s' % amp[1], file=outh)
        
        chrom_scaffolds = assemble_to_ref(ref_fa, os.path.join(tempdir, 'qry.fa'), tempdir)
        if len(chrom_scaffolds) == 1:
            scaffolds[asm_chrom.id] = list(chrom_scaffolds.values())[0]
        else:
            raise PipelineStepError("Assembled multiple scaffolds from one initial chromosome")
    
    chroms = sorted(scaffolds.keys())
    for chrom in chroms:
        print('>%s\n%s' % (chrom, wrap(scaffolds[chrom].get_imputed())), file=outfile)
    
    if out_aln is not None:
        with open(out_aln, 'w') as outh:
            for chrom in chroms:
                print('>%s' % chrom, file=outh)
                print('\n'.join('\t'.join(row) for row in scaffolds[chrom].as_matrix()), file=outh)

def old_impute_ref(assembly_fa=None, ref_fa=None, outfile=None,
               maxgap=30, out_aln=None,
               keep_tmp=False,
    ):
    # Check dependencies
    check_dependency('nucmer')
    check_dependency('delta-filter')
    check_dependency('show-tiling')

    # Temporary directory
    tempdir = create_tempdir('impute_ref')
    
    asm_dict = {s.id:s for s in SeqIO.parse(assembly_fa, 'fasta')}
    ref_dict = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}
    print("Max gap: %d" % maxgap, file=sys.stderr)
      
    scaffolds = {}
    for chrom in sorted(asm_dict.keys()):
        # Extract amplicons
        amps = extract_amplicons(chrom, str(asm_dict[chrom].seq), maxgap=maxgap)
        tmp_amplicons_fa = os.path.join(tempdir, '%s.amplicons.fasta' % chrom)
        with open(tmp_amplicons_fa, 'w') as outh:
            for i, amp in enumerate(amps):
                print('>%s.%d' % (amp[0].split()[0], i), file=outh)
                print('%s' % amp[1], file=outh)
        
        print('Found %d amplicons for %s' % (i+1, chrom), file=sys.stderr)
        
        # Align amplicons
        fil, til = align_nucmer(tmp_amplicons_fa, ref_fa, tempdir)
        
        # Parse tiling
        tr_list = [TilingRow(l) for l in open(til, 'rU')]
        
        # Get references
        refs = list(set([tr.ref for tr in tr_list]))
        if not len(refs) == 1:
            raise PipelineStepError("Amplicons align to multiple references")

        ref = refs[0]        
        ranked = sorted(tr_list, key=lambda x:x.qry_alen, reverse=True)

        # Extract nucmer alignments
        nuc_alns = []
        for tr in ranked:
            out1 = show_aligns(tr.ref, tr.qry, fil)
            nuc_alns.append(NucmerAlignment(out1.strip('\n').split('\n')))
        
        cur_scaf = EmptyAlignment(str(ref_dict[ref].seq).lower())
        for na in reversed(nuc_alns):
            cur_scaf = cur_scaf.combine_alignments(na)
        
        scaffolds[chrom] = cur_scaf
    
    for chrom in sorted(asm_dict.keys()):
        print('>%s\n%s' % (chrom, wrap(scaffolds[chrom].get_imputed())), file=outfile)
    
    if out_aln is not None:
        with open(out_aln, 'w') as outh:
            for chrom in sorted(asm_dict.keys()):
                print('>%s' % chrom, file=outh)
                print('\n'.join('\t'.join(row) for row in scaffolds[chrom].as_matrix()), file=outh)

def impute_ref(assembly_fa=None, ref_fa=None, outfile=None,
               maxgap=30, out_aln=None,
               keep_tmp=False,
    ):
    new_impute_ref(assembly_fa=assembly_fa, ref_fa=ref_fa, outfile=outfile, maxgap=maxgap, out_aln=out_aln, keep_tmp=keep_tmp)
    old_impute_ref(assembly_fa=assembly_fa, ref_fa=ref_fa, outfile=outfile, maxgap=maxgap, out_aln=out_aln, keep_tmp=keep_tmp)    


"""        
        for tr in tr_list:
            tr.ref_align(fil)
        
        _, imputed = overlapper(ref_dict, tr_list)
        if len(imputed) == 1:
            r, plist = imputed.popitem()
            all_imp.append((chrom, ''.join(plist)))
        else:
            print >>sys.stderr, 'WARNING: Assembly chromosome %s mapped to multiple references: %s' % (chrom, ', '.join(imputed.keys()))
            for i,t in enumerate(imputed.iteritems()):
                all_imp.append(('%s.%d' % (chrom,i), ''.join(t[1])))
"""    

'''
def impute_ref_parser(parser):
    parser.add_argument('--assembly_fa', type=existing_file, required=True,
                        help='File containing contigs, fasta format')
    parser.add_argument('--ref_fa', type=existing_file, required=True,
                        help='Prefix for reference fasta')
    parser.add_argument('--keep_tmp', action='store_true',
                        help='Additional options')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), 
                        default=sys.stdout, help='Output file')                        
    parser.set_defaults(func=impute_ref)                        
'''


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Impute ambiguous position from reference.')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
