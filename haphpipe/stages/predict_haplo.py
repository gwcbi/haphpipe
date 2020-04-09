#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import absolute_import
from builtins import map
import sys
import os
import re
import argparse
from glob import glob
import shutil

from Bio import SeqIO

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils
from haphpipe.utils.sysutils import PipelineStepError


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


config_template = '''### configuration file for the HIVhaplotyper
### prefix
%(prefix)s
### filename of reference sequence (FASTA)
%(ref_fasta)s
### do_visualize (1 = true, 0 = false)
%(do_visualize)s
### filname of the aligned reads (sam format)
%(alignment)s
### have_true_haplotypes  (1 = true, 0 = false)
%(have_true_haplotypes)s
### filname of the true haplotypes (MSA in FASTA format) (fill in any dummy filename if there is no "true" haplotypes)
%(file_true_haplotypes)s
### do_local_analysis  (1 = true, 0 = false) (must be 1 in the first run)
%(do_local_analysis)s
### max_reads_in_window;
%(max_reads)s
### entropy_threshold
%(entropy_threshold)s
###reconstruction_start
%(reconstruction_start)s
###reconstruction_stop
%(reconstruction_stop)s
###min_mapping_qual
%(min_mapping_qual)s
###min_readlength
%(min_readlength)s
###max_gap_fraction (relative to alignment length)
%(max_gap_fraction)s
###min_align_score_fraction (relative to read length)
%(min_align_score_fraction)s
###alpha_MN_local (prior parameter for multinomial tables over the nucleotides)
%(alpha_MN_local)s
###min_overlap_factor (reads must have an overlap with the local reconstruction window of at least this factor times the window size)
%(min_overlap_factor)s
###local_window_size_factor (size of  local reconstruction window relative to the median of the read lengths)
%(local_window_size_factor)s
### max number of clusters (in the truncated Dirichlet process)
%(max_number_of_clusters)s
### MCMC iterations
%(mcmc_iterations)s
### include deletions (0 = no, 1 = yes)
%(include_deletions)s
'''


DEFAULTS = {
             'do_visualize': 1,
             'have_true_haplotypes': 0,
             'file_true_haplotypes': 'dummy.fasta',
             'do_local_analysis': 1,
             'max_reads': 10000,
             'entropy_threshold': '4e-2',
             'min_mapping_qual': 25,
             'max_gap_fraction': 0.05,
             'min_align_score_fraction': 0.35,
             'alpha_MN_local': 25,
             'min_overlap_factor': 0.85,
             'local_window_size_factor': 0.7,
             'max_number_of_clusters': 25,
             'mcmc_iterations': 501,
             'include_deletions': 1,
}

def rename_best(d, rn):
    fasta = glob(os.path.join(d, '%s*global*.fas' % rn))
    coords = [re.match('\S+global_(\d+)_(\d+).fas', f) for f in fasta]
    coords = [list(map(int, m.groups())) if m else [0,0] for m in coords]
    sizes = [c[1]-c[0] for c in coords]
    if max(sizes) == 0: print("WARNING: Max sizes is 0.", file=sys.stderr)
    bestidx = sizes.index(max(sizes))
    bestfile = '%s.best_%d_%d.fas' % (rn, coords[bestidx][0], coords[bestidx][1])
    if os.path.exists(os.path.join(d, bestfile)):
        print("WARNING: File %s exists." % os.path.join(d, bestfile), file=sys.stderr)
        os.unlink(os.path.join(d, bestfile))
    shutil.copy(fasta[bestidx], os.path.join(d, bestfile))
    ret_fa = os.path.join(d, bestfile)

    # Copy output files to output directory
    html = glob(os.path.join(d, '%s*global*.html' % rn))
    coords = [re.match('\S+global_visuAlign_(\d+)_(\d+).html', f) for f in html]
    coords = [list(map(int, m.groups())) if m else [0,0] for m in coords]
    sizes = [c[1]-c[0] for c in coords]
    if max(sizes) == 0: print("WARNING: Max sizes is 0.", file=sys.stderr)
    bestidx = sizes.index(max(sizes))
    bestfile = '%s.best_%d_%d.html' % (rn, coords[bestidx][0], coords[bestidx][1])
    if os.path.exists(os.path.join(d, bestfile)):
        print("WARNING: File %s exists." % os.path.join(d, bestfile), file=sys.stderr)
        os.unlink(os.path.join(d, bestfile))
    shutil.copy(html[bestidx], os.path.join(d, bestfile))
    ret_html = os.path.join(d, bestfile)
    
    return ret_fa, ret_html


def stageparser(parser):
    """ Add stage-specific options to argparse parser

    Args:
        parser (argparse.ArgumentParser): ArgumentParser object

    Returns:
        None

    """
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--fq1', type=sysutils.existing_file, required=True,
                        help='Fastq file with read 1')
    group1.add_argument('--fq2', type=sysutils.existing_file, required=True,
                        help='Fastq file with read 2')
    group1.add_argument('--ref_fa', type=sysutils.existing_file, required=True,
                        help='Reference sequence used to align reads (fasta)')
    group1.add_argument('--region_txt', type=sysutils.existing_file,
                        help='''File with regions to perform haplotype reconstruction. 
                                Regions should be specified using the samtools region 
                                specification format: RNAME[:STARTPOS[-ENDPOS]]''')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('PredictHaplo Options')
    group2.add_argument('--min_readlength', type=int, default=36,
                        help='''Minimum readlength passed to PredictHaplo''')

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

    parser.set_defaults(func=predict_haplo)


def predict_haplo(
        fq1=None, fq2=None, ref_fa=None, region_txt=None, outdir='.',
        min_readlength=36,
        keep_tmp=False, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to assemble haplotypes

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        ref_fa (str): Path to reference fasta file
        region_txt (str): Path to region file
        outdir (str): Path to output directory
        min_readlength (int): Minimum readlength passed to PredictHaplo
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        best_fa (list): Path to best haplotype files (FASTA)

    """
    # Check dependencies
    sysutils.check_dependency('PredictHaplo-Paired')
    sysutils.check_dependency('bwa')

    # Temporary directory
    tempdir = sysutils.create_tempdir('predict_haplo', None, quiet, logfile)

    # Load reference fasta
    refs = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}

    # Identify reconstruction regions
    regions = []
    if region_txt:
        sysutils.log_message('Found regions file.\n', quiet, logfile)  
        for l in open(region_txt, 'r'):
            rname, spos, epos = sequtils.region_to_tuple(l.strip())
            if rname not in refs:
                raise PipelineStepError("ERROR: reference %s not valid" % rname)
            spos = 1 if spos is None else spos
            epos = len(refs[rname]) if epos is None else epos
            regions.append(('PH%02d' % (len(regions)+1), rname, spos, epos))
    else:
        for rname, s in refs.items():
            regions.append(('PH%02d' % (len(regions)+1), rname, 1, len(s)))
    
    sysutils.log_message('[--- Haplotype Reconstruction Regions ---]\n', quiet, logfile)
    for iv in regions:
        sysutils.log_message('%s -- %s:%d-%d\n' % iv, quiet, logfile)
    
    # Create alignment for each REFERENCE in the reconstruction regions
    alnmap = {}
    for ph, rname, spos, epos in regions:
        if rname not in alnmap:
            # Create alignment
            tmp_ref_fa = os.path.join(tempdir, 'ref.%d.fa' % len(alnmap))
            tmp_sam = os.path.join(tempdir, 'aligned.%d.sam' % len(alnmap))
            SeqIO.write(refs[rname], tmp_ref_fa, 'fasta')
            cmd1 = ['bwa', 'index', tmp_ref_fa, ]
            cmd2 = ['bwa', 'mem', tmp_ref_fa, fq1, fq2, '|', 'samtools', 'view', '-h', '-F', '12', '>', tmp_sam, ]
            cmd3 = ['rm', '-f', '%s.*' % tmp_ref_fa]
            sysutils.command_runner(
                [cmd1, cmd2, cmd3], 'predict_haplo:setup', quiet, logfile, debug
            )
            alnmap[rname] = (tmp_ref_fa, tmp_sam)
    
    best_fa = []
    # Run PredictHaplo for each REGION
    for ph, rname, spos, epos in regions:
        msg = "Reconstruction region %s:" % ph
        msg += " %s:%d-%d\n" % (rname, spos, epos)
        sysutils.log_message(msg, quiet, logfile)

        # rename the predicthaplo number (PH##) to include region (now: PH##_reg)
        ph = '%s_%s' % (ph, rname.split('|')[-2])

        # Construct params specific for region
        reg_params = dict(DEFAULTS)
        reg_params['min_readlength'] = min_readlength        
        reg_params['reconstruction_start'] = spos
        reg_params['reconstruction_stop'] = epos
        reg_params['prefix'] = '%s_out.' % ph


        # Lookup reference and alignment filename 
        reg_params['ref_fasta'] = os.path.basename(alnmap[rname][0])
        reg_params['alignment'] = os.path.basename(alnmap[rname][1])

        # Create config file for region
        config_file = '%s.config' % ph
        with open(os.path.join(tempdir, config_file), 'w') as outh:
            tmpconfig = config_template % reg_params
            print(tmpconfig.replace('###', '%'), file=outh)
        try:
            # Run PredictHaplo
            cmd1 = ['cd', tempdir, ]
            cmd2 = ['PredictHaplo-Paired', config_file, '&>', '%s.log' % config_file]
            sysutils.command_runner(
                [cmd1, cmd2,], 'predict_haplo:%s' % ph, quiet, logfile, debug
            )

            if not debug:
                # Copy files
                dest = os.path.join(outdir, ph)
                if not os.path.exists(dest):
                    os.makedirs(dest)
                shutil.copy(os.path.join(tempdir, '%s.config.log' % ph), dest)
                for f in glob(os.path.join(tempdir, '%s_out*global*.fas' % ph)):
                    shutil.copy(f, dest)
                for f in glob(os.path.join(tempdir, '%s_out*global*.html' % ph)):
                    shutil.copy(f, dest)
                bf, bh = rename_best(dest, ph)
                best_fa.append((ph, bf))
        except PipelineStepError as e:
            print(e, file=sys.stderr)
            if e.returncode == 139:
                print("PredictHaplo segfaulted", file=sys.stderr)
            best_fa.append((ph, None))
    
    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'predict_haplo', quiet, logfile)

    return best_fa


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='''Assemble haplotypes with PredictHaplo.''',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()