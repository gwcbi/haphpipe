# -*- coding: utf-8 -*-

from __future__ import print_function
from __future__ import absolute_import
from builtins import map
import sys
import os
import re
import argparse
from subprocess import check_call, CalledProcessError
from subprocess import Popen, PIPE
from glob import glob
import shutil
import time
from collections import defaultdict

from Bio import SeqIO

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils

# import PipelineStepError, check_dependency
# from ..utils.sysutils import existing_file, existing_dir, command_runner, args_params
# from ..utils.sysutils import create_tempdir, remove_tempdir
#from ..utils.sequtils import fastagen, unambig_intervals
# from . import post_assembly #import samtools_depth, get_covered_intervals

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
    
    return


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
    group1.add_argument('--interval_txt', type=sysutils.existing_file,
                        help='''File with intervals to perform haplotype
                                reconstruction''')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('PredictHaplo Options')
    group2.add_argument('--min_depth', type=int,
                        help='''Minimum depth to consider position covered.
                                Ignored if intervals are provided''')
    group2.add_argument('--max_ambig', type=int, default=200,
                        help='''Maximum size of ambiguous sequence within a
                                reconstruction region. Ignored if intervals are
                                provided.''')
    group2.add_argument('--min_interval', type=int, default=200,
                        help='Minimum size of reconstruction interval')
    group2.add_argument('--min_readlength', type=int, default=36,
                        help='''Minimum readlength passed to PredictHaplo''')

    group3 = parser.add_argument_group('Settings')
    # group3.add_argument('--ncpu', type=int,
    #                    help='Number of CPU to use')
    # group3.add_argument('--max_memory', type=int,
    #                     help='Maximum memory to use (in GB)')
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
        fq1=None, fq2=None, ref_fa=None, interval_txt=None, outdir='.',
        min_depth=None, max_ambig=200, min_interval=200, min_readlength=36,
        # ncpu=1,
        keep_tmp=False, quiet=False, logfile=None, debug=False,
    ):
    """ Pipeline step to assemble haplotypes

    Args:
        fq1 (str): Path to fastq file with read 1
        fq2 (str): Path to fastq file with read 2
        ref_fa (str): Path to reference fasta file
        interval_txt (str): Path to interval file
        outdir (str): Path to output directory
        min_depth (int): Minimum depth to consider position covered
        max_ambig (int): Maximum size of ambiguous sequence within a region
        min_interval (int): Minimum size of reconstruction region
        min_readlength (int): Minimum readlength passed to PredictHaplo
        ncpu (int): Number of CPUs to use
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out_aligned (str): Path to aligned BAM file
        out_bt2 (str): Path to bowtie2 report

    """
    """ Assemble haplotypes with predicthaplo
    """
    # Check dependencies
    sysutils.check_dependency('PredictHaplo-Paired')
    # sysutils.check_dependency('samtools')
    sysutils.check_dependency('bwa')

    # try:
    #     x = check_call('samtools 2>&1 >/dev/null | grep -q "collate"', shell=True)
    #     if debug: print("Using samtools collate", file=sys.stderr)
    #     has_collate = True
    # except CalledProcessError as e:
    #     has_collate = False
    #     if debug: print("Using samtools sort", file=sys.stderr)
    
    # Temporary directory
    tempdir = sysutils.create_tempdir('predict_haplo', None, quiet, logfile)

    # Create alignment
    tmp_ref_fa = os.path.join(tempdir, 'ref.fa')
    tmp_sam = os.path.join(tempdir, 'aligned.sam')
    shutil.copy(ref_fa, tmp_ref_fa)
    cmd1 = ['bwa', 'index', tmp_ref_fa, ]
    cmd2 = ['bwa', 'mem', tmp_ref_fa, fq1, fq2, '>', tmp_sam, ]
    sysutils.command_runner(
        [cmd1, cmd2, ], 'predict_haplo:setup', quiet, logfile, debug
    )

    # Set up parameters
    ph_params = dict(DEFAULTS)
    ph_params['min_readlength'] = min_readlength
    ph_params['alignment'] = 'aligned.sam'

    # Load reference fasta
    refs = {s.id:s for s in SeqIO.parse(ref_fa, 'fasta')}

    # Load intervals
    recon_intervals = []
    if interval_txt:
        print('Found intervals file...', file=sys.stderr)
        for l in open(interval_txt, 'rU'):
            chrom = l.split(':')[0]
            s, e = tuple(map(int, l.split(':')[1].split('-')))
            recon_intervals.append((chrom, s, e))
    else:
        for sid, s in refs.items():
            recon_intervals.append((sid, 1, len(s)))

    print('Haplotype Reconstruction Regions:', file=sys.stderr)
    for iv in recon_intervals:
        print('%s:%d-%d' % iv, file=sys.stderr)

    # Run in parallel when the number of runs is less than number of CPU
    # Since I don't want to handle queueing right now.
    # do_parallel = 1 < len(recon_intervals) < ncpu

    # Setup runs
    runnames = []
    for i, iv in enumerate(recon_intervals):
        runname = 'PH%02d' % (i+1)
        runnames.append(runname)
        msg = "Reconstruction region %s:" % runname
        msg += " %s:%d-%d" % (iv[0], iv[1], iv[2])
        sysutils.log_message(msg, quiet, logfile)

        # Construct params specific for region
        reg_params = dict(ph_params)
        reg_params['reconstruction_start'] = iv[1]
        reg_params['reconstruction_stop'] = iv[2]
        reg_params['prefix'] = '%s_out.' % runname

        # Create single reference FASTA
        _ref_fa = '%s_ref.fasta' % runname
        SeqIO.write(refs[iv[0]], os.path.join(tempdir, _ref_fa), 'fasta')
        reg_params['ref_fasta'] = _ref_fa

        # Create config file for region
        config_file = '%s.config' % runname
        with open(os.path.join(tempdir, config_file), 'w') as outh:
            tmpconfig = config_template % reg_params
            print(tmpconfig.replace('###', '%'), file=outh)

        # Create the output directory for this run
        rundir = os.path.join(outdir, runname)
        if not os.path.exists(rundir):
            os.makedirs(rundir)

        # Name for log file
        logfile = os.path.join(rundir, '%s.log' % runname)

        # Commands (to be run within temporary directory
        cmds = [
            ['cd', tempdir, ],
            ['PredictHaplo-Paired', config_file, '&>', logfile, ],
            ['cp', '%s*global*.fas' % runname, rundir, ],
            ['cp', '%s*global*.html' % runname, rundir, ],
        ]
        sysutils.command_runner(
            cmds, 'predict_haplo', quiet, logfile, debug
        )

    if not debug:
        for rn in runnames:
            rename_best(os.path.join(outdir, rn), rn)

    #if not keep_tmp:
    #    sysutils.remove_tempdir(tempdir, 'predict_haplo', quiet, logfile)

    return


"""

    # Run in parallel when the number of runs is less than number of CPU
    # Since I don't want to handle queueing right now.
    do_parallel = 1 < len(recon_intervals) < ncpu
    
    if do_parallel:
        print("Processing in parallel...", file=sys.stderr)
        processes = []
    else:
        print("Processing in order...", file=sys.stderr)
        cmds = [['cd', tempdir, ], ]
    



        
        # Construct params specific for region
        reg_params = dict(ph_params)
        reg_params['reconstruction_start'] = iv[1]
        reg_params['reconstruction_stop'] = iv[2]
        reg_params['prefix'] = '%s.' % runnames[-1]
        
        # Create config file for region
        config_file = '%s.config' % runnames[-1]        
        with open(os.path.join(tempdir, config_file), 'w') as outh:
            tmpconfig = config_template % reg_params
            print(tmpconfig.replace('###', '%'), file=outh)
        
        # Create the output directory for this run
        rundir = os.path.join(outdir, runnames[-1])
        if not os.path.exists(rundir):
            os.makedirs(rundir)
        rundir = os.path.abspath(rundir)
        
        # Name for log file    
        logfile = os.path.join(rundir, '%s.log' % runnames[-1])
        
        # Commands (to be run within temporary directory
        rcmds = [
            ['PredictHaplo-Paired', config_file, '&>', logfile, ],
            ['cp', '%s*global*.fas' % runnames[-1], rundir, ],
            ['cp', '%s*global*.html' %  runnames[-1], rundir, ],
        ]
        if do_parallel:
            print("Spawning process for %s" % runnames[-1], file=sys.stderr)
            rcmds = [['cd', tempdir, ]] + rcmds
            cmdstr = ' && '.join(' '.join(c) for c in rcmds)
            print(cmdstr, file=sys.stderr)
            if not debug:
                processes.append((runnames[-1], Popen(cmdstr, shell=True)))
        else:
            cmds.extend(rcmds)
    
    if do_parallel:
        while processes:
            for tup in processes[:]:
                rn, p = tup
                if p.poll() is not None:
                    print("PredictHaplo for %s is complete" % rn, file=sys.stderr)
                    rename_best(os.path.join(outdir, rn), rn)
                    processes.remove(tup)
            time.sleep(2)
    else:
        sysutils.command_runner(
            cmds, 'predict_haplo', quiet, logfile, debug
        )
        if not debug:
            for rn in runnames:
                rename_best(os.path.join(outdir, rn), rn)

    print('tempdir: %s' % tempdir)
    #if not keep_tmp:
    #    sysutils.remove_tempdir(tempdir, 'predict_haplo', quiet, logfile)
    
    return
"""

"""
    # Copy output files to output directory
    for rn in runnames:
        rundir = os.path.join(outdir, rn)
        if not os.path.exists(rundir):
            os.makedirs(rundir)
        # Copy the fasta files
        fasta = glob(os.path.join(tempdir, '%s*global*.fas' % rn))
        coords = [re.match('\S+global_(\d+)_(\d+).fas', f) for f in fasta]
        coords = [map(int, m.groups()) if m else [0,0] for m in coords]
        sizes = [c[1]-c[0] for c in coords]
        for i,f in enumerate(fasta):
            shutil.copy(f, rundir)
            if sizes[i] == max(sizes):
                bestfile = '%s.best_%d_%d.fas' % (rn, coords[i][0], coords[i][1])
                shutil.copy(f, os.path.join(rundir, bestfile))
        # Copy the html files
        html = glob(os.path.join(tempdir, '%s*global*.html' % rn))
        coords = [re.match('\S+global_visuAlign_(\d+)_(\d+).html', f) for f in html]
        coords = [map(int, m.groups()) if m else [0,0] for m in coords]
        sizes = [c[1]-c[0] for c in coords]
        for i,f in enumerate(html):
            shutil.copy(f, rundir)
            if sizes[i] == max(sizes):
                bestfile = '%s.best_%d_%d.html' % (rn, coords[i][0], coords[i][1])
                shutil.copy(f, os.path.join(rundir, bestfile))
        # Copy the log
        if os.path.isfile(os.path.join(tempdir, '%s.log' % rn)):
            shutil.copy(os.path.join(tempdir, '%s.log' % rn), rundir)
"""    


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='''Assemble haplotypes with PredictHaplo''',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()