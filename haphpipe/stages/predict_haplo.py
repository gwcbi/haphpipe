#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import re
import argparse
from subprocess import check_call, CalledProcessError
from subprocess import Popen, PIPE
from glob import glob
import shutil

from ..utils.sysutils import PipelineStepError, check_dependency
from ..utils.sysutils import existing_file, existing_dir, command_runner, args_params
from ..utils.sysutils import create_tempdir, remove_tempdir
from ..utils.sequtils import fastagen, unambig_intervals

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


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


DEFAULTS = { 'prefix': 'PH',
             'do_visualize': 1,
             'have_true_haplotypes': 0,
             'file_true_haplotypes': 'dummy.fasta',
             'do_local_analysis': 1,
             'max_reads': 10000,
             'entropy_threshold': '4e-2',
             'reconstruction_start': 6000,
             'reconstruction_stop': 8000,
             'min_mapping_qual': 25,
             'min_readlength': 10,
             'max_gap_fraction': 0.05,
             'min_align_score_fraction': 0.35,
             'alpha_MN_local': 25,
             'min_overlap_factor': 0.85,
             'local_window_size_factor': 0.7,
             'max_number_of_clusters': 25,
             'mcmc_iterations': 501,
             'include_deletions': 1,
}


def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--alignment', type=existing_file,
                        help='Alignment file (BAM)')
    group1.add_argument('--ref_fa', type=existing_file,
                        help='Reference sequence used to align reads (fasta)')
    group1.add_argument('--outdir', type=existing_dir,
                        help='Output directory')
    
    group2 = parser.add_argument_group('Predict Haplo Options')
    group2.add_argument('--min_interval', type=int,
                        help='Minimum size of reconstruction interval')
    group2.add_argument('--max_ambig', type=int,
                        help='''Maximum size of ambiguous sequence within a reconstruction
                                region''')
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--ncpu', type=int,
                        help='Number of CPU to use')
    group3.add_argument('--max_memory', type=int,
                        help='Maximum memory to use (in GB)')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Additional options')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    parser.set_defaults(func=predict_haplo)


def predict_haplo(alignment=None, ref_fa=None, outdir='.',
        min_interval=200, max_ambig=200,
        ncpu=1, max_memory=50, keep_tmp=False, debug=False,
    ):
    """ Assemble haplotypes with predicthaplo
    """
    # Check dependencies
    check_dependency('PredictHaplo-Paired')
    check_dependency('samtools')
    try:
        x = check_call('samtools 2>&1 >/dev/null | grep -q "collate"', shell=True)
        if debug: print >>sys.stderr, "Using samtools collate"
        has_collate = True
    except CalledProcessError as e:
        has_collate = False
        if debug: print >>sys.stderr, "Using samtools sort"        

    # Temporary directory
    tempdir = create_tempdir('predict_haplo')
    
    # Set up parameters        
    params = dict(DEFAULTS)
    
    # Load reference fasta
    with open(ref_fa, 'rU') as fh:
        seqs = [(n,s) for n,s in fastagen(fh)]
    
    assert len(seqs) == 1, 'ERROR: Reference must contain exactly one sequence'
    name, seq = seqs[0]
    
    # Identify reconstruction intervals
    recon_intervals = [iv for iv in unambig_intervals(seq, max_ambig)]
    recon_intervals = [iv for iv in recon_intervals if iv[1]-iv[0] > min_interval]
    if not recon_intervals:
        print >>sys.stderr, "No intervals larger than %d were found" % min_interval
        sys.exit()
    
    # Copy reference fasta
    with open(os.path.join(tempdir, 'reference.fasta'), 'w') as outh:
        print >>outh, '>%s\n%s' % (name, seq)
    params['ref_fasta'] = 'reference.fasta'
    
    # Collate or sort 
    if has_collate:
        cmd1a = ['samtools', 'collate',
            alignment, 
            os.path.join(tempdir, 'collated'),
        ]
    else:
        cmd1a = ['samtools', 'sort',
            '-n',
            alignment, 
            os.path.join(tempdir, 'collated'),
        ]
    cmd1b = ['samtools', 'view',
        '-h', os.path.join(tempdir, 'collated.bam'), 
        '>', 
        os.path.join(tempdir, 'collated.sam'),
    ]
    cmd1c = ['rm', '-f', os.path.join(tempdir, 'collated.bam')]

    # Create the SAM file
    command_runner([cmd1a, cmd1b, cmd1c ], 'predict_haplo:setup', debug)
    params['alignment'] = 'collated.sam'
    
    # Run in parallel when the number of runs is less than number of CPU
    # Since I don't want to handle queueing right now.
    do_parallel = 1 < len(recon_intervals) < ncpu
    
    if do_parallel:
        print >>sys.stderr, "Processing in parallel..."
        processes = []
    else:
        print >>sys.stderr, "Processing in order..."
        cmds = [['cd', tempdir, ], ]
    
    runnames = []
    for i,iv in enumerate(recon_intervals):
        runnames.append('PH%02d' % (i+1))
        print >>sys.stderr, "Reconstruction region %s: %d - %d" % (runnames[-1], iv[0], iv[1])
        config_file = '%s.config' % runnames[-1]
        tparams = dict(params)
        tparams['reconstruction_start'] = iv[0]
        tparams['reconstruction_stop'] = iv[1]
        tparams['prefix'] = '%s.' % runnames[-1]
        with open(os.path.join(tempdir, config_file), 'w') as outh:
            tmpconfig = config_template % tparams
            print >>outh, tmpconfig.replace('###', '%')
        
        if do_parallel:
            print >>sys.stderr, "spawning process for %s" % runnames[-1]
            processes.append((runnames[-1],
                Popen('cd %s && PredictHaplo-Paired %s' % (tempdir, config_file),
                      shell=True, stdout=PIPE)
            ))
        else:
            cmds.append(['PredictHaplo-Paired', config_file, ])
    
    if do_parallel:
        while processes:
            for tup in processes[:]:
                rn, p = tup
                if p.poll() is not None:
                    print >>sys.stderr, "PredictHaplo for %s is complete" % rn
                    with open(os.path.join(tempdir, '%s.log' % rn), 'w') as outh:
                        outh.write(p.stdout.read())
                    p.stdout.close()
                    processes.remove(tup)
    else:
        command_runner(cmds, 'predict_haplo', debug)
    
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
    
    if not keep_tmp:
        remove_tempdir(tempdir, 'predict_haplo')
    
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assembly haplotypes with PredictHaplo')
    stageparser(parser)
    args = parser.parse_args()
    args.func(**args_params(args))
