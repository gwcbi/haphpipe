from __future__ import print_function
import os
import sys
import argparse
import shutil
from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import MissingRequiredArgument

__author__ = 'Margaret C. Steiner'
__copyright__ = 'Copyright (C) 2020 Margaret C. Steiner'


def stageparser(parser):
    ### input options with argparse ###

    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--seqs', type=sysutils.existing_file, help='Alignment in FASTA or PHYLIP format')
    group1.add_argument('--run_id',type=str,help="Prefix for output files")
    group1.add_argument('--outname',type=str,help='Name for output file (Default: modeltest_results)')
    group1.add_argument('--outdir',type=sysutils.existing_dir,help='Output directory (Default: .)')

    group2 = parser.add_argument_group('ModelTest-NG Options')
    group2.add_argument('--data_type',type=str,help='Data type: nt or aa')
    group2.add_argument('--partitions',type=sysutils.existing_file,help='Partitions file')
    group2.add_argument('--seed',type=int,help='Seed for random number generator')
    group2.add_argument('--topology',type=str,help='Starting topology: ml, mp, fixed-ml-jc, fixed-ml-gtr, fixed-mp, random, or user')
    group2.add_argument('--utree',type=sysutils.existing_file,help='User-defined starting tree')
    group2.add_argument('--force',action='store_true',help='force output overriding')
    group2.add_argument('--asc_bias',type=str,help='Ascertainment bias correction: lewis, felsenstein, or stamatakis')
    group2.add_argument('--frequencies',type=str,help='Candidate model frequencies: e (estimated) or f (fixed)')
    group2.add_argument('--het',type=str,help='Set rate heterogeneity: u (uniform), i (invariant sites +I), g (gamma +G), or f (both invariant sites and gamma +I+G)')
    group2.add_argument('--models',type=sysutils.existing_file,help='Text file with candidate models, one per line')
    group2.add_argument('--schemes',type=int,help='Number of predefined DNA substitution schemes evaluated: 3, 5, 7, 11, or 203')
    group2.add_argument('--template',type=str,help='Set candidate models according to a specified tool: raxml, phyml, mrbayes, or paup')

    group3 = parser.add_argument_group('Options')
    group3.add_argument('--ncpu', type=int, default=1, help='Number of CPU to use')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                    (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Name for log file (output)')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    group3.add_argument('--keep_tmp',action='store_true',help='Keep temporary directory')

    parser.set_defaults(func=model_test)


def model_test(seqs=None,outname='modeltest_results',run_id=None, data_type='nt',partitions=None,seed=None,topology='ml',
               utree=None,force=None,asc_bias=None,frequencies=None,het=None,models=None,schemes=None,template=None,
               ncpu=1,quiet=False,logfile=None,debug=False,outdir='.',keep_tmp=False):

    # check dependency
    sysutils.check_dependency('modeltest-ng')

    # check required input & input options
    if seqs is None:
        msg="No alignment given"
        raise sysutils.MissingRequiredArgument(msg)
    if data_type not in ['nt','aa']:
        raise sysutils.PipelineStepError("Data type not valid")
    if topology not in ['ml','mp','fixed-ml-jc','fixed-ml-gtr','fixed-mp','random','user']:
        raise sysutils.PipelineStepError("Topology not valid")

    # make tempdir
    tempdir = sysutils.create_tempdir('model_test', None, quiet, logfile)

    # add prefix
    if run_id is not None:
        outname = run_id+'_'+outname

    # build command
    cmd1 = ['modeltest-ng -i %s' % seqs, '-t %s' % topology,'-o %s' % os.path.join(tempdir,outname),
            '-p %d' % ncpu, '-d %s' % data_type]

    if partitions is not None:
        cmd1 += ['-q %s' % partitions]

    if seed is not None:
        cmd1 += ['-r %d' % seed]

    if utree is not None:
        cmd1 += ['-u %s' % utree]

    if force is True:
        cmd1 += ['--force']

    if asc_bias is not None and asc_bias in ['lewis','felsenstein','stamatakis']:
        cmd1 += ['-a %s' % asc_bias]
    elif asc_bias is not None:
        raise sysutils.PipelineStepError("ASC bias correction not valid")

    if frequencies is not None and frequencies in ['e','f']:
        cmd1 += ['-f %s' % frequencies]
    elif frequencies is not None:
        raise sysutils.PipelineStepError("Frequencies not valid")

    if het is not None and het in ['u','i','g','f']:
        cmd1 += ['-h %s' % het]
    elif het is not None:
        raise sysutils.PipelineStepError("Rate heterogeneity not valid")

    if models is not None:
        with open(models,'r') as f:
            model_list = f.read().splitlines()
        for m in model_list:
            if data_type == 'nt' and m not in ['JC','HKY','TrN','TPM1','TPM2','TPM3','TIM1','TIM2','TIM3','TVM','GTR']:
                raise sysutils.PipelineStepError("At least one model is not valid")
            elif data_type == 'aa' and m not in ['DAYHOFF','LG','DCMUT','JTT','MTREV','WAG','RTREV','CPREV','VT','BLOSUM62',
                                                 'MTMAM','MTART','MTZOA','PMB','HIVB','HIVW','JTTDCMUT','FLU','SMTREV']:
                raise sysutils.PipelineStepError("At least one model is not valid")
        cmd1 += ['-m %s' % str(model_list)[1:-1]]

    if schemes is not None and schemes in [3,5,7,11,203]:
        cmd1 += ['-s %d' % schemes]
    elif schemes is not None:
        raise sysutils.PipelineStepError("Schemes not valid")

    if template is not None and template in ['raxml','phyml','mrbayes','paup']:
        cmd1 += ['-T %s' % template]
    elif template is not None:
        raise sysutils.PipelineStepError("Template not valid")

    # run command
    try:
        sysutils.command_runner([cmd1, ], 'model_test', quiet, logfile, debug)
    except sysutils.PipelineStepError as p:
        if p.returncode == -6:
            print("Warning: ignoring returncode -6")
        else:
            raise sysutils.PipelineStepError("Error in ModelTest-NG")

    # copy output file and delete tempdir
    if os.path.exists(os.path.join(tempdir, '%s.out' % outname)):
        shutil.copy(os.path.join(tempdir, '%s.out' % outname), os.path.abspath(outdir))
    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'model_test', quiet, logfile)

    # Parse .out file and write TSV summary file
    criteria = []
    bestmods = []
    with open(os.path.join(outdir,'%s.out' % outname)) as f1:
        for line in f1.read().splitlines():
            if "Best model according to" in line:
                criteria += line.split(' ')[-1:]
            if "Model: " in line:
                bestmods += line.split(' ')[-1:]
    with open(os.path.join(outdir,'%s_summary.tsv' % outname),'w') as f2:
        f2.write('File\tCriteria\tBest Model\n')
        for i in range(len(criteria)):
            f2.write('%s\t%s\t%s\n' % (seqs,criteria[i],bestmods[i]))

    # completion message
    cmd2 = ['echo', 'Stage completed. Output file is located here: %s\n' % os.path.abspath(os.path.join(outdir,'%s.out' % outname)),
            'echo','Summary TSV file is located here: %s\n' % os.path.abspath(os.path.join(outdir,'%s_summary.tsv' % outname))]
    sysutils.command_runner([cmd2,], 'model_test', quiet, logfile, debug)

    return


def console():
    parser = argparse.ArgumentParser(
        description='Determine best-fit evolutionary model with ModelTest-NG.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    try:
        args.func(**sysutils.args_params(args))
    except MissingRequiredArgument as e:
        parser.print_usage()
        print('error: %s' % e, file=sys.stderr)


if __name__ == '__main__':
    console()
