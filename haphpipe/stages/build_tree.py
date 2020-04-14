# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import os
import argparse
import shutil
import random

from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import MissingRequiredArgument

__author__ = "Uzma Rentia and Margaret C. Steiner"
__copyright__ = "Copyright (C) 2020 Uzma Rentia and Margaret C. Steiner"

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

def stageparser(parser):
    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--seqs', type=sysutils.existing_file,
                        help='Input alignment in PHYLIP or FASTA format')
    group1.add_argument('--data_type', type=str, help="type of data: NUC, AA, BIN, or MULTI (default: NUC)")
    group1.add_argument('--output_name', type=str,
                        help='Run name for trees')
    group1.add_argument('--model', type=str,
                        help="substitution model")
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')

    group2 = parser.add_argument_group('RAxML Settings')
    ### note: some RAxML options have been omitted if they are not applicable and/or not included in the current
    ### version of RAxML from Bioconda
    group2.add_argument('--run_full_analysis', action='store_true', help="Run bootstrap search and find best ML tree")
    group2.add_argument('--outgroup', type=str,
                        help='outgroup for tree')
    group2.add_argument('--parsimony_seed', type=int,
                        help='Parsimony Random Seed')
    group2.add_argument('--wgtFile', type=sysutils.existing_file,
                        help='Column weight file name to assign individual weights to each column of the alignment')
    group2.add_argument('--secsub', type=str,
                        help='Specify secondary structure substitution models, must also include file defining the secondary structure')
    group2.add_argument('--bootstrap', type=int, help='bootstrapRandomNumberSeed for non-parametric bootstrapping')
    group2.add_argument('--bootstrap_threshold', type=float, help='threshold for bootstopping criteria')
    group2.add_argument('--numCat', type=int,
                        help='number of distinct rate categories for RAxML when model of rate heterogeneity is set to CAT')
    group2.add_argument('--rand_starting_tree', action='store_true',
                        help='ML optimization from random starting tree')
    group2.add_argument('--convergence_criterion', action='store_true',
                        help='ML search convergence criterion')
    group2.add_argument('--likelihoodEpsilon', type=float,
                        help='set model optimization precision in log likelihood units for final optimization of tree topology')
    group2.add_argument('--excludeFileName', type=sysutils.existing_file,
                        help='file contains specifications of alignment positions to be excluded')
    group2.add_argument('--algo_option', type=str,
                        help='select what kind of algorithm RAxML shall execute')
    group2.add_argument('--cat_model', action='store_true',
                        help='enable ML tree searches under CAT model for very large trees')
    group2.add_argument('--groupingFile', type=sysutils.existing_file,
                        help='file name of a multifurcating constraint tree')
    group2.add_argument('--placementThreshold', type=float,
                        help='threshold value for ML­based evolutionary placement algorithm heuristics')
    group2.add_argument('--disable_pattern_compression', action='store_true',
                        help='disable pattern compression')
    group2.add_argument('--InitialRearrangement', type=int,
                        help='radius for pruned sub-tree re-insertion')
    group2.add_argument('--posteriori', type=str,
                        help='posteriori bootstopping analysis')
    group2.add_argument('--print_intermediate_trees', action='store_true',
                        help='print out a couple of intermediate trees')
    group2.add_argument('--majorityrule', type=str,
                        help='Compute majority rule consensus tree')
    group2.add_argument('--print_branch_length', action='store_true',
                        help='bootstrapped trees should be printed with branch lengths')
    group2.add_argument('--ICTCmetrics', type=str,
                        help='compute the TC and IC metrics on a consensus tree')
    group2.add_argument('--partition_branch_length', action='store_true',
                        help='Switch on estimation of individual per­partition branch lengths')
    group2.add_argument('--disable_check', action='store_true',
                        help='Disable check for completely undetermined sequence in alignment')
    group2.add_argument('--AAmodel', type=sysutils.existing_file,
                        help='Specify the file name of a user­defined AA (Protein) substitution model')
    group2.add_argument('--multiplemodelFile', type=sysutils.existing_file,
                        help='Specify the file name which contains the assignment of models to alignment partitions for multiple models of substitution')
    group2.add_argument('--binarytree', type=sysutils.existing_file,
                        help='Specify the file name of a binary constraint tree')
    group2.add_argument('--BinaryParameterFile', type=sysutils.existing_file,
                        help='Specify the file name of a binary model parameter file that has previously been generated with RAxML using the ­f e tree evaluation option.')
    group2.add_argument('--SecondaryStructure', type=sysutils.existing_file,
                        help='Specify the name of a secondary structure file')
    group2.add_argument('--UserStartingTree', type=sysutils.existing_file,
                        help='Specifies a user starting tree file name which must be in Newick format')
    group2.add_argument('--median_GAMMA', action='store_true',
                        help='use the median for the discrete approximation of the GAMMA model of rateheterogeneity')
    group2.add_argument('--version_info', action='store_true',
                        help='Display version information')
    group2.add_argument('--rate_heterogeneity', action='store_true',
                        help='Disable rate heterogeneity among site model and use one without rate heterogeneity instead')
    group2.add_argument('--directory', type=str,
                        help='Full directory of output file')
    group2.add_argument('--window', type=int,
                        help='Sliding window size for leave­one­out site­specific placement bias algorithm')
    group2.add_argument('--RapidBootstrapNumSeed', type=int,
                        help='Specify an integer number (random seed) and turn on rapid bootstrapping')
    group2.add_argument('--random_addition', action='store_true',
                        help='RAxML will only do a randomized stepwise addition order parsimony tree reconstruction without performing any additional SPRs')
    group2.add_argument('--starting_tree', action='store_true',
                        help='Only for computing parsimony')
    group2.add_argument('--quartetGroupingFileName', type=sysutils.existing_file,
                        help='Pass a quartet grouping file name defining four groups from which to draw quartets')
    group2.add_argument('--multipleTreeFile', type=sysutils.existing_file,
                        help='Specify the file name of a file containing multiple trees e.g. from a bootstrap that shall be used to draw bipartition values onto a tree provided with ­t.')
    group2.add_argument('--NumberofRuns', type=int,
                        help='Specify the number of alternative runs on distinct starting trees')
    group2.add_argument('--mesquite', action='store_true',
                        help='Print output files that can be parsed by Mesquite')
    group2.add_argument('--silent', action='store_true',
                        help='Disables printout of warnings related to identical sequences and entirely undetermined sites in the alignment')
    group2.add_argument('--noseqcheck', action='store_true',
                        help='Disables checking the input MSA for identical sequences and entirely undetermined sites')
    group2.add_argument('--nobfgs', action='store_true',
                        help='Disables automatic usage of BFGS method to optimize GTR rates on unpartitioned DNA datasets')
    group2.add_argument('--epaPlaceNum', type=int,
                        help='specify the number of potential placements you want to keep for each read in the EPA algorithm')
    group2.add_argument('--epaProbThreshold', type=float,
                        help='specify a percent threshold for including potential placements of a read depending on the maximum placement weight for this read')
    group2.add_argument('--epaLikelihood', type=float,
                        help='Specify an accumulated likelihood weight threshold')
    group2.add_argument('--HKY85', action='store_true',
                        help='specify that all DNA partitions will evolve under the HKY85 model')
    group2.add_argument('--BootstrapPerm', type=str,
                        help='specify the number of permutations to be conducted for the bootstopping/bootstrap convergence test; minimum 100')
    group2.add_argument('--option_help', action='store_true',
                        help='Display Help')

    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Keep temporary directory')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')

    parser.set_defaults(func=build_tree)


def build_tree(seqs=None, data_type='NUC', run_full_analysis=None, output_name='build_tree.tre', outdir='.',
               treedir='hp_build_tree', model='GTRGAMMAIX',
               outgroup=None, parsimony_seed=1234,
               wgtFile=None, secsub=None, bootstrap=None, bootstrap_threshold=None, numCat=None,
               rand_starting_tree=None, convergence_criterion=None,
               likelihoodEpsilon=None, excludeFileName=None,
               algo_option=None, cat_model=None, groupingFile=None, placementThreshold=None,
               disable_pattern_compression=None,
               InitialRearrangement=None,
               posteriori=None, print_intermediate_trees=None,
               majorityrule=None, print_branch_length=None, ICTCmetrics=None, partition_branch_length=None,
               disable_check=None, AAmodel=None, multiplemodelFile=None,
               binarytree=None, BinaryParameterFile=None,
               SecondaryStructure=None, UserStartingTree=None, median_GAMMA=None, version_info=None,
               rate_heterogeneity=None,
               window=None, RapidBootstrapNumSeed=None,
               random_addition=None, starting_tree=None, quartetGroupingFileName=None, multipleTreeFile=None,
               NumberofRuns=None, mesquite=None,
               silent=None, noseqcheck=None,
               nobfgs=None,
               epaPlaceNum=None, epaProbThreshold=None,
               epaLikelihood=None,
               HKY85=None, BootstrapPerm=None,
               quiet=False, logfile=None, debug=False, keep_tmp=False, option_help=None):

    # Check dependencies
    sysutils.check_dependency('raxmlHPC')

    cmd1 = []

    # check for required input

    if option_help is True:
        cmd4 = ['raxmlHPC', '-h']
        sysutils.command_runner([cmd4], 'build_tree', quiet, logfile, debug)
        return

    if version_info is True:
        cmd5 = ['raxmlHPC', '-v']
        sysutils.command_runner([cmd5], 'build_tree', quiet, logfile, debug)

    if seqs is None and option_help is None:
        msg = 'No alignment provided'
        raise sysutils.PipelineStepError(msg)

    # check model compatibility
    if data_type is not 'AA' and 'PROT' in model:
        msg = 'Protein model given for non-amino acid data'
        raise sysutils.PipelineStepError(msg)
    if data_type is not 'MULTI' and 'MULTI' in model:
        msg = 'Multi-state model given for non-multi-state data'
        raise sysutils.PipelineStepError(msg)
    if data_type is not 'BIN' and 'BIN' in model:
        msg = 'Binary model given for non-binary data'
        raise sysutils.PipelineStepError(msg)
    if data_type is not 'NUC':
        if data_type not in model:
            msg = 'model and data type not compatible'
            raise sysutils.PipelineStepError(msg)

    # Set Output Directory
    output_dir = os.path.join(outdir, treedir)
    cmd0 = ['mkdir -p %s' % output_dir]

    sysutils.command_runner([cmd0], 'build_tree', quiet, logfile, debug)

    # Temporary directory
    tempdir = sysutils.create_tempdir('build_tree', None, quiet, logfile)

    if run_full_analysis is True:
        # generate seeds
        seed1 = random.randint(10000, 99999)
        seed2 = random.randint(10000, 99999)
        cmd1 = ['echo', 'Using parsimony seed %s and bootstrap seed %s' % (seed1, seed2)]
        sysutils.command_runner([cmd1], 'build_tree', quiet, logfile, debug)
        # run raxml
        cmd2 = ['raxmlHPC', '-w %s' % os.path.abspath(tempdir), '-f a', '-p %d' % seed1, '-x %d' % seed2, '-# 100',
                '-m %s' % model, '-s %s' % os.path.abspath(seqs), '-n %s' % output_name]
        sysutils.command_runner([cmd2], 'build_tree', quiet, logfile, debug)


    else:
        # start raxml command
        cmd1 = ['raxmlHPC', '-w %s' % os.path.abspath(tempdir), '-p %d' % parsimony_seed, '-m %s' % model]

        if outgroup is not None:
            cmd1 += ['-o', '%s' % outgroup]
        if wgtFile is not None:
            cmd1 += ['-a', '%s' % os.path.join('.', wgtFile)]
        if secsub is not None and SecondaryStructure is not None:
            cmd1 += ['-A', '%s' % secsub]
            cmd1 += ['-S', '%s' % os.path.join('.', SecondaryStructure)]
        elif secsub is not None and SecondaryStructure is None:
            msg = 'Need to specify a file defining the secondary structure via the ­S option'
            raise sysutils.PipelineStepError(msg)
        if bootstrap is not None:
            cmd1 += ['-b', '%d' % bootstrap]
        if bootstrap_threshold is not None:
            cmd1 += ['-B', '%f' % bootstrap_threshold]
        if numCat is not None:
            cmd1 += ['-c', '%d' % numCat]
        if rand_starting_tree is True:
            cmd1 += ['-d']
        if convergence_criterion is True:
            cmd1 += ['-D']
        if likelihoodEpsilon is not None:
            cmd1 += ['-e', '%f' % likelihoodEpsilon]
        if excludeFileName is not None:
            cmd1 += ['-E', '%s' % os.path.join('.', excludeFileName)]
        if algo_option is not None and algo_option in ['a', 'A', 'b', 'B', 'c', 'C', 'd', 'D', 'e', 'E', 'F', 'g', 'G',
                                                       'h',
                                                       'H', 'i', 'I', 'j', 'J', 'k', 'm', 'n', 'N', 'o', 'p', 'q', 'r',
                                                       'R',
                                                       's',
                                                       'S', 't', 'T', 'U', 'v', 'V', 'w', 'W', 'x',
                                                       'y']:
            cmd1 += ['-f', '%s' % algo_option]
        if cat_model is True:
            cmd1 += ['-F']
        if groupingFile is not None:
            cmd1 += ['-g', '%s' % os.path.join('.', groupingFile)]
        if placementThreshold is not None:
            cmd1 += ['-G', '%f' % placementThreshold]
        if disable_pattern_compression is True:
            cmd1 += ['-H']
        if InitialRearrangement is not None:
            cmd1 += ['-i', '%d' % InitialRearrangement]
        if posteriori is not None:
            cmd1 += ['-I', '%s' % posteriori]
        if print_intermediate_trees is True:
            cmd1 += ['-j']
        if (majorityrule is not None) and (multipleTreeFile is not None):
            cmd1 += ['-J', '%s' % majorityrule]
            cmd1 += ['-z', '%s' % os.path.join('.', multipleTreeFile)]
        elif majorityrule is not None and multipleTreeFile is None:
            msg = 'Need to provide a tree file containing several UNROOTED trees via the ­z option'
            raise sysutils.PipelineStepError(msg)
        if print_branch_length is True:
            cmd1 += ['-k']
        if ICTCmetrics is not None:
            cmd1 += ['-L', '%s' % ICTCmetrics]
        if partition_branch_length is True:
            cmd1 += ['-M']
        if disable_check is True:
            cmd1 += ['-O']
        if AAmodel is not None:
            cmd1 += ['-P', '%s' % os.path.join('.', AAmodel)]
        if multiplemodelFile is not None:
            cmd1 += ['-q', '%s' % os.path.join('.', multiplemodelFile)]
        if binarytree is not None:
            cmd1 += ['-r', '%s' % os.path.join('.', binarytree)]
        if BinaryParameterFile is not None:
            cmd1 += ['-R', '%s' % os.path.join('.', BinaryParameterFile)]
        if SecondaryStructure is not None:
            cmd1 += ['-S', '%s' % os.path.join('.', SecondaryStructure)]
        if UserStartingTree is not None:
            cmd1 += ['-t', '%s' % os.path.join('.', UserStartingTree)]
        if median_GAMMA is True:
            cmd1 += ['-u']
        if rate_heterogeneity is True:
            cmd1 += ['-V']
        if window is not None:
            cmd1 += ['-W', '%d' % window]
        if RapidBootstrapNumSeed is not None:
            cmd1 += ['-x', '%d' % RapidBootstrapNumSeed]
        if random_addition is True:
            cmd1 += ['-X']
        if starting_tree is True:
            cmd1 += ['-y']
        if quartetGroupingFileName is not None:
            cmd1 += ['-Y', '%s' % os.path.join('.', quartetGroupingFileName)]
        if multipleTreeFile is not None:
            cmd1 += ['-z', '%s' % os.path.join('.', multipleTreeFile)]
        if NumberofRuns is not None:
            cmd1 += ['-N', '%d' % NumberofRuns]
        if mesquite is True:
            cmd1 += ['--mesquite']
        if silent is True:
            cmd1 += ['--silent']
        if noseqcheck is True:
            cmd1 += ['--no-seq-check']
        if nobfgs is True:
            cmd1 += ['--no-bfgs']
        if epaPlaceNum is not None:
            cmd1 += ['­­epa­keep­placements=%d' % epaPlaceNum]
        if epaProbThreshold is not None:
            cmd1 += ['­­epa­prob­threshold=%f' % epaProbThreshold]
        if epaLikelihood is not None:
            cmd1 += ['­­epa­accumulated­threshold=%f' % epaLikelihood]
        if HKY85 is True:
            cmd1 += ['--HKY85']
        if BootstrapPerm is not None:
            cmd1 += ['[­­bootstop­perms=%s' % BootstrapPerm]
        if option_help is True:
            cmd1 += ['-h']

        cmd1 += ['-s', '%s' % os.path.abspath(seqs), '-n',
                 '%s' % output_name]

        sysutils.command_runner([cmd1, ], 'build_tree', quiet, logfile, debug)

    # copy files from tmpdir to output directory
    if os.path.exists(os.path.join(tempdir, 'RAxML_bestTree.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_bestTree.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_info.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_info.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_perSiteLLs.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_perSiteLLs.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_bipartitionFrequencies.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_bipartitionFrequencies.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_bipartitionsBranchLabels.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_bipartitionsBranchLabels.%s' % output_name),
                    os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_bipartitions.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_bipartitions.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_bootstrap.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_bootstrap.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_checkpoint.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_checkpoint.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_randomTree.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_randomTree.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_parsimonyTree.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_parsimonyTree.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_result.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_result.%s' % output_name), os.path.abspath(output_dir))
    if os.path.exists(os.path.join(tempdir, 'RAxML_log.%s' % output_name)):
        shutil.copy(os.path.join(tempdir, 'RAxML_log.%s' % output_name), os.path.abspath(output_dir))

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'build_tree', quiet, logfile)

    cmd3 = ['echo', 'Stage completed. Output files are located here: %s\n' % os.path.abspath(output_dir)]
    sysutils.command_runner([cmd3, ], 'build_tree', quiet, logfile, debug)


def console():
    parser = argparse.ArgumentParser(
        description='Create Phylogenetic Tree with RAxML',
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
