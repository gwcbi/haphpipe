from __future__ import print_function
import os
import sys
import argparse
from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import MissingRequiredArgument
from Bio import SeqIO

__author__ = 'Margaret C. Steiner'
__copyright__ = 'Copyright (C) 2020 Margaret C. Steiner'


def stageparser(parser):
    ### input options with argparse ###

    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--seqs', type=sysutils.existing_file, help='FASTA file with sequences to be aligned')
    group1.add_argument('--dir_list', type=sysutils.existing_file,
                        help='List of directories which include either a final.fna or ph_haplotypes.fna file, one on each line')
    group1.add_argument('--haplotypes',action='store_true',help='files are haplotype files (default: False)')
    group1.add_argument('--ref_gtf', type=sysutils.existing_file, help='Reference GTF file')
    group1.add_argument('--out_align', help='Name for alignment file')
    group1.add_argument('--nuc', action='store_true', help='Assume nucleotide')
    group1.add_argument('--amino', action='store_true', help='Assume amino')
    group1.add_argument('--clustalout', action='store_true', help='Clustal output format')
    group1.add_argument('--phylipout', action='store_true', help='PHYLIP output format')
    group1.add_argument('--inputorder', action='store_true', help='Output order same as input')
    group1.add_argument('--reorder', action='store_true', help='Output order aligned')
    group1.add_argument('--treeout', action='store_true', help='Guide tree is output to the input.tree file')
    group1.add_argument('--quiet_mafft', action='store_true', help='Do not report progress')
    group1.add_argument('--outdir', type=sysutils.existing_dir, help='Output directory')

    group2 = parser.add_argument_group('MAFFT Algorithm Options')
    group2.add_argument('--algo',
                        help='Use different algorithm in command: linsi, ginsi, einsi, fftnsi, fftns, nwns, nwnsi')
    group2.add_argument('--auto', action='store_true', help='Automatically select algorithm')
    group2.add_argument('--sixmerpair', action='store_true',
                        help='Calculate distance based on shared 6mers, on by default')
    group2.add_argument('--globalpair', action='store_true', help='Use Needleman-Wunsch algorithm')
    group2.add_argument('--localpair', action='store_true', help='Use Smith-Waterman algorithm')
    group2.add_argument('--genafpair', action='store_true', help='Use local algorithm with generalized affine gap cost')
    group2.add_argument('--fastapair', action='store_true', help='Use FASTA for pairwise alignment')
    group2.add_argument('--weighti', type=float, help='Weighting factor for consistency term')
    group2.add_argument('--retree', type=int, help='Number of times to build guide tree')
    group2.add_argument('--maxiterate', type=int, help='Number of cycles for iterative refinement')
    group2.add_argument('--noscore', action='store_true', help='Do not check alignment score in iterative alignment')
    group2.add_argument('--memsave', action='store_true', help='Use Myers-Miller algorithm')
    group2.add_argument('--parttree', action='store_true', help='Use fast tree-building method with 6mer distance')
    group2.add_argument('--dpparttree', action='store_true', help='Use PartTree algorithm with distances based on DP')
    group2.add_argument('--fastaparttree', action='store_true',
                        help='Use PartTree algorithm with distances based on FASTA')
    group2.add_argument('--partsize', type=int, help='Number of partitions for PartTree')
    group2.add_argument('--groupsize', type=int, help='Max number of sequences for PartTree')

    group3 = parser.add_argument_group('MAFFT Parameters')
    group3.add_argument('--lop', type=float, help='Gap opening penalty')
    group3.add_argument('--lep', type=float, help='Offset value')
    group3.add_argument('--lexp', type=float, help='Gap extension penalty')
    group3.add_argument('--LOP', type=float, help='Gap opening penalty to skip alignment')
    group3.add_argument('--LEXP', type=float, help='Gap extension penalty to skip alignment')
    group3.add_argument('--bl', type=int, help='BLOSUM matrix: 30, 45, 62, or 80')
    group3.add_argument('--jtt', type=int, help='JTT PAM number >0')
    group3.add_argument('--tm', type=int, help='Transmembrane PAM number >0')
    group3.add_argument('--aamatrix', type=sysutils.existing_file, help='Path to user-defined AA scoring matrix')
    group3.add_argument('--fmodel', action='store_true', help='Incorporate AA/nuc composition info into scoring matrix')

    group4 = parser.add_argument_group('Options')
    group4.add_argument('--ncpu', type=int, default=1, help='Number of CPU to use')
    group4.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                    (silence stdout and stderr)''')
    group4.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Name for log file (output)')
    group4.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    group4.add_argument('--fastaonly', action='store_true',
                        help='Output fasta files separated by region but do not align')
    group4.add_argument('--alignall', action='store_true', help='Do not separate files by region, align entire file')

    parser.set_defaults(func=multiple_align)


def get_regions(ref_gtf):
    regions = []
    with open(ref_gtf, 'r') as gtffile:
        for line in gtffile:
            reg_name = line.split('\t')[8].split('\"')[1]
            if reg_name not in regions:
                regions.append(reg_name)
    return regions


def generate_fastas(dir_list=None, ref_gtf=None, seqs=None, msadir='.',name='final.fna'):
    ### function to format input data from individual fasta files ###

    ## open dir_list
    if dir_list is not None:
        with open(dir_list, 'r') as f:
            filenames = f.read().splitlines()
    else:
        filenames = []

    ## parse fasta files with SeqIO and append to allseqs list
    allseqs = []
    if seqs is not None:
        for record in SeqIO.parse(seqs, "fasta"):
            allseqs.append(record)
    for n in filenames:
        if len(n) > 0:
            x = os.path.join(n, name)
            for record in SeqIO.parse(x, "fasta"):
                allseqs.append(record)

    ## write all sequences to all_sequences.fasta file
    with open(os.path.join(msadir, 'all_sequences.fasta'), 'w') as outfile1:
        for record in allseqs:
            SeqIO.write(record, outfile1, 'fasta')

    if ref_gtf is not None:
        ## open gtf file and extract region names
        regions = get_regions(ref_gtf)

        ## create dictionary with entry for each region with empty array
        dict = {}
        for x in range(len(regions)):
            dict['region0{0}'.format(x)] = []

        ## append sequences to arrays by region
        for record in SeqIO.parse(os.path.join(msadir, 'all_sequences.fasta'), 'fasta'):
            for i in range(len(regions)):
                if record.id.split('|')[-2] == regions[i]:
                    dict['region' + '0%s' % str(i)].append(record)

        ## write fasta files for each region
        for reg in dict:
            outname = 'all_sequences_' + '%s' % reg + '.fasta'
            with open(os.path.join(msadir, outname), 'w') as outfile2:
                for record in dict[reg]:
                    SeqIO.write(record, outfile2, 'fasta')

        ## return number of regions for use in MAFFT step
        return len(regions)
    else:
        return 1


def run_mafft(inputseqs=None, out_align="alignment.fasta", auto=None, algo=None, sixmerpair=None, globalpair=None,
              localpair=None,
              genafpair=None, fastapair=None, weighti=None, retree=None,
              maxiterate=None, noscore=None, memsave=None, parttree=None, dpparttree=None, fastaparttree=None,
              partsize=None, groupsize=None, lop=None, lep=None, lexp=None, LOP=None,
              LEXP=None, bl=None, jtt=None, tm=None, aamatrix=None, fmodel=None, clustalout=None, inputorder=None,
              reorder=None, treeout=None, quiet_mafft=None, nuc=None, amino=None,
              quiet=False, logfile=None, debug=False, ncpu=1, msadir='.', phylipout=None):
    ### function to run MAFFT ###

    sysutils.check_dependency('mafft')

    ## create MAFFT command using input options
    if algo is None:
        cmd1 = ['mafft', '--thread', '%d' % ncpu, ]
    else:
        if algo not in ['linsi', 'ginsi', 'einsi', 'fftnsi', 'fftns', 'nwns', 'nwnsi']:
            msg = 'Algorithm not in MAFFT'
            raise sysutils.PipelineStepError(msg)
        else:
            cmd1 = ['%s' % algo]
    if clustalout is True:
        cmd1 += ['--clustalout']
    if inputorder is True:
        cmd1 += ['--inputourder']
    if reorder is True:
        cmd1 += ['--reorder']
    if treeout is True:
        cmd1 += ['--treeout']
    if quiet_mafft is True:
        cmd1 += ['--quiet']
    if nuc is True:
        cmd1 += ['--nuc']
    if amino is True:
        cmd1 += ['--amino']

    ### algorithm options
    if auto is True:
        cmd1 += ['--auto']
    if sixmerpair is True:
        cmd1 += ['--6merpair']
    if globalpair is True:
        cmd1 += ['--globalpair']
    if localpair is True:
        cmd1 += ['--localpair']
    if genafpair is True:
        cmd1 += ['--genafpair']
    if fastapair is True:
        cmd1 += ['--fastapair']
    if weighti is not None:
        cmd1 += ['--weighti', '%f' % weighti]
    if retree is not None:
        cmd1 += ['--retree', '%d' % retree]
    if maxiterate is not None:
        cmd1 += ['--maxiterate', '%d' % maxiterate]
    if noscore is True:
        cmd1 += ['--noscore']
    if memsave is True:
        cmd1 += ['--memsave']
    if parttree is True:
        cmd1 += ['--parttree']
    if dpparttree is True:
        cmd1 += ['--dpparttree']
    if fastaparttree is True:
        cmd1 += ['--fastaparttree']
    if partsize is not None:
        cmd1 += ['--partsize', '%d' % partsize]
    if groupsize is not None:
        cmd1 += ['--groupsize', '%d' % groupsize]

    ### parameters
    if lop is not None:
        cmd1 += ['--lop', '%f' % lop]
    if lep is not None:
        cmd1 += ['--lep', '%f' % lep]
    if lexp is not None:
        cmd1 += ['--lexp', '%f' % lexp]
    if LOP is not None:
        cmd1 += ['--LOP', '%f' % LOP]
    if LEXP is not None:
        cmd1 += ['--LEXP', '%f' % LEXP]
    if bl is not None:
        cmd1 += ['--bl', '%d' % bl]
    if jtt is not None:
        cmd1 += ['--jtt', '%d' % jtt]
    if tm is not None:
        cmd1 += ['--tm', '%d' % tm]
    if aamatrix is not None:
        cmd1 += ['--aamatrix', '%s' % aamatrix]
    if fmodel is True:
        cmd1 += ['--fmodel']

    # Outputs
    outName = os.path.join(msadir,
                           '%s' % os.path.basename(out_align))

    ## create command
    cmd1 += ['%s' % os.path.join(msadir, os.path.basename(inputseqs)), '>', '%s' % outName]

    ## run MAFFT command
    sysutils.command_runner([cmd1, ], 'multiple_align', quiet, logfile, debug)

    if phylipout is True:
        phyout = outName[:-6] + '.phy'
        SeqIO.convert(outName, 'fasta', phyout, 'phylip-relaxed')  # relaxed allows for long sequence names
        cmd2 = ['echo', 'Output converted to PHYLIP format from FASTA format.']
        sysutils.command_runner([cmd2, ], 'multiple_align', quiet, logfile, debug)

    if clustalout is True:
        clustout = outName[:-6] + '.aln'
        cmd3 = ['mv', outName, clustout]
        sysutils.command_runner([cmd3, ], 'multiple_align', quiet, logfile, debug)
        cmd4 = ['echo', 'Alignment output is in CLUSTAL format.']
        sysutils.command_runner([cmd4, ], 'multiple_align', quiet, logfile, debug)

    return



def multiple_align(seqs=None, dir_list=None, haplotypes=None,ref_gtf=None, out_align="alignment.fasta", auto=None, algo=None,
                   sixmerpair=None,
                   globalpair=None, localpair=None, genafpair=None, fastapair=None, weighti=None, retree=None,
                   maxiterate=None, noscore=None, memsave=None, parttree=None, dpparttree=None, fastaparttree=None,
                   partsize=None, groupsize=None, lop=None, lep=None, lexp=None, LOP=None,
                   LEXP=None, bl=None, jtt=None, tm=None, aamatrix=None, fmodel=None, clustalout=None, inputorder=None,
                   reorder=None, treeout=None, quiet_mafft=None, nuc=None, amino=None,
                   outdir=".", quiet=False, logfile=None, debug=False, fastaonly=False, alignall=False, ncpu=1,
                   phylipout=None):

    ### RUN MULTIPLE_ALIGN STAGE ###

    if seqs is None and dir_list is None:
        msg = '--seqs or --dir_list is required'
        raise sysutils.MissingRequiredArgument(msg)
    if ref_gtf is None and alignall is False and haplotypes is False:
        msg = 'No GTF file given'
        raise sysutils.MissingRequiredArgument(msg)

    ## create output directory
    msadir = os.path.join(outdir, 'multiple_align')

    cmd2 = ['mkdir -p', msadir]
    sysutils.command_runner([cmd2, ], 'multiple_align', quiet, None, debug)

    ### OPTION 1: aligning haplotype files ###

    # align haplotype fasta files (does not separate)

    if haplotypes is True:
        generate_fastas(dir_list=dir_list, ref_gtf=None, seqs=seqs, msadir=msadir,name='ph_haplotypes.fna')

        #### "Mafft stores the input sequences and other files in a temporary        ####
        #### directory, which by default is located in /tmp." - mafft documentation  ####
        ## Set the temporary directory for mafft
        # Temporary directory
        tempdir = sysutils.create_tempdir('multiple_align', None, quiet, logfile)
        os.environ["MAFFT_TMPDIR"] = tempdir
        run_mafft(inputseqs='all_sequences.fasta', out_align=out_align, auto=auto, algo=algo, sixmerpair=sixmerpair,
                  globalpair=globalpair,
                  localpair=localpair, genafpair=genafpair, fastapair=fastapair, weighti=weighti, retree=retree,
                  maxiterate=maxiterate, noscore=noscore, memsave=memsave, parttree=parttree, dpparttree=dpparttree,
                  fastaparttree=fastaparttree, partsize=partsize, groupsize=groupsize,
                  lop=lop, lep=lep, lexp=lexp, LOP=LOP, LEXP=LEXP, bl=bl, jtt=jtt, tm=tm, aamatrix=aamatrix,
                  fmodel=fmodel, clustalout=clustalout, inputorder=inputorder, reorder=reorder, treeout=treeout,
                  quiet_mafft=quiet_mafft, nuc=nuc, amino=amino, quiet=quiet, logfile=logfile,
                  debug=debug, ncpu=ncpu, msadir=msadir, phylipout=phylipout)
        return

    ### OPTION 2: only generate fasta files, do not align (FASTAONLY = TRUE) ###

    ## if fasta only option is entered, write separated fasta files and end stage

    if fastaonly is True:
        generate_fastas(dir_list, ref_gtf, seqs, msadir)
        return

    ### OPTION 3: do not separate by region before aligning (ALIGNALL = TRUE) ###

    ## if ONLY fasta file is given and alignall is true, run MAFFT and end stage
    if seqs is not None and dir_list is None and alignall is True:

        #### "Mafft stores the input sequences and other files in a temporary        ####
        #### directory, which by default is located in /tmp." - mafft documentation  ####
        ## Set the temporary directory for mafft
        # Temporary directory
        tempdir = sysutils.create_tempdir('multiple_align', None, quiet, logfile)
        os.environ["MAFFT_TMPDIR"] = tempdir

        run_mafft(inputseqs=seqs, out_align=out_align, auto=auto, algo=algo, sixmerpair=sixmerpair,
                  globalpair=globalpair,
                  localpair=localpair, genafpair=genafpair, fastapair=fastapair, weighti=weighti, retree=retree,
                  maxiterate=maxiterate, noscore=noscore, memsave=memsave, parttree=parttree, dpparttree=dpparttree,
                  fastaparttree=fastaparttree, partsize=partsize, groupsize=groupsize,
                  lop=lop, lep=lep, lexp=lexp, LOP=LOP, LEXP=LEXP, bl=bl, jtt=jtt, tm=tm, aamatrix=aamatrix,
                  fmodel=fmodel, clustalout=clustalout, inputorder=inputorder, reorder=reorder, treeout=treeout,
                  quiet_mafft=quiet_mafft, nuc=nuc, amino=amino, quiet=quiet, logfile=logfile,
                  debug=debug, ncpu=ncpu, msadir=msadir, phylipout=phylipout)
        return

    ## if an dir_list is given and alignall is true, generate and then align all_sequences.fasta
    elif alignall is True:
        generate_fastas(dir_list=dir_list, ref_gtf=None, seqs=seqs, msadir=msadir)

        #### "Mafft stores the input sequences and other files in a temporary        ####
        #### directory, which by default is located in /tmp." - mafft documentation  ####
        ## Set the temporary directory for mafft
        # Temporary directory
        tempdir = sysutils.create_tempdir('multiple_align', None, quiet, logfile)
        os.environ["MAFFT_TMPDIR"] = tempdir

        run_mafft(inputseqs='all_sequences.fasta', out_align=out_align, auto=auto, algo=algo, sixmerpair=sixmerpair,
                  globalpair=globalpair,
                  localpair=localpair, genafpair=genafpair, fastapair=fastapair, weighti=weighti, retree=retree,
                  maxiterate=maxiterate, noscore=noscore, memsave=memsave, parttree=parttree, dpparttree=dpparttree,
                  fastaparttree=fastaparttree, partsize=partsize, groupsize=groupsize,
                  lop=lop, lep=lep, lexp=lexp, LOP=LOP, LEXP=LEXP, bl=bl, jtt=jtt, tm=tm, aamatrix=aamatrix,
                  fmodel=fmodel, clustalout=clustalout, inputorder=inputorder, reorder=reorder, treeout=treeout,
                  quiet_mafft=quiet_mafft, nuc=nuc, amino=amino, quiet=quiet, logfile=logfile,
                  debug=debug, ncpu=ncpu, msadir=msadir, phylipout=phylipout)
        return

    ### OPTION 4 (default): separate by region and align each region individually ###

    ## generate separated fasta files from dir_list and seqs + run MAFFT on each region
    n = generate_fastas(dir_list, ref_gtf, seqs, msadir)  # n = number of regions
    for i in range(n):  # iterate over each region
        seqname = 'all_sequences_region' + '0%s' % i + '.fasta'  # consistent with generate_fastas output
        alignmentname = 'alignment_region' + '0%s' % i + '.fasta'

        #### "Mafft stores the input sequences and other files in a temporary        ####
        #### directory, which by default is located in /tmp." - mafft documentation  ####
        ## Set the temporary directory for mafft
        tempdir = sysutils.create_tempdir('multiple_align', None, quiet, logfile)
        os.environ["MAFFT_TMPDIR"] = tempdir

        ## run mafft
        run_mafft(inputseqs=seqname, out_align=alignmentname, auto=auto, algo=algo, sixmerpair=sixmerpair,
                  globalpair=globalpair,
                  localpair=localpair, genafpair=genafpair, fastapair=fastapair, weighti=weighti, retree=retree,
                  maxiterate=maxiterate, noscore=noscore, memsave=memsave, parttree=parttree, dpparttree=dpparttree,
                  fastaparttree=fastaparttree, partsize=partsize, groupsize=groupsize,
                  lop=lop, lep=lep, lexp=lexp, LOP=LOP, LEXP=LEXP, bl=bl, jtt=jtt, tm=tm, aamatrix=aamatrix,
                  fmodel=fmodel, clustalout=clustalout, inputorder=inputorder, reorder=reorder, treeout=treeout,
                  quiet_mafft=quiet_mafft, nuc=nuc, amino=amino, quiet=quiet, logfile=logfile,
                  debug=debug, ncpu=1, msadir=msadir, phylipout=phylipout)

    ## summary message at end
    cmd3 = ['echo', 'Stage completed. Output files are located here: %s\n' % os.path.abspath(msadir)]
    sysutils.command_runner([cmd3, ], 'multiple_align', quiet, logfile, debug)


def console():
    parser = argparse.ArgumentParser(
        description='Align sequences using MAFFT. Required input: --seqs AND/OR --dir_list AND --ref_gtf (unless --alignall option included).',
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
