from __future__ import print_function
import os
import argparse
import sys
from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import MissingRequiredArgument
from Bio import SeqIO

__author__ = 'Margaret C. Steiner, Keylie M. Gibson, and Matthew L. Bendall'
__copyright__ = 'Copyright (C) 2020 Margaret C. Steiner; (C) 2019 Keylie M. Gibson and Matthew L. Bendall'


def stageparser(parser):
    ### input options with argparse ###

    group1 = parser.add_argument_group('Input/Output')
    group1.add_argument('--dir_list', type=sysutils.existing_file,
                        help='List of directories which include the required files, one on each line')
    group1.add_argument('--ph_list', type=sysutils.existing_file,
                        help='List of directories which include haplotype summary files, one on each line')
    group1.add_argument('--amplicons', action='store_true', help='Amplicons used in assembly')
    group1.add_argument('--outdir', type=sysutils.existing_dir, help='Output directory')
    parser.set_defaults(func=summary_stats)

    group2 = parser.add_argument_group('Settings')
    group2.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console (silence stdout and stderr)''')
    group2.add_argument('--logfile', type=str, help='Name for log file (output)')
    group2.add_argument('--debug', action='store_true', help='Print commands but do not run')


### utility functions: search_file and parse_vcf_file

def search_file(infile, phrase):
    with open(infile, 'r') as infile:
        lines = infile.read().splitlines()
        for line in lines:
            if phrase in line:
                return line


def parse_vcf_file(file, reg):
    with open(file, 'r') as infile:
        count = 0
        lines = infile.read().splitlines()
        for line in lines:
            if line[0] is not '#' and len(line) > 0 and reg in line:
                count += 1
    return count


def summary_stats(dir_list=None, ph_list=None, quiet=False, logfile=None, debug=False, amplicons=False, outdir='.'):
    # check for samtools
    sysutils.check_dependency('samtools')

    # check for dir_list (required)
    if dir_list is not None:
        f = open(dir_list, 'r')
        filenames = f.read().splitlines()
    else:
        msg = 'no directory list given'
        raise MissingRequiredArgument(msg)

    # count number of samples
    numsamps = 0
    for f in filenames:
        if len(f) > 0:
            numsamps += 1

    # count number of PH files
    numph = 0
    if ph_list is not None:
        p = open(ph_list, 'r')
        phnames = p.read().splitlines()
        for f in phnames:
            if len(f) > 0:
                numph += 1

    tsv_header = []
    tsv_samps = []

    with open(os.path.join(outdir, 'summary_stats.txt'), 'w') as outfile:
        for i in range(numsamps):  # for each sample

            # set file names
            bowtiefile = os.path.join(filenames[i], 'final_bt2.out')
            trimfile = os.path.join(filenames[i], 'trimmomatic_summary.out')
            bamfile = os.path.join(filenames[i], 'final.bam')
            outidxstat = os.path.join(filenames[i], 'final.idxstat.txt')
            finalfina = os.path.join(filenames[i], 'final.fna')
            vcfzipped = os.path.join(filenames[i], 'final.vcf.gz')
            vcfunzipped = os.path.join(filenames[i], 'final.vcf')

            sampname = str(filenames[i])
            num_cols = sampname.count('/') + 1

            if i == 0:  # if the first iteration, create tsv_header
                for x in range(num_cols):
                    tsv_header += ['dir_%s' % str(x)]
                tsv_header += ['RAW', 'CLEAN', 'ALN_RATE']

            # output block 1
            outfile.write("SAMPLE " + "%s:\n" % sampname)
            outfile.write("\t Directory: %s\n" % str(os.path.abspath(filenames[i])))
            raw = search_file(trimfile, "Input Read Pairs").split(' ')[3]
            outfile.write("\t Number of raw read pairs: %s\n" % raw)
            cleaned = search_file(bowtiefile, "reads;").split(' ')[0]
            outfile.write("\t Number of cleaned read pairs: %s\n" % cleaned)
            aln_rate = search_file(bowtiefile, "overall alignment rate").split(' ')[0]
            outfile.write("\t Overall alignment rate: %s\n" % aln_rate)

            # create tsv line
            tsv_samp_temp = []
            tsv_samp_temp += sampname.split('/')
            tsv_samp_temp += [str(raw), str(cleaned), str(aln_rate)]

            # index bam file with samtools
            cmd0 = ["samtools index %s" % bamfile]
            sysutils.command_runner([cmd0, ], 'summary_stats', quiet, logfile, debug)

            # run idxstats with samtools
            cmd1 = ["samtools idxstats %s > %s" % (bamfile, outidxstat)]
            sysutils.command_runner([cmd1, ], 'summary_stats', quiet, logfile, debug)

            # unzip vcf file
            if os.path.isfile(vcfzipped):
                cmd2 = ["gunzip %s" % vcfzipped]
                sysutils.command_runner([cmd2, ], 'summary_stats', quiet, logfile, debug)

            # if amplicon assembly
            if amplicons is True:
                all_amplicons = []
                for record in SeqIO.parse(finalfina, 'fasta'):
                    reg_short = record.name.split('|')[5]
                    all_amplicons.append(str(reg_short))

                    # parse outidxstat and output
                    outfile.write("\t\t Amplicon %s:\n" % reg_short)
                    leng = search_file(outidxstat, reg_short).split('\t')[1]
                    outfile.write("\t\t\t Amplicon length: %s\n" % leng)
                    count = search_file(outidxstat, reg_short).split('\t')[2]
                    outfile.write("\t\t\t Amplicon read count: %s\n" % count)

                    # run depth with samtools for coverage
                    dep = os.path.join(filenames[i], 'final.depth.%s.txt' % reg_short)
                    cmd3 = ["samtools depth -r '%s' %s > %s" % (str(record.name), bamfile, dep)]
                    sysutils.command_runner([cmd3, ], 'summary_stats', quiet, logfile, debug)

                    # parse dep file from samtools
                    lines = 0
                    with open(dep) as depfile:
                        for line in depfile:
                            if len(line) > 0:
                                lines += 1
                    perc = (lines / int(leng)) * 100

                    # output coverage
                    outfile.write("\t\t\t Amplicon coverage: %s (%s percent)\n" % (lines, perc))
                    snps = parse_vcf_file(vcfunzipped, reg_short)
                    outfile.write("\t\t\t Number of SNPS: %s\n" % snps)
                    theta = float(snps) / float(leng)
                    outfile.write("\t\t\t Theta: %1.5f\n\n" % theta)

                    # add to tsv line
                    tsv_samp_temp += [str(leng), str(count), str(lines), str(perc), str(snps), str(theta)]

            # add line to list for tsv
            tsv_samps += [tsv_samp_temp]

        # HAPLOTYPE FILES
        if numph > 0:
            outfile.write("\n\nHAPLPOTYPE SUMMARY STATISTICS\n\n")

            ph_tsv_header = []
            ph_tsv_samps = []

            for i in range(numph):  # for each PH directory
                phfile = os.path.join(phnames[i], 'ph_summary.txt')
                num_cols = phnames[i].count('/') + 1

                if i == 0:  # if the first iteration, create ph_tsv_header
                    for x in range(num_cols):
                        ph_tsv_header += ['dir_%s' % str(x)]
                    ph_tsv_header += ['PH_NUM_HAP', 'PH_HAP_DIVERSITY', 'PH_SEQ_LEN']

                # output from parsing phfile
                outfile.write("PH OUTPUT FILE %s:\n" % phfile)
                num_hap = search_file(phfile, "PH_num_hap").split(' ')[1]
                outfile.write("\t Number of haplotypes: %s\n" % num_hap)
                div = search_file(phfile, "PH_hap_diversity").split(' ')[1]
                outfile.write("\t Haplotype diversity: %s\n" % div)
                seq_len = search_file(phfile, "PH_seq_len").split(' ')[1]
                outfile.write("\t Sequence length: %s\n" % seq_len)

                # create tsv line
                ph_tsv_samp_temp = []
                ph_tsv_samp_temp += phnames[i].split('/')
                ph_tsv_samp_temp += [str(num_hap), str(div), str(seq_len)]

                # add line to ph_tsv_samps
                ph_tsv_samps += [ph_tsv_samp_temp]

    # make summary_stats.tsv file
    with open(os.path.join(outdir, 'summary_stats.tsv'), 'w') as outfile:
        if amplicons is True:
            for amp in all_amplicons:
                tsv_header += ['%s_LEN' % amp, '%s_RC' % amp, '%s_COV_NUM' % amp,
                               '%s_COV_PERC' % amp, '%s_SNPS' % amp, '%s_THETA' % amp]
        outfile.write(('\t').join(tsv_header) + '\n')
        for samp in tsv_samps:
            outfile.write(('\t').join(samp) + '\n')

    # make PH_summary_stats.tsv file
    if ph_list is not None:
        with open(os.path.join(outdir, 'PH_summary_stats.tsv'), 'w') as outfile:
            outfile.write(('\t').join(ph_tsv_header) + '\n')
            for samp in ph_tsv_samps:
                outfile.write(('\t').join(samp) + '\n')

    # ending summary message
    cmd3 = ['echo', 'Stage completed. Summary stats are located here: %s\n' % os.path.abspath('summary_stats.txt')]
    if amplicons is True:
        cmd3 += ['echo', 'Amplicons: %s\n' % (', ').join(all_amplicons)]
    sysutils.command_runner([cmd3, ], 'summary_stats', quiet, logfile, debug)


def console():
    parser = argparse.ArgumentParser(
        description='Generate summary statistics. Required input: dir_list and/or ph_list and ref_gtf if amplicons used.',
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
