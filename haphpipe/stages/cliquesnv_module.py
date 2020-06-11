from __future__ import print_function
import os
import sys
import argparse
import shutil

from Bio import SeqIO

from haphpipe.utils import sysutils
from haphpipe.utils.sysutils import MissingRequiredArgument

__author__ = 'Margaret C. Steiner and Matthew L. Bendall'
__copyright__ = 'Copyright (C) 2020 Margaret C. Steiner and (C) 2019 Matthew L. Bendall'

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
    group1.add_argument('--fq1',type=sysutils.existing_file,
                        help='Fastq file with read 1 OR single read fasta file')
    group1.add_argument('--fq2', type=sysutils.existing_file,
                        help='Fastq file with read 2')
    group1.add_argument('--ref_fa',type=sysutils.existing_file,help="Reference FASTA file")
    group1.add_argument('--outdir',type=sysutils.existing_dir,default='.',help='Output directory')
    group1.add_argument('--single', action='store_true', help="Single read date (Default: FALSE, e.g. paired-end)")

    group2 = parser.add_argument_group('CliqueSNV Options')
    group2.add_argument('--jardir',type=str,default='.',help='Path to clique-snv.jar (existing) (Default: current directory)')
    group2.add_argument('--O22min',type=float,help="minimum threshold for O22 value")
    group2.add_argument('--O22minfreq',type=float,help="minimum threshold for O22 frequency relative to read coverage")
    group2.add_argument('--printlog',action='store_true',help="Print log data to console")
    group2.add_argument('--merging',type=str,help='Cliques merging algorithm: accurate or fast')
    group2.add_argument('--fasta_format',type=str,default='extended4',help="Fasta defline format: short or extended, add number at end to adjust precision of frequency")
    group2.add_argument('--outputstart',type=int,help="Output start position")
    group2.add_argument('--outputend', type=int, help="Output end position")

    group3 = parser.add_argument_group('HAPHPIPE Options')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                        (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Name for log file (output)')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')
    group3.add_argument('--ncpu', type=int, default=1, help='Number of CPU to use')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Do not delete temporary directory')

    parser.set_defaults(func=cliquesnv)

def cliquesnv(fq1=None,fq2=None,ref_fa=None,outdir='.',jardir='.',O22min=None,O22minfreq=None,printlog=None,single=False,
              merging=None,fasta_format='extended4',outputstart=None,outputend=None,keep_tmp=False, quiet=False, logfile=None, debug=False,ncpu=1):

    # check dependencies and required arguments
    if single == False and (fq1 is None or fq2 is None):
        raise MissingRequiredArgument("Either fq1 or fq2 missing.")
    if single == True and fq1 is None:
        raise MissingRequiredArgument("Read file fq1 missing.")
    if ref_fa is None:
        raise MissingRequiredArgument("Reference FASTA missing.")

    sysutils.check_dependency('samtools')
    sysutils.check_dependency('bwa')

    if(os.path.isfile(os.path.join(jardir,"clique-snv.jar"))):
        print("CliqueSNV JAR file found.")
    else:
        raise MissingRequiredArgument("No JAR file found.")

    # Temporary directory
    tempdir = sysutils.create_tempdir('clique_snv', None, quiet, logfile)

    # Load reference fasta
    refs = {s.id: s for s in SeqIO.parse(ref_fa, 'fasta')}

    # Identify reconstruction regions
    regions = []
    for rname, s in refs.items():
        regions.append(('cs%02d' % (len(regions) + 1), rname, 1, len(s)))

    sysutils.log_message('[--- Haplotype Reconstruction Regions ---]\n', quiet, logfile)
    for iv in regions:
        sysutils.log_message('%s -- %s:%d-%d\n' % iv, quiet, logfile)

    if single == False: #paired end
        # remove .1 and .2 from read names
        fq1_c = os.path.join(tempdir,"fq1_corrected.fastq")
        fq2_c = os.path.join(tempdir, "fq2_corrected.fastq")
        cmd01 = ["cat %s | sed 's/\.1 / /' > %s" % (fq1,fq1_c)]
        cmd02 = ["cat %s | sed 's/\.2 / /' > %s" % (fq2,fq2_c)]
        sysutils.command_runner([cmd01,cmd02],'clique_snv:setup',quiet,logfile,debug)

        # Create alignment for each REFERENCE in the reconstruction regions
        alnmap = {}
        for cs, rname, spos, epos in regions:
            if rname not in alnmap:
                # Create alignment
                tmp_ref_fa = os.path.join(tempdir, 'ref.%d.fa' % len(alnmap))
                tmp_sam = os.path.join(tempdir, 'aligned.%d.sam' % len(alnmap))
                SeqIO.write(refs[rname], tmp_ref_fa, 'fasta')
                cmd1 = ['bwa', 'index', tmp_ref_fa, ]
                cmd2 = ['bwa', 'mem', tmp_ref_fa, fq1_c, fq2_c, '|', 'samtools', 'view', '-h', '-F', '12', '>', tmp_sam, ]
                cmd3 = ['rm', '-f', '%s.*' % tmp_ref_fa]
                sysutils.command_runner(
                    [cmd1, cmd2, cmd3], 'clique_snv:setup', quiet, logfile, debug
                )
                alnmap[rname] = (tmp_ref_fa, tmp_sam)

    else: #single read

        # Create alignment for each REFERENCE in the reconstruction regions
        alnmap = {}
        for cs, rname, spos, epos in regions:
            if rname not in alnmap:
                # Create alignment
                tmp_ref_fa = os.path.join(tempdir, 'ref.%d.fa' % len(alnmap))
                tmp_sam = os.path.join(tempdir, 'aligned.%d.sam' % len(alnmap))
                SeqIO.write(refs[rname], tmp_ref_fa, 'fasta')
                cmd1 = ['bwa', 'index', tmp_ref_fa, ]
                cmd2 = ['bwa', 'mem', tmp_ref_fa, fq1, '|', 'samtools', 'view', '-h', '-F', '12', '>',
                        tmp_sam, ]
                cmd3 = ['rm', '-f', '%s.*' % tmp_ref_fa]
                sysutils.command_runner(
                    [cmd1, cmd2, cmd3], 'clique_snv:setup', quiet, logfile, debug
                )
                alnmap[rname] = (tmp_ref_fa, tmp_sam)


    # Run CliqueSNV for each region
    cmd4 = ['mkdir -p %s' % os.path.join(outdir, 'clique_snv')]
    sysutils.command_runner([cmd4, ], stage='cliquesnv', quiet=quiet, logfile=logfile, debug=debug)
    i=0 #index for filenames
    for cs, rname, spos, epos in regions:
        msg = "Reconstruction region %s:" % cs
        msg += " %s:%d-%d\n" % (rname, spos, epos)
        sysutils.log_message(msg, quiet, logfile)

        # rename the cliquesnv number (cs##) to include region (now: cs##_reg)
        cs = '%s_%s' % (cs, rname.split('|')[-2])

        samfile = os.path.join(tempdir,'aligned.%d.sam' % i)
        method = 'snv-illumina'
        cmd5 = ['java -jar %s -m %s -in %s -threads %d -outDir %s -fdf %s'
                % (os.path.join(jardir,'clique-snv.jar'),method,samfile,
                   ncpu,tempdir,fasta_format)]
        if O22min is not None:
            cmd5 += ['-t %f' % O22min]
        if O22minfreq is not None:
            cmd5 += ['-tf %f' % O22minfreq]
        if printlog is not None:
            cmd5 += ['-log']
        if merging is not None:
            cmd5 += ['-cm %s' % merging]
        if outputstart is not None:
            cmd5 += ['-os %d' % outputstart]
        if outputend is not None:
            cmd5 += ['-oe %d' % outputend]
        sysutils.command_runner([cmd5, ], stage='clique_snv', quiet=quiet, logfile=logfile, debug=debug)

        # copy output file and delete tempdir
        outname1 = 'aligned.%d.txt' % i
        outname2 = 'aligned.%d.fasta' % i

        os.makedirs(os.path.join(outdir, 'clique_snv/%s' % cs), exist_ok=True)
        if os.path.exists(os.path.join(tempdir, '%s' % outname1)):
            shutil.copy(os.path.join(tempdir, '%s' % outname1), os.path.join(outdir,'clique_snv/%s/%s.txt' % (cs,cs)))
        if os.path.exists(os.path.join(tempdir, '%s' % outname2)):
            shutil.copy(os.path.join(tempdir, '%s' % outname2), os.path.join(outdir,'clique_snv/%s/%s.fasta' % (cs,cs)))


        # parse output file
        with open(os.path.join(outdir,'clique_snv/%s/%s_summary.txt' % (cs,cs)),'w') as sumfile, open(os.path.join(outdir,'clique_snv/%s/%s.txt' % (cs,cs)),'r') as infile:
            l = infile.readlines()
            print(l)
            freqs = []
            haps = []
            tempnum=''
            for line in l:
                if "SNV got" in line:
                    print(line)
                    tempnum = line.split(' ')[2]
                if "frequency" in line:
                    freqs += [line.split(' ')[2][:-2]]
                if "haplotype=" in line:
                    haps += [line.split('=')[1][1:-2]]
            sumfile.write('snv_num_haps\t%s\n' % tempnum)
            for k in range(len(freqs)):
                sumfile.write('freq_hap_%d\t%s\n' % (k,freqs[k]))
                sumfile.write('len_hap_%d\t%s\n' % (k,len(haps[k])))

        with open(os.path.join(outdir, 'clique_snv/%s/%s.fasta' % (cs, cs)), 'r') as fastafile:
            fastadata=fastafile.read().replace('aligned.%d' % i,rname)
            with open(os.path.join(outdir, 'clique_snv/%s/%s.fasta' % (cs, cs)), 'w') as newfastafile:
                newfastafile.write(fastadata)

        i += 1

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'clique_snv', quiet, logfile)

    return

def console():
    parser = argparse.ArgumentParser(
        description='Haplotype reconstruction with CliqueSNV.',
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
