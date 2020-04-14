# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import argparse

from haphpipe.utils import sysutils
from haphpipe.utils import sequtils
from haphpipe.utils import alignutils


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"

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
    group1.add_argument('--contigs_fa', type=sysutils.existing_file,
                        required=True,
                        help='Fasta file with assembled contigs')
    group1.add_argument('--ref_fa', type=sysutils.existing_file,
                        required=True,
                        help='''Fasta file with reference genome to scaffold
                                against''')
    group1.add_argument('--outdir', type=sysutils.existing_dir, default='.',
                        help='Output directory')
    
    group2 = parser.add_argument_group('Scaffold options')
    group2.add_argument('--seqname', default='sample01',
                        help='Name to append to scaffold sequence.')
    
    group3 = parser.add_argument_group('Settings')
    group3.add_argument('--keep_tmp', action='store_true',
                        help='Additional options')
    group3.add_argument('--quiet', action='store_true',
                        help='''Do not write output to console
                                (silence stdout and stderr)''')
    group3.add_argument('--logfile', type=argparse.FileType('a'),
                        help='Append console output to this file')
    group3.add_argument('--debug', action='store_true',
                        help='Print commands but do not run')

    parser.set_defaults(func=assemble_scaffold)


def assemble_scaffold(
        contigs_fa=None, ref_fa=None, outdir='.',
        seqname='sample01',
        keep_tmp=False, quiet=False, logfile=None, debug=False
    ):
    """ Pipeline step to assemble contigs to reference scaffold

    Args:
        contigs_fa (str): Path to fasta file with assembled contigs
        ref_fa (str): Path to reference fasta file
        outdir (str): Path to output directory
        seqname (str): Name to append to scaffold sequence
        keep_tmp (bool): Do not delete temporary directory
        quiet (bool): Do not write output to console
        logfile (file): Append console output to this file
        debug (bool): Print commands but do not run

    Returns:
        out_scaffold (str): Path to scaffold FASTA. Reference positions that
                            were not covered have 'n'
        out_imputed (str):  Path to imputed FASTA. Reference positions that
                            were not covered have reference base.
        out_aln (str):      Path to FASTA alignment between scaffold and
                            reference.
        out_padded (str):   Path to output with all contigs aligned to
                            reference.
    """
    # Check dependencies
    sysutils.check_dependency('nucmer')
    sysutils.check_dependency('delta-filter')
    sysutils.check_dependency('show-tiling')
    
    # Outputs
    out_scaffold = os.path.join(outdir, 'scaffold_assembly.fa')
    out_imputed = os.path.join(outdir, 'scaffold_imputed.fa')
    out_aln = os.path.join(outdir, 'scaffold_aligned.fa')
    out_padded = os.path.join(outdir, 'scaffold_padded.out')
    
    # Temporary directory
    tempdir = sysutils.create_tempdir(
        'assemble_scaffold', None, quiet, logfile
    )

    # Create fasta file with sequence IDs only (remove decription)
    tmp_contigs_fa = sequtils.clean_seqnames_file(contigs_fa, tempdir)

    with open(out_padded, 'w') as pad_fh:
        scaffolds = alignutils.assemble_to_ref(
            tmp_contigs_fa, ref_fa, tempdir, pad_fh=pad_fh,
            quiet=quiet, logfile=logfile, debug=debug
        )

    # Output scaffolds as FASTA
    with open(out_scaffold, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            s = scaffolds[ref].scaffold()
            print('>%s\n%s' % (n, sequtils.wrap(s)), file=outh)

    # Output imputed as FASTA
    with open(out_imputed, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            s = scaffolds[ref].imputed()
            print('>%s\n%s' % (n, sequtils.wrap(s)), file=outh)

    # Output alignments for other pipeline stages
    with open(out_aln, 'w') as outh:
        for ref in sorted(scaffolds.keys()):
            n = '%s.%s' % (ref.split('.')[0], seqname)
            print('>REF|%s\n%s' % (n, scaffolds[ref].raln()), file=outh)
            print('>%s\n%s' % (n, scaffolds[ref].qaln()), file=outh)

    if not keep_tmp:
        sysutils.remove_tempdir(tempdir, 'assemble_scaffold', quiet, logfile)

    return out_scaffold, out_imputed, out_aln, out_padded


def console():
    """ Entry point

    Returns:
        None

    """
    parser = argparse.ArgumentParser(
        description='Assemble contigs to genome.',
        formatter_class=sysutils.ArgumentDefaultsHelpFormatterSkipNone,
    )
    stageparser(parser)
    args = parser.parse_args()
    args.func(**sysutils.args_params(args))


if __name__ == '__main__':
    console()
