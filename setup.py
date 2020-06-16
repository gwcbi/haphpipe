# -*- coding: utf-8 -*-

# Copyright (C) 2019 Matthew L. Bendall and (C) 2020 Keylie M. Gibson

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

import setuptools

from haphpipe._version import VERSION

setuptools.setup(name='haphpipe',
      version=VERSION,
      description='HAplotype and PHylodynamics pipeline for viral assembly, population genetics, and phylodynamics.',
      url='https://github.com/gwcbi/haphpipe',
      author='Matthew L. Bendall, Keylie M. Gibson, Maggie C. Steiner',
      author_email='bendall@gwu.edu and kmgibson@gwu.edu',
      packages=setuptools.find_packages(),
      package_data={'': ['bin/refs/HIV_B.K03455.HXB2*']},
      include_package_data=True,
      scripts=[
          'bin/haphpipe_assemble_01',
          'bin/haphpipe_assemble_02',
          'bin/haphpipe_demo',
      ],
      entry_points={
          'console_scripts': [
              'haphpipe=haphpipe.haphpipe:console',
              # miscellaneous subcommands
              'hp_demo=haphpipe.stages.demo:console',
              # hp_reads subcommands
              'hp_sample_reads=haphpipe.stages.sample_reads:console',
              'hp_trim_reads=haphpipe.stages.trim_reads:console',
              'hp_join_reads=haphpipe.stages.join_reads:console',
              'hp_ec_reads=haphpipe.stages.ec_reads:console',
              # hp_assemble subcommands
              'hp_assemble_denovo=haphpipe.stages.assemble_denovo:console',
              'hp_assemble_amplicons=haphpipe.stages.assemble_amplicons:console',
              'hp_assemble_scaffold=haphpipe.stages.assemble_scaffold:console',
              'hp_align_reads=haphpipe.stages.align_reads:console',
              'hp_call_variants=haphpipe.stages.call_variants:console',
              'hp_vcf_to_consensus=haphpipe.stages.vcf_to_consensus:console',
              'hp_refine_assembly=haphpipe.stages.refine_assembly:console',
              'hp_finalize_assembly=haphpipe.stages.finalize_assembly:console',
              # hp_annotate subcommands
              'hp_pairwise_align=haphpipe.stages.pairwise_align:console',
              'hp_extract_pairwise=haphpipe.stages.extract_pairwise:console',
              'hp_summary_stats=haphpipe.stages.summary_stats:console',
              # hp_haplotype subcommands
              'hp_predict_haplo=haphpipe.stages.predict_haplo:console',
              'hp_ph_parser=haphpipe.stages.ph_parser:console',
              'hp_cliquesnv=haphpipe.stages.cliquesnv:console',
              # hp_phylo subcommands
              'hp_multiple_align=haphpipe.stages.multiple_align:console',
              'hp_model_test=haphpipe.stages.model_test:console',
              'hp_build_tree_NG=haphpipe.stages.build_tree_NG:console',
          ],
      },
      zip_safe=False,
)
