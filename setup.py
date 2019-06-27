import setuptools

setuptools.setup(name='haphpipe',
      version='0.8',
      description='HAplotype and PHylodynamics pipeline for viral assembly, population genetics, and phylodynamics.',
      url='https://github.com/gwcbi/haphpipe',
      author='Matthew L. Bendall',
      author_email='bendall@gwu.edu',
      packages=setuptools.find_packages(),
      scripts=[
          'bin/haphpipe_assemble_01',
          'bin/haphpipe_assemble_02',
      ],
      entry_points={
          'console_scripts': [
              # Stage groups
              'hp_reads=haphpipe.reads:console',
              'hp_assemble=haphpipe.assemble:main',
              'hp_annotate=haphpipe.annotate:main',
              'hp_haplotype=haphpipe.haplotype:main',
              # hp_reads subcommands
              'hp_sample_reads=haphpipe.stages.sample_reads:console',
              'hp_trim_reads=haphpipe.stages.trim_reads:console',
              'hp_join_reads=haphpipe.stages.join_reads:console',
              'hp_ec_reads=haphpipe.stages.ec_reads:console',
              #hp_assemble subcommands
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
              # hp_haplotype subcommands
              'hp_predict_haplo=haphpipe.stages.predict_haplo:console',
          ],
      },
      zip_safe=False,
)
