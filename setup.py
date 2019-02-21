import setuptools

setuptools.setup(name='haphpipe',
      version='0.8',
      description='HAplotype and PHylodynamics pipeline for viral assembly, population genetics, and phylodynamics.',
      url='https://github.com/gwcbi/haphpipe',
      author='Matthew L. Bendall',
      author_email='bendall@gwu.edu',
      packages=setuptools.find_packages(),
      scripts=[
          'bin/haphpipe_asm.sh',
          'bin/haphpipe_predicthaplo.sh',
          'bin/haphpipe_amp.sh',
      ],
      entry_points={
          'console_scripts': [
              # Stage groups
              'hp_reads=haphpipe.reads:console',
              'hp_assemble=haphpipe.assembly:main',
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
              'hp_finalize_consensus=haphpipe.stages.fix_consensus:console'
          ],
      },
      zip_safe=False,
)
