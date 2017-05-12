from setuptools import setup

setup(name='haphpipe',
      version='0.1',
      description='HAplotype and PHylodynamics pipeline for viral assembly, population genetics, and phylodynamics.',
      url='https://github.com/gwcbi/haphpipe',
      author='Matthew L. Bendall',
      author_email='bendall@gwu.edu',
      packages=['haphpipe'],
      scripts=[
          'bin/haphpipe_asm.sh',
          'bin/haphpipe_predicthaplo.sh',
          'bin/haphpipe_amp.sh',
      ],
      entry_points={
          'console_scripts': [
              'hp_assemble=haphpipe.assembly:main',
              'hp_haplotype=haphpipe.haplotype:main',
          ],
      },
      zip_safe=False,
)
