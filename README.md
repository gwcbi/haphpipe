# HAPHPIPE

*__HA__plotype and __PH__ylodynamics pipeline for viral assembly, population genetics, and phylodynamics.*


## Installation

__1. Create a conda environment with the following dependencies__

```bash
conda create -n haphpipe \
    python=2.7 \
    future \
    biopython \
    bowtie2 \
    flash \
    freebayes \
    mummer \
    picard \
    samtools \
    trimmomatic \
    trinity \
    gatk
```

__2. Activate the environment__

```
conda activate haphpipe
```

__3. Install GATK__

Due to license restrictions, bioconda cannot distribute
and install GATK directly. To fully install GATK, you must
download a licensed copy of GATK (version 3.8-0) from the Broad Institute:
[https://software.broadinstitute.org/gatk/download/archive](https://software.broadinstitute.org/gatk/download/archive).

Register the package using gatk3-register:

```
gatk3-register /path/to/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2
```

This will copy GATK into your conda environment.

NOTE: HAPHPIPE was developed and tested using GATK 3.8.

__4. Install HAPHPIPE__

```
pip install git+git://github.com/gwcbi/haphpipe.git
```

## Pipelines

Describe haphpipe amplicon, assembly, and other pipelines

## Stages

Describe each haphpipe stage

