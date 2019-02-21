# HAPHPIPE

_**HA**plotype and **PH**ylodynamics pipeline for viral assembly, population genetics, and phylodynamics._


## Installation

__1. Create a conda environment with the following dependencies__

```bash
conda create -n haphpipe2 \
    python=2.7 \
    future \
    pyyaml \
    biopython \
    seqtk \
    bowtie2 \
    flash \
    freebayes \
    mummer \
    picard \
    trimmomatic \
    samtools=1.9 \
    gatk=3.8 \
    spades \
    blast

```

```
conda install biopython=1.73 bowtie2=2.3.4.3 flash=1.2.11 freebayes=1.2.0 gatk=3.8 mummer=3.23 picard=2.18.26 samtools=1.9 trimmomatic=0.38 trinity=date.2011_11_26
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

### hp_reads

##### hp_sample_reads

##### hp_trim_reads

##### hp_join_reads

##### hp_ec_reads

### hp_assemble

####



