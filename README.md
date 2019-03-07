# HAPHPIPE

_**HA**plotype and **PH**ylodynamics pipeline for viral assembly, population genetics, and phylodynamics._


## Installation

__1. Create a conda environment with the following dependencies__

```bash
conda create -n haphpipe \
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


[TODO] describe pipelines

##### `haphpipe_assemble_01`




## Stages

Each stage can be run on its own. Stages are grouped into categories...

[TODO] Describe each haphpipe stage

### hp_reads

Manipulate reads. Input is reads in fastq format, output is modified reads in fastq format.

##### hp_sample_reads

Subsample reads using seqtk. Input is reads in fastq format. Options are the number of reads to sample, fraction of reads to sample (between 0 and 1), and seed. Output is sampled reads in fastq format.

##### hp_trim_reads

Trim reads using Trimmomatic. Input is reads in fastq format. Options are adapter file, trimmers, and encoding. Output is trimmed reads in fastq format.

##### hp_join_reads

Join reads using FLASH. Input is reads in fastq format. Output is joined reads in fastq format.

##### hp_ec_reads

Error correction using spades. Input is reads in fastq format. Output is error-corrected reads in fastq format.

### hp_assemble

Assemble consensus sequence(s). Input reads (in fastq format) are assembled 
using either denovo assembly or reference-based alignment. 
Resulting consensus can be further refined.

##### hp_assemble_denovo

##### hp_assemble_amplicons

##### hp_assemble_scaffold

##### hp_align_reads

##### hp_call_variants

##### hp_vcf_to_consensus

##### hp_refine_assembly

##### hp_finalize_assembly


### hp_annotate

Annotation

### hp_haplotype

Haplotype
