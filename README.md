# HAPHPIPE

_**HA**plotype and **PH**ylodynamics pipeline for viral assembly, population genetics, and phylodynamics._


## Installation

__1. Create a conda environment with the following dependencies__

```bash
conda create -n haphpipe \
    python \
    future \
    pyyaml \
    biopython \
    seqtk \
    bowtie2 \
    bwa \
    flash \
    freebayes \
    mummer \
    picard \
    trimmomatic \
    samtools=1.9 \
    gatk=3.8 \
    spades \
    blast \
    sierrapy

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

##### `haphpipe_assemble_01`

This pipeline implements amplicon assembly using a denovo approach. Reads are
error-corrected and used to refine the initial assembly, with up to 5
refinement steps.

##### `haphpipe_assemble_02`

This pipeline implements amplicon assembly using a reference-based mapping 
approach. Reads are error-corrected and used to refine the initial assembly,
with up to 5 refinement steps.

See more information regarding the pipelines at the [wiki](https://github.com/gwcbi/haphpipe/wiki/Example-Pipelines).



## Stages

Each stage can be run on its own. Stages are grouped into 4 categories: hp_reads, hp_assemble, hp_haplotype, and hp_annotate.
More detailed description of command line options for each stage are available in the [wiki](https://github.com/gwcbi/haphpipe/wiki).


### hp_reads

Manipulate reads. Input is reads in fastq format, output is modified reads in fastq format.

##### sample_reads

Subsample reads using seqtk. Input is reads in fastq format. Options are the number of reads to sample, fraction of reads to sample (between 0 and 1), and seed. Output is sampled reads in fastq format.

##### trim_reads

Trim reads using Trimmomatic. Input is reads in fastq format. Options are adapter file, trimmers, and encoding. Output is trimmed reads in fastq format.

##### join_reads

Join reads using FLASH. Input is reads in fastq format. Output is joined reads in fastq format.

##### ec_reads

Error correction using spades. Input is reads in fastq format. Output is error-corrected reads in fastq format.


### hp_assemble

Assemble consensus sequence(s). Input reads (in fastq format) are assembled 
using either denovo assembly or reference-based alignment. 
Resulting consensus can be further refined.

##### assemble_denovo

Assemble reads using denovo assembly. Input is reads in FASTQ format. Output is contigs in FNA format.

##### assemble_amplicons

Assemble contigs using reference sequence and amplicon regions. Input is contigs and reference sequence in FASTA format and amplicon regions in GTF format.

##### assemble_scaffold

Scaffold contigs using reference sequence. Input is contigs in FASTA format and reference sequence in FASTA format. Output is scaffold assembly, alligned scaffold, imputed scaffold, and padded scaffold in FASTA format.

##### align_reads

Map reads to reference sequence. Input is reads in FASTQ format and reference sequence in FASTA format. 

##### call_variants

Variant calling from alignment. Input is alignment file in BAM format and reference sequence in FASTA format. Output is a VCF format file. 

##### vcf_to_consensus

Generate a consensus sequence from a VCF file. Input is a VCF file. Output is the consensus sequence in FASTA format. 

##### refine_assembly

Map reads to a denovo assembly or reference alignment. Assembly or alignment is iteratively updated. Input is reads in FASTQ format and reference sequence (assembly or reference alignment) in FASTA format. Output is refined assembly in 

##### finalize_assembly

Finalize consensus, map reads to consensus, and call variants. Input is reads in FASTQ format and reference sequence in FASTA format. Output is finalized reference sequence, alignment, and variants (in FASTA, BAM, and VCF formats, respectively).


### hp_haplotype

Reconstructing haplotypes. This is an active area of research. The user is required to install PredictHaplo on their own, prior to running these stages.

##### predict_haplo

Assemble haplotypes with PredictHaplo. Input are the paired reads in FASTQ format and a reference sequence in FASTA format that is used to align the reads to. Output is reconstructed haplotypes in PredicHaplo's formats. The best (aka longest reconstructed haplotype) is used as the input for the complementary ph_parser stage below. 

##### ph_parser

Parse output from PredictHaplo. Input is the best.fsa file from the previous step, and the output includes a FASTA format with the reconstructed haplotypes and a summary text file that contains the number of haplotypes reconstructed, a Waterson's genetic diversity estimate, and the length of the reconstructed haplotypes. 


### hp_annotate

Annotate consensus sequences.

##### pairwise_align 

##### extract_pairwise

##### annotate_from_ref


