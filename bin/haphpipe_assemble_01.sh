#! /bin/bash

###############################################################################
# Example pipeline
#
# Usage:
###############################################################################

# Variables
raw1=PL-12_S3_L001_R1_001.fastq.gz
raw2=PL-12_S3_L001_R2_001.fastq.gz
refFA=ref/HIV_B.K03455.HXB2.fasta
refGTF=ref/HIV_B.K03455.HXB2.gtf
outdir=testpipe001
quiet='--quiet'

mkdir -p ${outdir}

hp_reads sample_reads \
    --fq1 ${raw1} \
    --fq2 ${raw2} \
    --nreads 10000 \
    --seed 999 \
    ${quiet} --logfile ${outdir}/haphpipe.out \
    --outdir ${outdir}

hp_reads trim_reads \
    --fq1 ${outdir}/sample_1.fastq \
    --fq2 ${outdir}/sample_2.fastq \
    ${quiet} --logfile ${outdir}/haphpipe.out \
    --outdir ${outdir}

hp_reads ec_reads \
    --fq1 ${outdir}/trimmed_1.fastq \
    --fq2 ${outdir}/trimmed_2.fastq \
    ${quiet} --logfile ${outdir}/haphpipe.out \
    --outdir ${outdir}

hp_assemble_denovo \
    --fq1 ${outdir}/trimmed_1.fastq \
    --fq2 ${outdir}/trimmed_2.fastq \
    ${quiet} --logfile ${outdir}/haphpipe.out \
    --outdir ${outdir}

hp_assemble_amplicons \
    --contigs_fa ${outdir}/denovo_contigs.fna \
    --ref_fa ${refFA} \
    --ref_gtf ${refGTF} \
    ${quiet} --logfile ${outdir}/haphpipe.out \
    --outdir ${outdir}

hp_refine_assembly \
    --fq1 ${outdir}/corrected_1.fastq \
    --fq2 ${outdir}/corrected_2.fastq \
    --ref_fa ${outdir}/amplicon_assembly.fna \
    --max_step 2 \
    ${quiet} --logfile ${outdir}/haphpipe.out \
    --outdir ${outdir}

hp_finalize_assembly \
    --fq1 ${outdir}/corrected_1.fastq \
    --fq2 ${outdir}/corrected_2.fastq \
    --sample_id PL01 \
    --ref_fa ${outdir}/refined.fna \
    ${quiet} --logfile ${outdir}/haphpipe.out \
    --outdir ${outdir}