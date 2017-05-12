#! /bin/bash

SN="haphpipe_amp.sh"
read -r -d '' USAGE <<EOF
$SN [sample_dir] [reference_fasta] [reference_gtf] <adapters_fasta>

Viral amplicon assembly from fastq files

EOF

#--- Read command line args
[[ -n "$1" ]] && [[ "$1" == '-h' ]] && echo "$USAGE" && exit 0

[[ -n "$1" ]] && samp="$1"
[[ -n "$2" ]] && ref="$2"
[[ -n "$3" ]] && refgtf="$3"
[[ -n "$4" ]] && adapters="$4"

[[ -z ${samp+x} ]] && echo "FAILED: sample_dir is not set" && echo "$USAGE" && exit 0
[[ -z ${ref+x} ]] && echo "FAILED: reference_fasta is not set" && echo "$USAGE" && exit 0
[[ -z ${refgtf+x} ]] && echo "FAILED: ref_gtf is not set" && echo "$USAGE" && exit 0

module unload python
module load miniconda3
source activate haphpipe

[[ -e /scratch ]] && export TMPDIR=/scratch
# MAXPROC=$(nproc)
MAXPROC=8

echo "[---$SN---] ($(date)) Starting $SN"

#--- Check that fastq files exist
[[ ! -e "$samp/00_raw/original_1.fastq" ]] &&\
 echo "[---$SN---] ($(date)) FAILED: file $samp/00_raw/original_1.fastq does not exist" &&\
 exit 1

[[ ! -e "$samp/00_raw/original_2.fastq" ]] &&\
 echo "[---$SN---] ($(date)) FAILED: file $samp/00_raw/original_2.fastq does not exist" &&\
 exit 1

echo "[---$SN---] ($(date)) Sample:    $samp"

#--- Check that reference exists
[[ ! -e "$ref" ]] && echo "[---$SN---] ($(date)) FAILED: reference index $ref does not exist" && exit 1
echo "[---$SN---] ($(date)) Reference: $ref"

#--- Check adapters if provided
if [[ -n "$adapters" ]]; then
    [[ ! -e "$adapters" ]] && echo "[---$SN---] ($(date)) FAILED: adapters $adapters does not exist" && exit 1
    echo "[---$SN---] ($(date)) Adapters:  $adapters"
    aparam="--adapter_file $adapters"
else
    echo "[---$SN---] ($(date)) Adapters:  $adapters"
    aparam=""
fi

echo "[---$SN---] ($(date)) Ad param:  $aparam"

#--- Start the timer
t1=$(date +"%s")

##########################################################################################
# Step 1: Trim reads.
##########################################################################################
stage="trim_reads"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/00_trim

if [[ -e $samp/00_trim/read_1.fq && -e $samp/00_trim/read_2.fq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage read_1.fq,read_2.fq"
else
    read -r -d '' cmd <<EOF
hp_assemble trim_reads --ncpu $MAXPROC\
 $aparam\
 --fq1 $samp/00_raw/original_1.fastq\
 --fq2 $samp/00_raw/original_2.fastq\
 --outdir $samp/00_trim"
EOF

    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd
    
    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
    
    # Symlink to rename files
    cd $samp/00_trim &&\
      ln -fs trimmed_1.fastq read_1.fq && ln -fs trimmed_2.fastq read_2.fq &&\
      cd ../..
fi

##########################################################################################
# Step 2a: Error correction using BLESS2
##########################################################################################
stage="ec_reads"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/00_bless

if [[ -e $samp/00_bless/read_1.fq && -e $samp/00_bless/read_2.fq ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage read_1.fq,read_2.fq"
else
    # bless is not currently available in bioconda
    module load bless
    
    read -r -d '' cmd <<EOF
hp_assemble ec_reads --ncpu $MAXPROC\
 --fq1 $samp/00_trim/read_1.fq\
 --fq2 $samp/00_trim/read_2.fq\
 --kmerlength 31\
 --outdir $samp/00_bless
EOF
    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd
    
    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
    
    # Symlink to rename files
    cd $samp/00_bless &&\
      ln -fs bless.1.corrected.fastq read_1.fq && ln -fs bless.2.corrected.fastq read_2.fq &&\
      cd ../..
fi

# Calculate number of error corrected reads
numec=$(( $(wc -l < $samp/00_bless/read_1.fq) / 4 ))
echo -e "[---$SN---] ($(date)) Number of error corrected reads: $numec"

##########################################################################################
# Step 2b: Convert FASTQ to unaligned BAM file
##########################################################################################
stage="fq_to_bam"
echo "[---$SN---] ($(date)) Stage: $stage"

if [[ -e $samp/00_bless/reads.bam ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage reads.bam"
else
    read -r -d '' cmd <<EOF
picard FastqToSam SM=$samp RG=$samp\
 F1=$samp/00_bless/read_1.fq F2=$samp/00_bless/read_2.fq\
 O=$samp/00_bless/reads.bam
EOF

    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd
    
    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

##########################################################################################
# Step 4: Denovo assembly using Trinity
##########################################################################################
stage="assemble_denovo"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/04_assembly

if [[ -e $samp/04_assembly/contigs.fa ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage contigs.fa"
else
    # Subsample reads - 50 percent, seed=9
    samtools view -bs 9.5 $samp/00_bless/reads.bam > $samp/04_assembly/sub.bam
    picard SamToFastq I=$samp/04_assembly/sub.bam F=$samp/04_assembly/sub_1.fq F2=$samp/04_assembly/sub_2.fq

    read -r -d '' cmd <<EOF
hp_assemble assemble_denovo\
 --fq1 $samp/04_assembly/sub_1.fq\
 --fq2 $samp/04_assembly/sub_2.fq\
 --ncpu $(( $MAXPROC - 3 )) --max_memory 50\
 --outdir $samp/04_assembly
EOF

    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

rm -f $samp/04_assembly/sub.bam
rm -f $samp/04_assembly/sub_1.fq
rm -f $samp/04_assembly/sub_2.fq

##########################################################################################
# Step 5: Assign contigs to subtypes
##########################################################################################
# stage="subtype"
# echo "[---$SN---] ($(date)) Stage: $stage"
# echo "[---$SN---] ($(date)) Skipping $stage"

##########################################################################################
# Step 6: Assemble contigs to amplicons
##########################################################################################
stage="assemble_amplicons"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/06_amplicons

if [[ -e $samp/06_amplicons/assembly.fa ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage assembly.fa"
else
    read -r -d '' cmd <<EOF 
hp_assemble assemble_amplicons\
 --contigs_fa $samp/04_assembly/contigs.fa\
 --ref_fa $ref\
 --ref_gtf $refgtf\
 --seqname $samp \
 --outdir $samp/06_amplicons
EOF

    echo -e "[---$SN---] ($(date)) $stage command:\n\n$cmd\n"
    eval $cmd
    
    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
fi

##########################################################################################
# Step 7: Refine alignment 1
##########################################################################################
stage="refine_alignment_1"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/07_refine1

if [[ -e $samp/07_refine1/refined.fasta ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage refined.fasta"
else
    # Subsample reads - 20 percent, seed=1
    samtools view -bs 1.2 $samp/00_bless/reads.bam > $samp/07_refine1/sub.bam
    picard SamToFastq I=$samp/07_refine1/sub.bam F=$samp/07_refine1/sub_1.fq F2=$samp/07_refine1/sub_2.fq

    # Use the amplicon assembly
    cp $samp/06_amplicons/assembly.fa $samp/07_refine1/initial.fa
    
    read -r -d '' cmd <<EOF 
hp_assemble align_reads --ncpu $MAXPROC\
 --ref_fa $samp/07_refine1/initial.fa\
 --fq1 $samp/07_refine1/sub_1.fq\
 --fq2 $samp/07_refine1/sub_2.fq\
 --rgid $samp\
 --bt2_preset fast-local\
 --outdir $samp/07_refine1
EOF

    echo -e "[---$SN---] ($(date)) $stage align_reads command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage align_reads" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )

    read -r -d '' cmd <<EOF 
hp_assemble call_variants --ncpu $MAXPROC\
 --ref_fa $samp/07_refine1/initial.fa\
 --aln_bam $samp/07_refine1/aligned.bam\
 --emit_all\
 --outdir $samp/07_refine1
EOF

    echo -e "[---$SN---] ($(date)) $stage call_variants command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage call_variants" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )

    module load viral-ngs

    read -r -d '' cmd <<EOF
assembly.py vcf_to_fasta\
 --min_coverage 1 \
 $samp/07_refine1/variants.vcf.gz\
 $samp/07_refine1/refined.fasta
EOF

    echo -e "[---$SN---] ($(date)) $stage vcf_to_fasta command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage vcf_to_fasta" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )

    module unload viral-ngs
fi

rm -f $samp/07_refine1/sub.bam
rm -f $samp/07_refine1/sub_1.fq
rm -f $samp/07_refine1/sub_2.fq
rm -f $samp/07_refine1/aligned.bam
rm -f $samp/07_refine1/aligned.bam.bai

##########################################################################################
# Step 8: Refine alignment 2
##########################################################################################
stage="refine_alignment_2"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/08_refine2

if [[ -e $samp/08_refine2/refined.fasta ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage refined.fasta"
else
    # Subsample reads - 50 percent, seed=3
    samtools view -bs 3.5 $samp/00_bless/reads.bam > $samp/08_refine2/sub.bam
    picard SamToFastq I=$samp/08_refine2/sub.bam F=$samp/08_refine2/sub_1.fq F2=$samp/08_refine2/sub_2.fq
    
    # Use the refined assembly
    cp $samp/07_refine1/refined.fasta $samp/08_refine2/initial.fa

    read -r -d '' cmd <<EOF 
hp_assemble align_reads --ncpu $MAXPROC\
 --ref_fa $samp/08_refine2/initial.fa\
 --fq1 $samp/08_refine2/sub_1.fq\
 --fq2 $samp/08_refine2/sub_2.fq\
 --rgid $samp\
 --bt2_preset very-sensitive-local\
 --outdir $samp/08_refine2
EOF

    echo -e "[---$SN---] ($(date)) $stage align_reads command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage align_reads" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )

    read -r -d '' cmd <<EOF 
hp_assemble call_variants --ncpu $MAXPROC\
 --ref_fa $samp/08_refine2/initial.fa\
 --aln_bam $samp/08_refine2/aligned.bam\
 --emit_all\
 --outdir $samp/08_refine2
EOF

    echo -e "[---$SN---] ($(date)) $stage call_variants command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage call_variants" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )

    module load viral-ngs

    read -r -d '' cmd <<EOF
assembly.py vcf_to_fasta\
 --min_coverage 1 \
 $samp/08_refine2/variants.vcf.gz\
 $samp/08_refine2/refined.fasta
EOF

    echo -e "[---$SN---] ($(date)) $stage vcf_to_fasta command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage vcf_to_fasta" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )

    module unload viral-ngs
fi

rm -f $samp/08_refine2/sub.bam
rm -f $samp/08_refine2/sub_1.fq
rm -f $samp/08_refine2/sub_2.fq
rm -f $samp/08_refine2/aligned.bam
rm -f $samp/08_refine2/aligned.bam.bai

##########################################################################################
# Step 9a: Fix consensus - pairwise_align
##########################################################################################
stage="fix_consensus"
echo "[---$SN---] ($(date)) Stage: $stage"
mkdir -p $samp/09_fixed

if [[ -e $samp/09_fixed/consensus.fa ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage consensus.fa"
else
    # Use the refined assembly
    cp $samp/08_refine2/refined.fasta $samp/09_fixed/amplicons.fa
    
    module load blast+
    
    # Align refined assembly to reference
    read -r -d '' cmd <<EOF
hp_assemble pairwise_align \
 --amplicons_fa $samp/09_fixed/amplicons.fa \
 --ref_fa $ref \
 --ref_gtf $refgtf \
 --outdir $samp/09_fixed
EOF

    echo -e "[---$SN---] ($(date)) $stage pairwise_align command:\n\n$cmd\n"
    eval $cmd

    # Extract padded consensus
    hp_assemble extract_pairwise \
        --align_json $samp/09_fixed/alignments.json \
        --outfmt nuc_fa > $samp/09_fixed/consensus.fa
    
    # Extract GTF with amplicon regions
    hp_assemble extract_pairwise \
        --align_json $samp/09_fixed/alignments.json \
        --outfmt amp_gtf > $samp/09_fixed/consensus.gtf
    
    # Extract full alignment
    hp_assemble extract_pairwise \
        --align_json $samp/09_fixed/alignments.json \
        --outfmt aln_fa > $samp/09_fixed/refalign.fa

fi

##########################################################################################
# Step 9b: Fix consensus - align_reads
##########################################################################################
if [[ -e $samp/09_fixed/final.bam ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage final.bam"
else
    read -r -d '' cmd <<EOF
hp_assemble align_reads --ncpu $MAXPROC\
 --ref_fa $samp/09_fixed/consensus.fa\
 --fq1 $samp/00_bless/read_1.fq\
 --fq2 $samp/00_bless/read_2.fq\
 --rgid $samp\
 --bt2_preset very-sensitive-local\
 --no_markdup\
 --outdir $samp/09_fixed
EOF

    echo -e "[---$SN---] ($(date)) $stage align_reads command:\n\n$cmd\n"
    eval $cmd
    
    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage align_reads" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )

    # Filter the final alignment
    # Keep 2 (PROPER_PAIR)
    # Remove 3852 (UNMAP,MUNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY)
    samtools view -b -f 2 -F 3852 $samp/09_fixed/aligned.bam > $samp/09_fixed/final.bam
    samtools index $samp/09_fixed/final.bam
fi

##########################################################################################
# Step 9c: Fix consensus - call_variants
##########################################################################################
if [[ -e $samp/09_fixed/final.vcf.gz ]]; then
    echo "[---$SN---] ($(date)) EXISTS: $stage final.vcf.gz"
else
    # Call variants
    read -r -d '' cmd <<EOF 
hp_assemble call_variants --ncpu $MAXPROC\
 --ref_fa $samp/09_fixed/consensus.fa\
 --aln_bam $samp/09_fixed/final.bam\
 --outdir $samp/09_fixed
EOF

    echo -e "[---$SN---] ($(date)) $stage call_variants command:\n\n$cmd\n"
    eval $cmd

    [[ $? -eq 0 ]] && echo "[---$SN---] ($(date)) COMPLETED: $stage call_variants" || \
        (  echo "[---$SN---] ($(date)) FAILED: $stage" && exit 1 )
    
    mv $samp/09_fixed/variants.vcf.gz $samp/09_fixed/final.vcf.gz
    mv $samp/09_fixed/variants.vcf.gz.tbi $samp/09_fixed/final.vcf.gz.tbi
fi


##########################################################################################
# Step 10: Post-Assembly
##########################################################################################
mkdir -p $samp/10_summary

# Summary file
echo -e "SAMPID\t$samp" > $samp/10_summary/summary.txt
[[ -z ${numraw+x} ]] && numraw=$(( $(wc -l < $samp/00_raw/original_1.fastq) / 4 ))
echo -e "NREADS_RAW\t$numraw" >> $samp/10_summary/summary.txt
[[ -z ${numec+x} ]] && numec=$(( $(wc -l < $samp/00_bless/read_1.fq) / 4 ))
echo -e "NREADS_EC\t$numec" >> $samp/10_summary/summary.txt
[[ -z ${numfin+x} ]] && numfin=$(samtools view -f 66 -F 268 $samp/09_fixed/final.bam | wc -l)
echo -e "NREADS_FINAL\t$numfin" >> $samp/10_summary/summary.txt
alnrate=$(grep 'overall alignment rate' $samp/09_fixed/bowtie2.out | sed 's/% overall alignment rate//')
echo -e "ALN_RATE\t$alnrate" >> $samp/10_summary/summary.txt
cat $samp/06_amplicons/summary.txt >> $samp/10_summary/summary.txt

# Annotate sequence
cat $samp/09_fixed/consensus.gtf > $samp/10_summary/annotation.gtf

hp_assemble annotate_from_ref \
 --align_json $samp/09_fixed/alignments.json\
 --ref_gtf  ref/HIV_B.K03455.HXB2.gtf \
  >> $samp/10_summary/annotation.gtf

# Summarize targets
module load bcftools/1.3.1
hp_assemble post_assembly \
    --consensus_fa $samp/09_fixed/consensus.fa\
    --annotations_gtf $samp/10_summary/annotation.gtf\
    --align_bam $samp/09_fixed/aligned.bam\
    --variants_vcf $samp/09_fixed/final.vcf.gz > $samp/10_summary/region_summary.txt

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."

# 
# # Check files
# for method in "trim" "bless" "flash"; do
#     [[ -e $samp/09_fixed/$method/consensus.fasta ]] && \
#     [[ -e $samp/09_fixed/$method/final.bam ]] && \
#     [[ -e $samp/09_fixed/$method/variants.ug.vcf.gz ]] && echo "[---$SN---] ($(date)) $samp $method SUCCESS" || echo "[---$SN---] ($(date)) $samp $method FAILED"
# done | tee $samp/x.log
# 
# [[ $(grep -c 'FAILED' $samp/x.log)  == 0 ]] && echo "[---$SN---] ($(date)) $samp SUCCESS" || echo "[---$SN---] ($(date)) $samp FAILED"
# rm -f $samp/x.log
