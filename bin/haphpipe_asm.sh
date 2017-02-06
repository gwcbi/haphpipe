#! /bin/bash

SN="haphpipe_asm.sh"
read -r -d '' USAGE <<EOF
$SN
Viral genome assembly from fastq files

EOF

#--- Read command line args
[[ -n "$1" ]] && [[ "$1" == '-h' ]] && echo $USAGE && exit 0

[[ -n "$1" ]] && samp="$1"
[[ -n "$2" ]] && ref="$2"
[[ -n "$3" ]] && adapters="$3" || adapters=""

module unload python
module load miniconda3
source activate haphpipe

[[ -e /scratch ]] && export TMPDIR=/scratch

# vppath=/lustre/groups/cbi/Projects/viralpop.git/viralpop
# export PYTHONPATH=$vppath:$PYTHONPATH
# rroot=$(git rev-parse --show-toplevel)

echo "[---$SN---] ($(date)) Starting $SN"

#--- Check that fastq files exist
[[ ! -e "$samp/00_raw/original_1.fastq" ]] && echo "[---$SN---] ($(date)) FAILED: file $samp/00_raw/original_1.fastq does not exist" && exit 1
[[ ! -e "$samp/00_raw/original_2.fastq" ]] && echo "[---$SN---] ($(date)) FAILED: file $samp/00_raw/original_2.fastq does not exist" && exit 1
echo "[---$SN---] ($(date)) Sample:    $samp"

#--- Check that reference exists
[[ ! -e "$ref" ]] && echo "[---$SN---] ($(date)) FAILED: reference index $ref does not exist" && exit 1
echo "[---$SN---] ($(date)) Reference: $ref"

#--- Check adapters if provided
if [[ -n $adapters ]];
    [[ ! -e "$adapters" ]] && echo "[---$SN---] ($(date)) FAILED: adapters $adapters does not exist" && exit 1
    echo "[---$SN---] ($(date)) Adapters:  $adapters"
    aparam="--adapter_file $adapters"
else
    aparam=""
fi

#--- Start the timer
t1=$(date +"%s")

##########################################################################################
# Step 1: Trim reads.
##########################################################################################
echo "[---$SN---] ($(date)) Stage: trim_reads"
mkdir -p $samp/00_trim

hp_assemble trim_reads --ncpu $(nproc) \
    --adapter_file $aparam \
    --fq1 $samp/00_raw/original_1.fastq \
    --fq2 $samp/00_raw/original_2.fastq \
    --outdir $samp/00_trim

cd $samp/00_trim && ln -s trimmed_1.fastq read_1.fq && ln -s trimmed_2.fastq read_2.fq && cd ../..

##########################################################################################
# Step 2: Join reads with FLASh.
##########################################################################################
echo "[---$SN---] ($(date)) Stage: join_reads"
mkdir -p $samp/00_flash

hp_assemble join_reads --ncpu $(nproc) \
    --fq1 $samp/00_trim/read_1.fq \
    --fq2 $samp/00_trim/read_2.fq \
    --max_overlap 150 \
    --outdir $samp/00_flash

cd $samp/00_flash && ln -s flash.extendedFrags.fastq read_U.fq && cd ../..

##########################################################################################
# Step 3: Error correction using BLESS2
##########################################################################################
echo "[---$SN---] ($(date)) Stage: ec_reads"
mkdir -p $samp/00_bless

# bless is not currently available in bioconda
module load bless

hp_assemble ec_reads --ncpu $(nproc) \
    --fq1 $samp/00_trim/read_1.fq \
    --fq2 $samp/00_trim/read_2.fq \
    --kmerlength 31 \
    --outdir $samp/00_bless

cd $samp/00_bless && ln -s bless.1.corrected.fastq read_1.fq && ln -s bless.2.corrected.fastq read_2.fq && cd ../..

##########################################################################################
# Step 4: Denovo assembly using Trinity
##########################################################################################

for method in "trim" "bless" "flash"; do
    echo "[---$SN---] ($(date)) Using method $method"
    
    f1arg=$([[ -e $samp/00_${method}/read_1.fq ]] && echo "--fq1 $samp/00_${method}/read_1.fq" || echo "")
    f2arg=$([[ -e $samp/00_${method}/read_2.fq ]] && echo "--fq2 $samp/00_${method}/read_2.fq" || echo "")
    fUarg=$([[ -e $samp/00_${method}/read_U.fq ]] && echo "--fqU $samp/00_${method}/read_U.fq" || echo "")

    echo "[---$SN---] ($(date)) Stage: assemble_denovo $method"
    mkdir -p $samp/04_assembly/$method
        
    hp_assemble assemble_denovo \
        $f1arg $f2arg $fUarg \
        --ncpu $(( $(nproc) - 3 )) --max_memory 50 \
        --outdir $samp/04_assembly/$method

##########################################################################################
# Step 5: Assign contigs to subtypes
##########################################################################################
    echo "[---$SN---] ($(date)) Stage: subtype"
    echo "[---$SN---] ($(date)) Skipping subtyping stage"
    # Skip this part now

##########################################################################################
# Step 6: Scaffold contigs
##########################################################################################
    echo "[---$SN---] ($(date)) Stage: assemble_scaffold $method"
    mkdir -p $samp/06_scaffold/$method
    
    hp_assemble assemble_scaffold --keep_tmp \
            --contigs_fa $samp/04_assembly/$method/contigs.fa \
            --ref_fa $rroot/HIV/references/subtypes/HIV_B.K03455.HXB2_LAI_IIIB_BRU.fasta \
            --seqname $samp \
            --outdir $samp/06_scaffold/$method

##########################################################################################
# Step 7: Refine alignment 1
##########################################################################################
    echo "[---$SN---] ($(date)) Stage: refine_alignment (1) $method"
    mkdir -p $samp/07_refine1/$method
    
    # Use the imputed version
    cp $samp/06_scaffold/$method/imputed.fa $samp/07_refine1/$method/initial.fa
    
    hp_assemble refine_assembly --ncpu $(nproc) \
        $f1arg $f2arg $fUarg \
        --assembly_fa $samp/07_refine1/$method/initial.fa \
        --rgid $samp \
        --min_dp 1 \
        --bt2_preset very-fast \
        --outdir $samp/07_refine1/$method

##########################################################################################
# Step 8: Refine alignment 2
##########################################################################################
    echo "[---$SN---] ($(date)) Stage: refine_alignment (2) $method"
    mkdir -p $samp/08_refine2/$method
    
    # Impute refined assembly using first assembly
    # hp_assemble impute_ref \
    #     --assembly_fa $samp/07_refine1/$method/refined.fa \
    #     --ref_fa $samp/07_refine1/$method/initial.fa \
    #         > $samp/08_refine2/$method/initial.fa
    # Or not
    cp $samp/07_refine1/$method/refined.fa $samp/08_refine2/$method/initial.fa
    
    hp_assemble refine_assembly --ncpu $(nproc) \
        $f1arg $f2arg $fUarg \
        --assembly_fa $samp/08_refine2/$method/initial.fa \
        --rgid $samp \
        --min_dp 2 \
        --bt2_preset very-sensitive \
        --outdir $samp/08_refine2/$method

##########################################################################################
# Step 9: Fix consensus
##########################################################################################
    echo "[---$SN---] ($(date)) Stage: fix_consensus $method"
    mkdir -p $samp/09_fixed/$method
    
    hp_assemble fix_consensus --ncpu $(nproc) \
        $f1arg $f2arg $fUarg \
        --assembly_fa $samp/08_refine2/$method/refined.fa \
        --ref_fa $rroot/HIV/references/subtypes/HIV_B.K03455.HXB2_LAI_IIIB_BRU.fasta \
        --rgid $samp \
        --bt2_preset very-sensitive \
        --outdir $samp/09_fixed/$method

done

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
