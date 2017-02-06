#! /bin/bash

samp=$1

module unload python
module load miniconda3
source activate viralpop

vppath=/lustre/groups/cbi/Projects/viralpop.git/viralpop
export PYTHONPATH=$vppath:$PYTHONPATH
export TMPDIR=/scratch
rroot=$(git rev-parse --show-toplevel)

##########################################################################################
# Step 1: Trim reads.
##########################################################################################
mkdir -p $samp/00_trim

$vppath/assembly.py trim_reads --ncpu $(nproc) \
    --adapter_file $rroot/adapters/adapters.fasta \
    --fq1 $samp/00_raw/original_1.fastq \
    --fq2 $samp/00_raw/original_2.fastq \
    --outdir $samp/00_trim

cd $samp/00_trim && ln -s trimmed_1.fastq read_1.fq && ln -s trimmed_2.fastq read_2.fq && cd ../..

##########################################################################################
# Step 2: Join reads with FLASh.
##########################################################################################
mkdir -p $samp/00_flash

$vppath/assembly.py join_reads --ncpu $(nproc) \
    --fq1 $samp/00_trim/read_1.fq \
    --fq2 $samp/00_trim/read_2.fq \
    --max_overlap 150 \
    --outdir $samp/00_flash

cd $samp/00_flash && ln -s flash.extendedFrags.fastq read_U.fq && cd ../..

##########################################################################################
# Step 3: Error correction using BLESS2
##########################################################################################
mkdir -p $samp/00_bless

# bless is not currently available in bioconda
module load bless

$vppath/assembly.py ec_reads --ncpu $(nproc) \
    --fq1 $samp/00_trim/read_1.fq \
    --fq2 $samp/00_trim/read_2.fq \
    --kmerlength 31 \
    --outdir $samp/00_bless

cd $samp/00_bless && ln -s bless.1.corrected.fastq read_1.fq && ln -s bless.2.corrected.fastq read_2.fq && cd ../..

##########################################################################################
# Step 4: Denovo assembly using Trinity
##########################################################################################

for method in "bless" "flash"; do
    echo "Assembling $method"
    mkdir -p $samp/04_assembly/$method
    f1arg=$([[ -e $samp/00_${method}/read_1.fq ]] && echo "--fq1 $samp/00_${method}/read_1.fq" || echo "")
    f2arg=$([[ -e $samp/00_${method}/read_2.fq ]] && echo "--fq2 $samp/00_${method}/read_2.fq" || echo "")
    fUarg=$([[ -e $samp/00_${method}/read_U.fq ]] && echo "--fqU $samp/00_${method}/read_U.fq" || echo "")

    $vppath/assembly.py assemble_denovo \
        $f1arg $f2arg $fUarg \
        --ncpu $(( $(nproc) - 3 )) --max_memory 50 \
        --outdir $samp/04_assembly/$method

##########################################################################################
# Step 5: Assign contigs to subtypes
##########################################################################################
    echo "Skip subtyping"
# Skip this part now

##########################################################################################
# Step 6: Scaffold contigs
##########################################################################################
    echo "Scaffolding $method"
    mkdir -p $samp/06_scaffold/$method
    
    $vppath/assembly.py assemble_scaffold --keep_tmp \
            --contigs_fa $samp/04_assembly/$method/contigs.fa \
            --ref_fa $rroot/HIV/references/subtypes/HIV_B.K03455.HXB2_LAI_IIIB_BRU.fasta \
            --seqname $samp \
            --outdir $samp/06_scaffold/$method

##########################################################################################
# Step 7: Refine alignment 1
##########################################################################################
    echo "Refinment 1 $method"
    mkdir -p $samp/07_refine1/$method
    
    # Use the imputed version
    cp $samp/06_scaffold/$method/imputed.fa $samp/07_refine1/$method/initial.fa
    
    $vppath/assembly.py refine_assembly --ncpu $(nproc) \
        $f1arg $f2arg $fUarg \
        --assembly_fa $samp/07_refine1/$method/initial.fa \
        --rgid $samp \
        --min_dp 1 \
        --bt2_preset very-fast \
        --outdir $samp/07_refine1/$method

##########################################################################################
# Step 8: Refine alignment 2
##########################################################################################
    echo "Refinment 2 $method"
    mkdir -p $samp/08_refine2/$method
    
    # Impute refined assembly using first assembly
    # $vppath/assembly.py impute_ref \
    #     --assembly_fa $samp/07_refine1/$method/refined.fa \
    #     --ref_fa $samp/07_refine1/$method/initial.fa \
    #         > $samp/08_refine2/$method/initial.fa
    # Or not
    cp $samp/07_refine1/$method/refined.fa $samp/08_refine2/$method/initial.fa
    
    $vppath/assembly.py refine_assembly --ncpu $(nproc) \
        $f1arg $f2arg $fUarg \
        --assembly_fa $samp/08_refine2/$method/initial.fa \
        --rgid $samp \
        --min_dp 2 \
        --bt2_preset very-sensitive \
        --outdir $samp/08_refine2/$method

##########################################################################################
# Step 9: Fix consensus
##########################################################################################
    echo "Fix consensus $method"
    mkdir -p $samp/09_fixed/$method
    
    $vppath/assembly.py fix_consensus --ncpu $(nproc) \
        $f1arg $f2arg $fUarg \
        --assembly_fa $samp/08_refine2/$method/refined.fa \
        --ref_fa $rroot/HIV/references/subtypes/HIV_B.K03455.HXB2_LAI_IIIB_BRU.fasta \
        --rgid $samp \
        --bt2_preset very-sensitive \
        --outdir $samp/09_fixed/$method

done

# for method in "trim" "bless" "flash"; do

