#! /bin/bash

# Copyright (C) 2019 Matthew L. Bendall and (C) 2020 Keylie M. Gibson

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

SN="haphpipe_predicthaplo.sh"

if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
    # This is an array job
    [[ ! -n "$samplist" ]] && echo "Array job requires sample list as environment variable \"\$samplist\"" && exit 1
    samp=$(sed -n "$SLURM_ARRAY_TASK_ID"p $samplist)
else
    [[ -n "$1" ]] && samp="$1"
fi

[[ -z "$samp" ]] && echo "\"\$samp\" was not provided" && exit 1
[[ ! -d "$samp" ]] && echo "Directory \"$samp\" was not found" && exit 1

module unload python
module load miniconda3
source activate haphpipe

[[ -e /scratch ]] && export TMPDIR=/scratch

module load PredictHaplo/0.4
module load samtools/1.3.1

echo "[---$SN---] ($(date)) Starting $SN"

#--- Start the timer
t1=$(date +"%s")

for method in "trim" "bless"; do
    echo "[---$SN---] ($(date)) Stage: PredictHaplo, $method"
    mkdir -p $samp/10_predicthaplo/${method}
    
    [[ -e $samp/09_fixed/${method}/covered.intervals ]] \
      && intargs="--interval_txt $samp/09_fixed/${method}/covered.intervals" \
      || intargs="--min_depth 20"
    
    hp_haplotype predict_haplo --ncpu $(nproc) \
        --alignment $samp/09_fixed/${method}/final.bam \
        --ref_fa $samp/09_fixed/${method}/consensus.fasta \
        --outdir $samp/10_predicthaplo/${method} \
        $intargs \
        --min_interval 200 &
done

# Wait for all processes to finsh
wait

#---Complete job
t2=$(date +"%s")
diff=$(($t2-$t1))
echo "[---$SN---] ($(date)) $(($diff / 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "[---$SN---] ($(date)) $SN COMPLETE."
