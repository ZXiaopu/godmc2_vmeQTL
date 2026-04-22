#!/bin/bash

source ../../config

cpg_index=$1
chr=$2
cohort=`cat ${cohort_results_path}/chr${chr}_info_full | cut -f 1 -d ',' | sort | uniq`

for c in $cohort;
do
sbatch 03.adjust_meta_input.sh BF $c $1 $chr
sbatch 03.adjust_meta_input.sh SVLM $c $1 $chr
sbatch 03.adjust_meta_input.sh DRM $c $1 $chr
done
