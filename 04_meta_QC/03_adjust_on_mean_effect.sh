#!/bin/bash

source config

echo "Extract genetic variants remained after LD clumping and conditional independent analysis"

exec &> >(tee ${vmeQTL_03_logfile}_chr${1})

mkdir -p ${processed_data_dir}/adjust_mean_effect
i=${1}

${plink} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
         --extract ${processed_data_dir}/LD_COJO_input/chr${i}/vmeQTL_vCpG_pair_after_LD_COJO.snplist \
         --recode A-transpose \
         --out ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.filtered

${R_directory}/Rscript ${vmeQTL_script_dir}/adjust_mean_effects.R \
        ${untransformed_methylation_adjusted_pcs}.RData \
        ${processed_data_dir}/LD_COJO_input/chr${i}/vmeQTL_vCpG_pair_after_LD_COJO.txt \
        ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.filtered.traw \
        ${processed_data_dir}/adjust_mean_effect/vmeQTL_vCpG_adjust_mean_effect_chr${i}.txt

echo "successfully run 03 adjust addtive genetic effects of chr$i"
