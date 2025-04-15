#!/bin/bash

source config
exec &> >(tee ${vmeQTL_04_logfile}_chr${1})

mkdir -p ${processed_data_dir}/interaction_analysis

i=${1}

${R_directory}/Rscript ${vmeQTL_script_dir}/interaction_analysis.R \
        ${untransformed_methylation_adjusted_pcs}.RData \
        ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.filtered.traw \
        ${processed_data_dir}/adjust_mean_effect/vmeQTL_vCpG_adjust_mean_effect_chr${i}.txt \
        ${covariates_combined}.txt \
        ${processed_data_dir}/interaction_analysis \
        ${i}
