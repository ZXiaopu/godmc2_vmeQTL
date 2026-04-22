#!/bin/bash
#SBATCH --job-name=generate_meta_input
#SBATCH --partition cpu,interruptible_cpu
#SBATCH --mem=256GB
#SBATCH --ntasks=8
#SBATCH --time=12:0:0

exec &> >(tee ../process_reports/generate_meta_chr${1}_${SLURM_ARRAY_TASK_ID}.out)

source ../../config

chr=$1
info=`head -n ${SLURM_ARRAY_TASK_ID} ${vQTL_path}/chr${chr}_info_full | tail -n1`
cohort=`echo $info | cut -d ',' -f 1`
cpg_chr=`echo $info | cut -d ',' -f 2`
g_chunk=`echo $info | cut -d ',' -f 3`
sample_size=`echo $info | cut -d ',' -f 4`

input_path="${vQTL_path}/$cohort/$cpg_chr"
home_path=${meta_input_path}

${Rscript} 02.generate_meta_input.R ${input_path} ${cohort} ${g_chunk} ${cpg_chr} ${home_path} ${sample_size}

echo "Meta-analysis input for $cohort has been done"
