#!/bin/bash
#SBATCH --job-name=get_info
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=128GB
#SBATCH --ntasks=8
#SBATCH --time=2:0:0
#SBATCH --output=../process_reports/get_info_%A_%a
#SBATCH --array=1-22

source ./../config
cd ${cohort_results_path}
for i in `ls */chr${SLURM_ARRAY_TASK_ID}/*output`;do echo $i >> chr${SLURM_ARRAY_TASK_ID}_info;done

${Rscript} ${script_path}/01_DataProcess/scripts/01b.process_info_file.R ${SLURM_ARRAY_TASK_ID} ${vQTL_path}

sort chr${SLURM_ARRAY_TASK_ID}_info_full | uniq > chr${SLURM_ARRAY_TASK_ID}_info_full_uniq
mv chr${SLURM_ARRAY_TASK_ID}_info_full_uniq chr${SLURM_ARRAY_TASK_ID}_info_full
