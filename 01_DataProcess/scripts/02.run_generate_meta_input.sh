#!/bin/bash
#SBATCH --job-name=generate_meta_input
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=16GB
#SBATCH --ntasks=8
#SBATCH --time=8:0:0
#SBATCH --output=../process_reports/generate_meta_input_%A_%a

cohort=$1
cpg_chunk=22
g_chunk=${SLURM_ARRAY_TASK_ID}
input_path="/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/Stage2/vQTL_results/Cohorts_with_all_results/$cohort/chr$cpg_chunk"
home_path="/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/01_DataProcess"

/users/k20125144/.conda/envs/R/bin/Rscript 02.generate_meta_input.R ${input_path} ${cohort} ${g_chunk} ${cpg_chunk} ${home_path}
echo "Meta-analysis input for $cohort has been done"
