#!/bin/bash
#SBATCH --job-name=process_07_res
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=128GB
#SBATCH --ntasks=8
#SBATCH --time=2:0:0
#SBATCH --output=../process_reports/process_07_res_%A_%a
#SBATCH --array=22

# $1 is the cohort name
echo $1
vQTL_path="/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/Stage2/vQTL_results/Cohorts_with_all_results"
bash convert_besd.sh ${vQTL_path}/$1/chr${SLURM_ARRAY_TASK_ID}
echo "successfully convert besd files to txt file format of chr"${SLURM_ARRAY_TASK_ID}

#/users/k20125144/.conda/envs/R/bin/Rscript combine_summary_stats.R ${vQTL_path}/$1/chr${SLURM_ARRAY_TASK_ID}
#echo "successfully combine results by three different methods of chr"${SLURM_ARRAY_TASK_ID}
