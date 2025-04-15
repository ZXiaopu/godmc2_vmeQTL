#!/bin/bash -l
#SBATCH --job-name=GoDMC_04
#SBATCH --output=job_reports/GoDMC_04_%j
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=2:0:0

bash 04_run_interaction_analysis.sh ${SLURM_ARRAY_TASK_ID}
