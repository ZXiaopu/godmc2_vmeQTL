#!/bin/bash -l
#SBATCH --job-name=GoDMC_03
#SBATCH --output=job_reports/GoDMC_03_%j
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=2:0:0

bash 03_adjust_on_mean_effect.sh ${SLURM_ARRAY_TASK_ID}
