#!/bin/bash -l
#SBATCH --job-name=GoDMC_02b
#SBATCH --output=job_reports/GoDMC_02b_%j
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=12:0:0

bash 02_LDclumping_cojo.sh ${SLURM_ARRAY_TASK_ID}
