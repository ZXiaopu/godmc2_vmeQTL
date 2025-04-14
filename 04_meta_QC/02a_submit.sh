#!/bin/bash -l
#SBATCH --job-name=GoDMC_02a
#SBATCH --output=job_reports/GoDMC_02a_%j
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=12:0:0

bash 02a_LDclumping_cojo.sh ${SLURM_ARRAY_TASK_ID} ${1}
