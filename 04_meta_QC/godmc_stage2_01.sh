#!/bin/bash -l
#SBATCH --job-name=GoDMC_01
#SBATCH --output=job_reports/GoDMC_01_%j
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=2:0:0

bash 01_observe_vmeQTL_snps.sh ${SLURM_ARRAY_TASK_ID}
