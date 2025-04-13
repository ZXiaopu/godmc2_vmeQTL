#!/bin/bash -l
#SBATCH --job-name=GoDMC_02
#SBATCH --output=job_reports/GoDMC_02_%j
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=12:0:0

sbatch 02_LDclumping_cojo.sh ${SLURM_ARRAY_TASK_ID} BF
sbatch 02_LDclumping_cojo.sh ${SLURM_ARRAY_TASK_ID} SVLM
sbatch 02_LDclumping_cojo.sh ${SLURM_ARRAY_TASK_ID} DRM
