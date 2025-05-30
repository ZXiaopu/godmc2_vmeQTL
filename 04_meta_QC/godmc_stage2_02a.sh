#!/bin/bash -l
#SBATCH --job-name=GoDMC_02a
#SBATCH --output=job_reports/GoDMC_02a_%j
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=12:0:0

sbatch --array 1-22 02a_submit.sh BF
sbatch --array 1-22 02a_submit.sh SVLM
sbatch --array 1-22 02a_submit.sh DRM
