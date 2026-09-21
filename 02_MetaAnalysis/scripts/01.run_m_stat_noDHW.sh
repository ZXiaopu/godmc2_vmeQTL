#!/bin/bash
#SBATCH --job-name=m_stat
#SBATCH --partition interruptible_cpu
#SBATCH --mem=256GB
#SBATCH --ntasks=8
#SBATCH --time=6:0:0

source ../../config

${Rscript} 01.run_m_stat_noDHW.R
