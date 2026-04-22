#!/bin/bash
#SBATCH --job-name=run_m_stat
#SBATCH --partition interruptible_cpu
#SBATCH --mem=256GB
#SBATCH --ntasks=8
#SBATCH --time=6:0:0

for c in $(seq 2 22);
do
        sbatch 00.get_mstat_list.sh ${c} BF
        sbatch 00.get_mstat_list.sh ${c} DRM
        sbatch 00.get_mstat_list.sh ${c} SVLM
done
