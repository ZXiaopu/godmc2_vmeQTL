#!/bin/bash
#SBATCH --job-name=process_07_res
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=128GB
#SBATCH --ntasks=8
#SBATCH --time=2:0:0
#SBATCH --array=6

# $1 is the cohort name
exec &> >(tee ../process_reports/process_07_${1}_chr${SLURM_ARRAY_TASK_ID}.out)

source ../../config
echo $1
bash convert_besd.sh ${vQTL_path}/$1/chr${SLURM_ARRAY_TASK_ID}
echo "successfully convert besd files to txt file format of chr"${SLURM_ARRAY_TASK_ID}

