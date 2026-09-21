#!/bin/bash
#SBATCH --job-name=filter_BF_only
#SBATCH --partition cpu,interruptible_cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=2
#SBATCH --time=1:0:0

exec &> >(tee ../process_reports/filter_BF_only_asso_${SLURM_ARRAY_TASK_ID}.out)

info=`head -n ${SLURM_ARRAY_TASK_ID} ../data/Meta_file_summary_index.txt | tail -n1`
chr=`echo $info | cut -f 1 -d ' '`
index=`echo $info | cut -f 2 -d ' '`

echo $chr
echo $index

conda activate R
${Rscript} 07a.filter_BF_only_asso_new.R $chr $index
echo "association filtration for chr$chr index$index has been done"
