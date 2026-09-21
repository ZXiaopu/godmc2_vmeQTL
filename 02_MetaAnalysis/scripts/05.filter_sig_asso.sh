#!/bin/bash
#SBATCH --job-name=filter_sig_meta
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --time=12:0:0
#SBATCH --mem=32G
#SBATCH --output=filter_sig_meta.%A.%a

file=`head -n ${SLURM_ARRAY_TASK_ID} ../data/Meta_file_summary.txt | tail -n1`
file1="${file%.txt}"

cd ../data/
awk 'BEGIN{FS="\t"}{if($10<5e-8) print $0}' $file > ${file1}_fixed_meta_pval5e-8.txt
cat ./Meta_Results/title ${file1}_fixed_meta_pval5e-8.txt > ${file1}_fixed_meta_pval5e-8.txt1
mv ${file1}_fixed_meta_pval5e-8.txt1 ${file1}_fixed_meta_pval5e-8.txt

echo "filtering significant associations has been done in ${file}"
