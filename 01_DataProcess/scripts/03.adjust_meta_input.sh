#!/bin/bash
#SBATCH --job-name=adjust_meta_input
#SBATCH --partition cpu,interruptible_cpu
#SBATCH --mem=32GB
#SBATCH --ntasks=1
#SBATCH --time=2:0:0

exec &> >(tee ../process_reports/adjust_meta_chr${4}_index${3}_${1}_${2}.out)

source ../../config

method=$1
cohort=$2
index=$3
chr=$4

home_dir=${meta_input_path}
file="${method}_${cohort}_cpg_index${index}_meta_input.txt"

cd ${home_dir}
mkdir -p chr${chr}
cat meta_input_title ${method}_${cohort}_*chr${chr}_cpg_index${index}_meta_input.txt > chr${chr}/${file}
echo "adjust meta-analysis for $file has been successful"
