#!/bin/bash
#SBATCH --job-name=random_meta
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=256GB
#SBATCH --time=12:0:0

source ../../config
cd ../data/Meta_Input_Data/chr${1}/${2}

${metal} ${3}
#mv ${1}_random_meta_results_cpg_index${2}* ../../../02_MetaAnalysis/data

#echo "meta-analysis of vmeQTLs on cpg index - method ${1} cpg index ${2} finished"
