#!/bin/bash
#SBATCH --job-name=meta_check
#SBATCH --partition interruptible_cpu,cpu
#SBATCH --mem=32GB
#SBATCH --ntasks=1
#SBATCH --time=1:0:0

exec &> >(tee ../process_reports/meta_res_check_chr${1}_index${2}.out)

cd ../data/Meta_Results/chr$1

chr=$1
index=$2

BF=`wc -l BF_random_meta_results_cpg_index${index}1.txt`
DRM=`wc -l DRM_random_meta_results_cpg_index${index}1.txt`
SVLM=`wc -l SVLM_random_meta_results_cpg_index${index}1.txt`

BFline=`echo $BF | cut -f 1 -d ' '`
DRMline=`echo $DRM | cut -f 1 -d ' '`
SVLMline=`echo $SVLM | cut -f 1 -d ' '`

BFcpg=`cut -f 1 BF_random_meta_results_cpg_index${index}1.txt  | cut -f 4 -d '_' | sort | uniq | wc -l`
DRMcpg=`cut -f 1 DRM_random_meta_results_cpg_index${index}1.txt  | cut -f 4 -d '_' | sort | uniq | wc -l`
SVLMcpg=`cut -f 1 SVLM_random_meta_results_cpg_index${index}1.txt  | cut -f 4 -d '_' | sort | uniq | wc -l`

echo "meta-analysis in chr$chr index${index}: line - $BFline; cpg - $BFcpg"
echo "meta-analysis in chr$chr index${index}: line - $DRMline; cpg - $DRMcpg"
echo "meta-analysis in chr$chr index${index}: line - $SVLMline; cpg - $SVLMcpg"
