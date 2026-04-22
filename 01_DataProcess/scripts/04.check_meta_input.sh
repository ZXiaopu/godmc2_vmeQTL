#!/bin/bash
#SBATCH --job-name=check_meta_input
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --mem=16GB
#SBATCH --ntasks=1
#SBATCH --time=2:0:0

exec &> >(tee ../process_reports/check_meta_input_${1}_chr${2}_index${3}.out)

source ../../config

cohort=$1
chr=$2
index=$3

cd {meta_input_path}/data/chr${chr}
BF=`wc -l BF_${cohort}_cpg_index${index}_meta_input.txt`
DRM=`wc -l DRM_${cohort}_cpg_index${index}_meta_input.txt`
SVLM=`wc -l SVLM_${cohort}_cpg_index${index}_meta_input.txt`

BFline=`echo $BF | cut -f 1 -d ' '`
DRMline=`echo $DRM | cut -f 1 -d ' '`
SVLMline=`echo $SVLM | cut -f 1 -d ' '`

BFcpg=`cut -f 7 BF_${cohort}_cpg_index${index}_meta_input.txt | sort | uniq | grep -v Probe |wc -l`
DRMcpg=`cut -f 7 DRM_${cohort}_cpg_index${index}_meta_input.txt | sort | uniq | grep -v Probe | wc -l`
SVLMcpg=`cut -f 7 SVLM_${cohort}_cpg_index${index}_meta_input.txt | sort | uniq | grep -v Probe | wc -l`

echo "meta-analysis in ${cohort} chr$chr index${index}: line - $BFline; cpg - $BFcpg"
echo "meta-analysis in ${cohort} chr$chr index${index}: line - $DRMline; cpg - $DRMcpg"
echo "meta-analysis in ${cohort} chr$chr index${index}: line - $SVLMline; cpg - $SVLMcpg"
