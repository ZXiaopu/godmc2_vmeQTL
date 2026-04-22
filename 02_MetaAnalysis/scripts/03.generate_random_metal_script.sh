#!/bin/bash

chr=$1
index=$2

cd ../../01_DataProcess/data/chr${chr}
cp ../../../02_MetaAnalysis/scripts/metal_BF_run_temple.txt random_metal_BF_cpg_index${index}_run.txt
cp ../../../02_MetaAnalysis/scripts/metal_DRM_run_temple.txt random_metal_DRM_cpg_index${index}_run.txt
cp ../../../02_MetaAnalysis/scripts/metal_SVLM_run_temple.txt random_metal_SVLM_cpg_index${index}_run.txt

for i in `ls BF*_cpg_index${index}_meta_input.txt`;
do
echo "PROCESS $i" >> random_metal_BF_cpg_index${index}_run.txt
done
echo "OUTFILE BF_random_meta_results_cpg_index${index} .txt" >> random_metal_BF_cpg_index${index}_run.txt
echo "ANALYZE RANDOM" >> random_metal_BF_cpg_index${index}_run.txt

for i in `ls DRM*_cpg_index${index}_meta_input.txt`;
do
echo "PROCESS $i" >> random_metal_DRM_cpg_index${index}_run.txt
done
echo "OUTFILE DRM_random_meta_results_cpg_index${index} .txt" >> random_metal_DRM_cpg_index${index}_run.txt
echo "ANALYZE RANDOM" >> random_metal_DRM_cpg_index${index}_run.txt

for i in `ls SVLM*_cpg_index${index}_meta_input.txt`;
do
echo "PROCESS $i" >> random_metal_SVLM_cpg_index${index}_run.txt
done
echo "OUTFILE SVLM_random_meta_results_cpg_index${index} .txt" >> random_metal_SVLM_cpg_index${index}_run.txt
echo "ANALYZE RANDOM" >> random_metal_SVLM_cpg_index${index}_run.txt

