#!/bin/bash

cd ../data
for i in `cat list1_3methods_significant_association_file_name.txt`;
do
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$4,$27}' $i | sort | uniq > ${i}.pair_uniq
done

for i in `cat list1_3methods_significant_association_file_name.txt`;
do
cat ${i}.pair_uniq | grep -v SNP_Probe >> list1_3methods_significant_results_Sep2026
done

for i in `cat list1_BF_significant_only_file_name.txt`;
do
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$4,$27}' $i > ${i}.pair_uniq
done

for i in `cat list1_BF_significant_only_file_name.txt`;
do
cat ${i}.pair_uniq | grep -v SNP_Probe >> list1_BF_significant_results_Sep2026
done

for i in `cat list2_DRM_SVLM_missing_CpGs_file_name.txt`;
do
awk 'BEGIN{FS="\t";OFS="\t"}{print $4}' $i | sort | uniq > ${i}.cpg_uniq
done

for i in `cat list2_DRM_SVLM_missing_CpGs_file_name.txt`;
do
cat ${i}.cpg_uniq | grep -v CpG >> list2_DRM_SVLM_missing_CpGs
done

sort list2_DRM_SVLM_missing_CpGs | uniq > list2_DRM_SVLM_missing_CpGs_unique
