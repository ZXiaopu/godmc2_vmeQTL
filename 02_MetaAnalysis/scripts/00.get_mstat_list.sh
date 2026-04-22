#!/bin/bash

exec &> >(tee ../process_reports/getmstat_list_cohort${1}_${2}_cpgindex${3}.out)

source ../../config

Cohort=`cat ${cohort_results_path}/Cohort_Information_short.csv | head -n ${1} | tail -n1 | cut -f 1 -d ','`

vQTL_res=${cohort_results_path}
for i in $(seq 1 22);
do
cat ${vQTL_res}/${Cohort}/chr${i}/vQTL_BF*pval5e-8_output | grep -v SNP >> ../data/BF_${Cohort}_all_pval5e-8_output
cat ${vQTL_res}/${Cohort}/chr${i}/vQTL_drm*pval5e-8_output | grep -v SNP >> ../data/DRM_${Cohort}_all_pval5e-8_output
cat ${vQTL_res}/${Cohort}/chr${i}/vQTL_svlm*pval5e-8_output | grep -v SNP >> ../data/SVLM_${Cohort}_all_pval5e-8_output
done
