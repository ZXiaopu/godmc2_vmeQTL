#!/bin/bash

osca="/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/GoDMC_NSHD99/godmc_phase2/resources/bin/osca"
besd_file_path=$1
ls ${besd_file_path}/*besd > ${besd_file_path}/temp
awk 'BEGIN{FS=".besd"}{print $1}' ${besd_file_path}/temp > ${besd_file_path}/filetable.txt

for file in `cat ${besd_file_path}/filetable.txt`;
do
$osca \
        --beqtl-summary ${file} \
        --query 1 \
        --out ${file}_pval1_output
done


