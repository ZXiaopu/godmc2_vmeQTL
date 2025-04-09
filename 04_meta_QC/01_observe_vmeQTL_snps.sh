#!/bin/bash

source config
exec &> >(tee ${vmeQTL_01_logfile}_chr${1})
i=${1}

mkdir -p ${processed_data_dir}
mkdir -p ${processed_data_dir}/genetic_data
mkdir -p ${processed_data_dir}/LD_COJO_input/chr${i}

# observe vmeQTLs in each chromosome
echo "Generating vmeQTL list of chr${i}"
cut -f 1 | grep -v SNP ${input_data_dir}/vmeQTL_vCpG_pairs_1e-5_three_methods_chr$i.txt > ${processed_data_dir}/genetic_data/snplist.chr$i

# combine vmeQTL SNPs from genetic chunks
tabfile_index=$(awk -v chr=$i 'BEGIN{FS=" "}{if ($2==chr) print $1}' ${section_07_dir}/tabfile.info1)

for index in ${tabfile_index};
do
    ${plink} --bfile ${tabfile}/data.tab.${index} \
             --extract ${processed_data_dir}/genetic_data/snplist.chr$i \
             --make-bed \
             --out ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.tab${index}
    echo "vmeQTL.chr${i}.tab${index}" >> ${processed_data_dir}/genetic_data/merge_list_chr${i}
done

sort ${processed_data_dir}/genetic_data/merge_list_chr${i} | uniq > ${processed_data_dir}/genetic_data/merge_list_chr${i}.txt
rm ${processed_data_dir}/genetic_data/merge_list_chr${i}

if [ $(wc -l ${processed_data_dir}/genetic_data/merge_list_chr${i}.txt | cut -d ' ' -f 1) -gt 1 ];
then
    cd ${processed_data_dir}/genetic_data
    
    f=$(head -n1 ${processed_data_dir}/genetic_data/merge_list_chr${i}.txt)
    tail -n+2 ${processed_data_dir}/genetic_data/merge_list_chr${i}.txt > ${processed_data_dir}/genetic_data/combine_file_chr${i}.txt
    
    ${plink} --bfile ${processed_data_dir}/genetic_data/$f \
             --merge-list ${processed_data_dir}/genetic_data/combine_file_chr${i}.txt \
             --make-bed \
             --out ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged
    cd ${vmeQTL_QC_dir}
fi
echo "Merged vmeQTL genotypes from chromosome ${i}"

# generate LD clumping input file for each vCpG
${R_directory}/Rscript \
       ${vmeQTL_script_dir}/separate_vCpGs_LD_COJO_input.R \
       $i \
       ${input_data_dir}/vmeQTL_vCpG_pairs_1e-5_three_methods_chr${i}.txt \
       ${cohort_vQTL_dir}/chr${i}/vQTL_cis_BF_results_summary_${cohort_name}_chr${i}.txt \
       ${cohort_vQTL_dir}/chr${i}/vQTL_cis_drm_results_summary_${cohort_name}_chr${i}.txt \
       ${cohort_vQTL_dir}/chr${i}/vQTL_cis_svlm_results_summary_${cohort_name}_chr${i}.txt \
       ${processed_data_dir}/LD_COJO_input
                        
echo "Successfully generate genotype data and vCpG data for LD clumping and COJO"
