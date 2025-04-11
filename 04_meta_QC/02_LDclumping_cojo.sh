#!/bin/bash/

source config

exec &> >(tee ${vmeQTL_02_logfile}_chr${1})
i=${1}

cut -f 2 ${input_data_dir}/vmeQTL_vCpG_pairs_1e-5_three_methods_chr${i}.txt | while read -r line; do
${plink} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
      --clump ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_SVLM.LDinput \
      --clump-r2 0.2 \
      --clump-kb 2000 \
      --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_SVLM.LDinput

${plink} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
      --clump ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_DRM.LDinput \
      --clump-r2 0.2 \
      --clump-kb 2000 \
      --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_DRM.LDinput

${plink} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
      --clump ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_BF.LDinput \
      --clump-r2 0.2 \
      --clump-kb 2000 \
      --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_BF.LDinput

awk 'BEGIN{FS=" "}{print $3}' ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_SVLM.clumped > ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_SVLM.snplist
awk 'BEGIN{FS=" "}{print $3}' ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_DRM.clumped > ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_DRM.snplist
awk 'BEGIN{FS=" "}{print $3}' ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_BF.clumped > ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_BF.snplist

${gcta} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
        --cojo-file ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_SVLM.cojo.ma \
        --extract ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_SVLM.snplist \
        --cojo-slct \
        --cojo-p ${cojo_p} \
        --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_SVLM.cojo

${gcta} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
        --cojo-file ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_DRM.cojo.ma \
        --extract ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_DRM.snplist \
        --cojo-slct \
        --cojo-p ${cojo_p} \
        --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_DRM.cojo

${gcta} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
        --cojo-file ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_BF.cojo.ma \
        --extract ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_BF.snplist \
        --cojo-slct \
        --cojo-p ${cojo_p} \
        --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_BF.cojo

echo "LD and COJO have been successfully run of ${line}"
done

${R_directory}/Rscript ${vmeQTL_script_dir}/combine_COJO_output.R ${processed_data_dir}/LD_COJO_input/chr${i}
echo "LD and COJO results have been successfully summarized for chr${i} vmeQTLs"
