#!/bin/bash/

source config

exec &> >(tee ${vmeQTL_02_logfile}_chr${1}_${method})
i=${1}
method=${2}

cat ${processed_data_dir}/LD_COJO_input/chr${i}/cpglist_${method} | while read -r line; do
${plink} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
      --clump ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.LDinput \
      --clump-r2 0.2 \
      --clump-kb 2000 \
      --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.LDinput

awk 'BEGIN{FS=" "}{print $3}' ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.clumped > ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.snplist

${gcta} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
        --cojo-file ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.cojo.ma \
        --extract ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.snplist \
        --cojo-slct \
        --cojo-p ${cojo_p} \
        --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.cojo

done

echo "LD and COJO have been successfully done"

${R_directory}/Rscript ${vmeQTL_script_dir}/combine_COJO_output.R ${processed_data_dir}/LD_COJO_input/chr${i}
echo "LD and COJO results have been successfully summarized for chr${i} vmeQTLs"
