#!/bin/bash/

source config

exec &> >(tee ${vmeQTL_02a_logfile}_chr${1}_${2})
i=${1}
method=${2}

cat ${processed_data_dir}/LD_COJO_input/chr${i}/cpglist_${method} | while read -r line; do
${plink} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
      --clump ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.LDinput \
      --clump-r2 0.2 \
      --clump-kb 2000 \
      --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.LDinput

if [ -f ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.LDinput.clumped ]; then
    awk 'BEGIN{FS=" "}{print $2}' ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.LDinput.clumped > ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.snplist
else
    awk 'BEGIN{FS=" "}{print $2}' ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.LDinput > ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.snplist
fi

${gcta} --bfile ${processed_data_dir}/genetic_data/vmeQTL.chr${i}.merged \
        --cojo-file ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.cojo.ma \
        --extract ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.snplist \
        --cojo-slct \
        --cojo-p ${cojo_p} \
        --out ${processed_data_dir}/LD_COJO_input/chr${i}/${line}_${method}.cojo

done

echo "LD and COJO have been successfully done"
