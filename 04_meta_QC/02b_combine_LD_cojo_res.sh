#!/bin/bash/

source config

exec &> >(tee ${vmeQTL_02b_logfile}_chr${1}_${method})
i=${1}

${R_directory}/Rscript ${vmeQTL_script_dir}/combine_COJO_output.R ${processed_data_dir}/LD_COJO_input/chr${i}
echo "LD and COJO results have been successfully summarized for chr${i} vmeQTLs"
