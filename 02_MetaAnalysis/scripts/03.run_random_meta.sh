#!/bin/bash
#SBATCH --job-name=random_meta
#SBATCH --partition interruptible_cpu,gpu,cpu
#SBATCH --time=12:0:0
#SBATCH --mem=128G

exec &> >(tee ../process_reports/random_meta_${1}_chr${3}_index${2}.out)
# $1 is BF/SVLM/DRM $2 is index
echo "running meta-analysis of vmeQTLs - method ${1}, cpg index ${2}"

#cd ../../01_DataProcess/data/chr${3}/
cd /scratch/prj/dtr/recovered/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/Xiaopu/GoDMC_vQTL/data/chr${3}

/scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/random-metal/executables/metal random_metal_${1}_cpg_index${2}_run.txt
mv ${1}_random_meta_results_cpg_index${2}* /scratch/prj/bell/recovered/epigenetics/Analysis/subprojects/xiaopu/GoDMC/godmc2_vmeQTL/02_MetaAnalysis/data/Meta_Results/chr${3}

echo "meta-analysis of vmeQTLs on cpg index - method ${1} cpg index ${2} finished"
