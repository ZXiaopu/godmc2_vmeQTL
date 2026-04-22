#!/bin/bash

source ../../config
cohort=`cat ${cohort_list_file}`
for c in $cohort;
do
    sbatch 01a.run_process_module_07.sh $c
done
