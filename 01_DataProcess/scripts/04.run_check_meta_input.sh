#!/bin/bash

source ../../config

chr=$1
index=$2
cohort=`cat ${cohort_list_file}`

for c in $cohort;
do
bash 04.check_meta_input.sh $c $chr $index
done
