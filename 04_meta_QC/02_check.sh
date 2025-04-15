#!/bin/bash

source config

cd processed_data/LD_COJO_input
for i in $(seq 1 22);
do 
count1=`cut -f 2,3 chr$i/vmeQTL_vCpG_pair_after_LD_COJO.txt | grep -v cpg |sort|uniq |wc -l`
ref=`wc -l chr$i/cpglist* | tail -n1 | awk '{print $1}'`
if [ $count1 == $ref ]; then 
    echo "chr$i done" 
fi 
done
