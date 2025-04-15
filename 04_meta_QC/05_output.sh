#!/bin/bash

mkdir -p output

for i in $(seq 1 22)
do
mkdir -p output/LD_COJO_output/chr$i
cp ./processed_data/LD_COJO_input/chr$i/vmeQTL* ./output/LD_COJO_output/chr$i
done

cp -r ./processed_data/adjust_mean_effect ./output/
cp -r ./processed_data/interaction_analysis ./output/

