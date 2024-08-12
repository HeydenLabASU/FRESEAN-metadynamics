#!/bin/bash

for ((j=0;j<20;j++))
do
cp -r single_rw 07-reweight_$j
sed -i "7s/0/$j/" 07-reweight_$j/run.sh
sed -i "14s/0/$j/" 07-reweight_$j/run.sh
cd 07-reweight_$j
sbatch run.sh
cd ..
done
