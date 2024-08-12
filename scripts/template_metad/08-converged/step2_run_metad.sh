#!/bin/bash

# For each replica
for ((i=0;i<20;i++));
do

# Copy template
cp -r single_metad 06-metadyn_$i
cd 06-metadyn_$i

# Change the line that points to what snapshot to use
echo "Run replica $i"
sed -i "7s|0|$i|" run.sh
sed -i "15s|0|$i|" run.sh
sbatch run.sh
cd ..
done
