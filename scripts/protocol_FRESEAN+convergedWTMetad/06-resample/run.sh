#!/bin/bash
#SBATCH -p general
#SBATCH -G a100:1
#SBATCH -N 1
#SBATCH -c 12
#SBATCH -t 2-00:00
#SBATCH -J RESAMPLE

next=1
nextDir=../07-metadyn

# Working directory for where the topology and starting structure are
workingDir=..
inpTOP=$workingDir/00-prep/topol.top
inpGRO=$workingDir/01-em+equi/equi/equi.gro
gmx=gmx_plumed

# Create a new directory to run 100 ns unbiased MD
mkdir run-NPT
cd run-NPT
$gmx grompp -f ../sample-states.mdp -c ${inpGRO} -p ${inpTOP} -o sample-states.tpr >& grompp.out
$gmx mdrun -v -deffnm sample-states -cpi >& mdrun.out
cd ..
 
# Pull configuration every 5 ns
mkdir snapshots
$gmx trjconv -s run-NPT/sample-states.tpr -f run-NPT/sample-states.trr -pbc mol -sep -dt 5000.0 -ndec 8 -o snapshots/state_.gro << STOP >& trjconv.out
0                                                                               
STOP 

#Start next part of the project if next=1
if [ ${next} -eq 1 ]; then
  if [ -f snapshots/state_19.gro ]; then
    if [ -d ${nextDir} ]; then
      curDir=`pwd`
      cd ${nextDir}
      if [ -f run.sh ]; then
        ./run.sh
      fi
      cd ${curDir}
    fi
  fi
fi

