#!/bin/bash

#BEGIN INPUT
pdb=complex.pdb
ff=6
wat=1
box=10.0
salt=0.15
ionGrp=13
#END INPUT

#Generate first topology file (protein + crystal water/ions)
gmx_plumed pdb2gmx -ignh -f ${pdb} -p topol_prot.top -o prot.gro << STOP
${ff}
${wat}
STOP

#Define simulation box for preiodic boundary conditions (PBC)
gmx_plumed editconf -f prot.gro -box ${box} -o box.gro

#Make copy of protein topology file and delete crystal water/ions from original
cp topol_prot.top topol.top
awk '{if($1!="SOL"&&$1!="NA"&&$1!="CL") {printf("%s\n",$0);}}' topol.top >& topol_prot.top

#Solvate the protein and update new topology file (keep crystal water/ions)
gmx_plumed solvate -cp box.gro -cs -p -o solv.gro
rm \#topol.top.1\#

#Add ions to the solution (requires *.tpr input file)
gmx_plumed grompp -f em.mdp -c solv.gro -p -o tmp.tpr -maxwarn 1
gmx_plumed genion -s tmp.tpr -p -o prep.gro -neutral -conc ${salt} << STOP
${ionGrp}
STOP
rm \#topol.top.1\#
