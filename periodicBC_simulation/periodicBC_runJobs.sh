#!/bin/bash

echo "Lowest Emin, Highest Emin, Spacing:"
read lowE highE dE

echo "Simulation Particle Density:"
read n

echo "Box Length:"
read L

echo "Output Destination:"
read pathOut

echo "Save interval:"
read dNsave



for E in $(seq ${lowE} ${dE} ${highE})
do
    file="periodicBC_Emin${E}_density${n}_boxL${L}.txt"

    echo "[DEFAULT]" > $file
    echo "" >> $file
    echo "" >> $file
    echo "Minimum Scattering Energy = ${E}" >> $file
    echo "Experimental Particle Density = 1e12" >> $file
    echo "Simulation Particle Density = ${n}" >> $file
    echo "Box Length = ${L}" >> $file
    echo "" >> $file
    echo "Number of CPUs = 40" >> $file
    echo "time steps = 10000000" >> $file
    echo "save interval = ${dNsave}" >> $file
    echo "base output path = ./${pathOut}/" >> $file
    echo "inital conditions = periodicBC_Emin${E}_density${n}_boxL${L}" >> $file

    sed -i '$d' periodicBC.sbatch
    echo "python ../runPeriodicBCscript.py ${file}" >> "periodicBC.sbatch"
    echo sbatch periodicBC.sbatch
    
done
