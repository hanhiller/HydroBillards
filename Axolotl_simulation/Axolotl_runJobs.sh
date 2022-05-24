#!/bin/bash

echo "Lowest Temperature, Highest Tempature, Spacing:"
read lowT highT dT

echo "Experimental Particle Density:"
read partDensity

echo "Save interval:"
read dNsave

for T in $(seq ${lowT} ${dT} ${highT})
do
    file="Axolotl_T${T}_DE.txt"

    echo "[DEFAULT]" > $file
    echo "" >> $file
    echo "Temperature = $T" >> $file
    echo "Experimental Particle Density = $partDensity" >> $file
    echo "Simulation Particle Density = $partDensity" >> $file
    echo "Number of CPUs = 40" >> $file
    echo "time steps = 10000000" >> $file
    echo "" >> $file
    echo "Scattering Probability = 0" >> $file
    echo "Source Drain Ratio = 1.2" >> $file
    echo "save interval = $dNsave" >> $file
    echo "" >> $file
    echo "diffusive edges? = yes" >> $file
    echo "reinject partilces with unity energy? = no" >> $file
    echo "add probe tip? = no" >> $file
    echo "probe tip location X,Y = 0,0" >> $file
    echo "" >> $file
    echo "base output path = ./Axolotl/" >> $file
    echo "inital conditions = Axolotl_T${T}_DE_easy" >> $file
    echo "" >> $file
    echo "generate reinjection probabilities? = no" >> $file
    echo "set probs method = by length" >> $file
    echo "generated probs by simulation = []" >> $file
    
    sed -i '$d' teslaValve.sbatch
    echo "python ../runAxolotlScript.py ${file}" >> "Axolotl.sbatch"
    sbatch Axolotl.sbatch
    
done
