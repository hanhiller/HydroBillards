#!/bin/bash

echo "Lowest Temperature, Highest Tempature, Spacing:"
read lowT highT dT

echo "Experimental Particle Density:"
read expPartDensity

echo "Simulation Particle Density:"
read simPartDensity

echo "PC width:"
read constrWidth

echo "Save interval:"
read dNsave

for T in $(seq ${lowT} ${dT} ${highT})
do
    file="constriction_SD1.1_T${T}_n${simPartDensity}_noScatter_compute_13hrs.txt"

    echo "[DEFAULT]" > $file
    echo "" >> $file
    echo "Constriction Width and Thickness = ${constrWidth},.05" >> $file
    echo "" >> $file
    echo "Temperature = ${T}" >> $file
    echo "Experimental Particle Density = ${expPartDensity}" >> $file
    echo "Simulation Particle Density = ${simPartDensity}" >> $file
    echo "Scattering Probability = 0" >> $file
    echo "Number of CPUs = 1" >> $file
    echo "time steps = 10000000" >> $file
    echo "" >> $file
    echo "Source Drain Ratio = 1.1" >> $file
    echo "save interval = ${dNsave}" >> $file
    echo "" >> $file
    echo "diffusive edges? = yes" >> $file
    echo "add probe tip? = no" >> $file
    echo "probe tip location X,Y = 0,0" >> $file
    echo "" >> $file
    echo "base output path = ./Tsweep_n100_scatterNoScatter_compute_interupts/" >> $file
    echo "inital conditions = constriction_SD1.1_T${T}_n${simPartDensity}_noScatter_compute_6.5hrs" >> $file
    echo "" >> $file
    echo "generate reinjection probabilities? = no" >> $file
    echo "set probs method = by length" >> $file
    echo "generated probs by simulation = []" >> $file
    echo "reinject partilces with unity energy? = yes" >> $file
    
    sed -i '$d' constriction.sbatch
    echo "python ../runConstrictionScript.py ${file}" >> "constriction.sbatch"
    sbatch constriction.sbatch
    
done
