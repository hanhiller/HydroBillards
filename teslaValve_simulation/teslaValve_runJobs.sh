#!/bin/bash

echo "Lowest Temperature, Highest Tempature, Spacing:"
read lowT highT dT

echo "Experimental Particle Density:"
read expPartDensity

echo "Simulation Particle Density:"
read simPartDensity

echo "Save interval:"
read dNsave

for T in $(seq ${lowT} ${dT} ${highT})
do
    file="teslaValve_T${T}_n${simPartDensity}_compute_easy.txt"

    echo "[DEFAULT]" > $file
    echo "" >> $file
    echo "Temperature = $T" >> $file
    echo "Experimental Particle Density = $expPartDensity" >> $file
    echo "Simulation Particle Density = $simPartDensity" >> $file
    echo "Number of CPUs = 1" >> $file
    echo "time steps = 50000000" >> $file
    echo "" >> $file
    echo "current direction? = easy" >> $file
    echo "Scattering Probability = 0" >> $file
    echo "Source Drain Ratio = 1.14" >> $file
    echo "save interval = $dNsave" >> $file
    echo "" >> $file
    echo "diffusive edges? = yes" >> $file
    echo "reinject partilces with unity energy? = yes" >> $file
    echo "add probe tip? = no" >> $file
    echo "probe tip location X,Y = 0,0" >> $file
    echo "" >> $file
    echo "base output path = ./compute_partition_SD1.14/" >> $file
    echo "inital conditions = teslaValve_T${T}_n${simPartDensity}_compute_easy" >> $file
    echo "" >> $file
    echo "generate reinjection probabilities? = no" >> $file
    echo "set probs method = by length" >> $file
    echo "generated probs by simulation = []" >> $file

    sed -i '$d' teslaValve.sbatch
    echo "python ../runTeslaValveScript.py ${file}" >> "teslaValve.sbatch"
    sbatch teslaValve.sbatch


   file="teslaValve_T${T}_n${simPartDensity}_compute_hard.txt"

    echo "[DEFAULT]" > $file
    echo "" >> $file
    echo "Temperature = $T" >> $file
    echo "Experimental Particle Density = $expPartDensity" >> $file
    echo "Simulation Particle Density = $simPartDensity" >> $file
    echo "Number of CPUs = 1" >> $file
    echo "time steps = 50000000" >> $file
    echo "" >> $file
    echo "current direction? = hard" >> $file
    echo "Scattering Probability = 0" >> $file
    echo "Source Drain Ratio = 1.14" >> $file
    echo "save interval = $dNsave" >> $file
    echo "" >> $file
    echo "diffusive edges? = yes" >> $file
    echo "reinject partilces with unity energy? = yes" >> $file
    echo "add probe tip? = no" >> $file
    echo "probe tip location X,Y = 0,0" >> $file
    echo "" >> $file
    echo "base output path = ./compute_partition_SD1.14/" >> $file
    echo "inital conditions = teslaValve_T${T}_n${simPartDensity}_compute_hard" >> $file
    echo "" >> $file
    echo "generate reinjection probabilities? = no" >> $file
    echo "set probs method = by length" >> $file
    echo "generated probs by simulation = []" >> $file
    
    sed -i '$d' teslaValve.sbatch
    echo "python ../runTeslaValveScript.py ${file}" >> "teslaValve.sbatch"
    sbatch teslaValve.sbatch
    
done
