#!/bin/bash

echo "Lowest Temperature, Highest Tempature, Spacing:"
read lowT highT dT

echo "Experimental Particle Density:"
read partDensity

echo "Probe Tip Location (x y):"
read x y

echo "Number of CPUs:"
read nCPU

echo "Time Steps:"
read timeSteps

for T in $(seq ${lowT} ${dT} ${highT})
do
    file="./fullCircle_T${T}_Probe${x},${y}.txt"

    echo "[DEFAULT]" > $file
    echo "" >> $file
    echo "Temperature = $T" >> $file
    echo "Experimental Particle Density = $partDensity" >> $file
    echo "Number of CPUs = $nCPU" >> $file
    echo "time steps = $timeSteps" >> $file
    echo "" >> $file
    echo "Constriction Width = 0.3" >> $file
    echo "Scattering Probability = 0" >> $file
    echo "Source Drain Ratio = 1.2" >> $file
    echo "save interval = 10000" >> $file
    echo "diameter = 10" >> $file
    echo "injector height,width = 1.5,0.6" >> $file
    echo "diffusive edges? = no" >> $file
    echo "add probe tip? = yes" >> $file
    echo "probe tip location X,Y = 0,0" >> $file
    echo "" >> $file
    echo "base output path = ./data/" >> $file
    echo "inital conditions = fullCircle_T${T}_Probe${x},${y}" >> $file
    echo "" >> $file
    echo "generate reinjection probabilities? = no" >> $file
    echo "set probs method = by length" >> $file
    echo "generated probs by simulation = []" >> $file
    
done
