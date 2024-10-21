#!/bin/bash

rm *.res post* *.dat

cp convection_diffusion2.py convection_diffusion_ori2.py

for con in 100000 1.0 0.5 0.25 0.1 0.01
# for con in 0.01
do
    cp convection_diffusion_ori2.py convection_diffusion2.py
    Pe_L=$(bc -l <<< "scale=10; 1.0 / $con") # Bash cannot perform floating point arithmetic natively. We use bc.
    echo "run the simulation with Pe_L = ${Pe_L}"
    sed -i "s/fp(1000000) # control Pe_L/fp($con) # control Pe_L/g" convection_diffusion2.py
    python3 convection_diffusion2.py
done

mv convection_diffusion_ori2.py convection_diffusion2.py