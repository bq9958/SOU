#!/bin/bash

#rm *.res post* *.dat

for uf in -1.0
do
    for con in 0.1 0.1111111111111111 0.125 0.14285714285714285 0.16666666666666666 0.2 0.25 0.3333333333333333 0.5 1.0
    do
        cp convection_diffusion.py convection_diffusion_ori.py
        sed -i "s/fp(1.0) # uf/fp($uf) # uf/g" convection_diffusion.py
        Pe_L=$(bc -l <<< "scale=10; $uf / $con") # Bash cannot perform floating point arithmetic natively. We use bc.
        echo "run the simulation with Pe_L = ${Pe_L}"
        sed -i "s/fp(1000000) # control Pe_L/fp($con) # control Pe_L/g" convection_diffusion.py
        python3 convection_diffusion.py
        mv convection_diffusion_ori.py convection_diffusion.py
    done

done

for uf in 1.0
do
    for con in 1.0 0.5 0.3333333333333333 0.25 0.2 0.16666666666666666 0.14285714285714285 0.125 0.1111111111111111 0.1
    # for con in 1.0 0.5
    do
        cp convection_diffusion.py convection_diffusion_ori.py
        sed -i "s/fp(1.0) # uf/fp($uf) # uf/g" convection_diffusion.py
        Pe_L=$(bc -l <<< "scale=10; $uf / $con") # Bash cannot perform floating point arithmetic natively. We use bc.
        echo "run the simulation with Pe_L = ${Pe_L}"
        sed -i "s/fp(1000000) # control Pe_L/fp($con) # control Pe_L/g" convection_diffusion.py
        python3 convection_diffusion.py
        mv convection_diffusion_ori.py convection_diffusion.py
    done

done