for uf in 1.0
do 
    for con in 1.0
    do
        for ncx in 5 10 15 20 25 30 35
        do
            cp convection_diffusion.py convection_diffusion_ori.py
            sed -i "s/ncx = 5 # ncx, use odd number/ncx = $ncx # ncx, use odd number/g" convection_diffusion.py
            sed -i "s/fp(1000000) # control Pe_L/fp($con) # control Pe_L/g" convection_diffusion.py
            Pe_L=$(bc -l <<< "scale=10; $uf / $con")
            echo "run the simulation with ncx = ${ncx}, Pe_L = ${Pe_L}"
            python3 convection_diffusion.py
            mv convection_diffusion_ori.py convection_diffusion.py
            echo "Original convection_diffusion.py has been restored." 
        done
    done
done 
