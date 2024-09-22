set title 'center value of 1D Convection-Diffusion along x direction'
set grid
set xlabel 'Pe_L'
set ylabel 'Non-dimensional phi_C'
set key right top
p 'center_temp_x_center.dat' u 1:3 w lp lw 1.5 lc rgb 'red' ti 'Exact Solution', 'center_temp_x_center.dat' u 1:2 w lp lw 1.5 lc rgb 'blue' ti 'Central difference scheme solution'
rep 'center_temp_x_upwind.dat' u 1:2 w lp lw 1.5 lc rgb 'orange' ti 'Upwind scheme solution'