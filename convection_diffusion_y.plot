set title '1D Convection-Diffusion along y direction'
set grid
set xlabel 'y-coord'
set ylabel 'Scalar'
set key left top
p 'temp_y_0.0.dat' u 1:2 w l lc rgb 'red' ti 'Pe = 0', 'temp_y_1.0.dat' u 1:2 w l lc rgb 'blue' ti 'Pe = 1', \
'temp_y_2.0.dat' u 1:2 w l lc rgb 'green' ti 'Pe = 2', 'temp_y_4.0.dat' u 1:2 w l lc rgb 'orange' ti 'Pe = 4', \
'temp_y_10.0.dat' u 1:2 w l lc rgb 'purple' ti 'Pe = 10', 'temp_y_100.0.dat' u 1:2 w l lc rgb 'brown' ti 'Pe = 100'
rep 'analytical_temp_y_0.0.dat' u 1:2 w p lc rgb 'red' ti 'ref, Pe = 0', \
'analytical_temp_y_1.0.dat' u 1:2 w p lc rgb 'blue' ti 'ref, Pe = 1', \
'analytical_temp_y_2.0.dat' u 1:2 w p lc rgb 'green' ti 'ref, Pe = 2', \
'analytical_temp_y_4.0.dat' u 1:2 w p lc rgb 'orange' ti 'ref, Pe = 4', \
'analytical_temp_y_10.0.dat' u 1:2 w p lc rgb 'purple' ti 'ref, Pe = 10', \
'analytical_temp_y_100.0.dat' u 1:2 w p lc rgb 'brown' ti 'ref, Pe = 100'