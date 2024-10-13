set title 'center value of 1D Convection-Diffusion along x direction'
set grid
set xlabel 'log(dx)'
set ylabel 'log(Error)'
set key right top
set xrange [-4:0]

stats 'convergence_test_SOU.dat' using 1:2 nooutput
x1 = STATS_min_x     
y1 = STATS_min_y 
slope = 1.90
print "First point: (", x1, ",", y1, ")"
f(x) = slope * (x - x1) + y1

p 'convergence_test_SOU.dat' u 1:2 with points pointsize 5 lc rgb 'red' ti 'SOU Numerical'
rep f(x) w l lw 2 lc rgb 'blue' ti sprintf('Line with slope %.2f', slope)
