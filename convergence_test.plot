set title 'Convergence test of different schemes (PeL = 1)'
set grid
set xlabel 'log(dx)'
set ylabel 'log(Error)'
set key left top
set xrange [-4:0]

stats 'convergence_test_SOU.dat' using 1:2 nooutput
x1 = STATS_min_x     
y1 = STATS_min_y 
slope1 = 2
print "First point: (", x1, ",", y1, ")"
f(x) = slope1 * (x - x1) + y1

stats 'convergence_test_upwind.dat' using 1:2 nooutput
x2 = STATS_min_x     
y2 = STATS_min_y 
slope2 = 1
print "First point: (", x2, ",", y2, ")"
g(x) = slope2 * (x - x2) + y2

p   'convergence_test_SOU.dat' u 1:2 with points pointsize 5 lc rgb 'red' ti 'SOU Numerical', \
    'convergence_test_upwind.dat' u 1:2 with points pointsize 5 lc rgb 'green' ti 'upwind Numerical', \
    'convergence_test_center.dat' u 1:2 with points pointsize 5 lc rgb 'purple' ti 'CD Numerical', \
    f(x) w l lw 2 lc rgb 'blue' ti sprintf('Line with slope %.2f', slope1), \
    g(x) w l lw 2 lc rgb 'yellow' ti sprintf('Line with slope %.2f', slope2)