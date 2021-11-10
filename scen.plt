
set terminal png size 1000, 450 font 'Verdana, 10'
set output 'Integration.png' 
set title "MPI Integration"
set xlabel 'Number of processes'
set ylabel 'Speedup(relative 4)'

set key left top

set xzeroaxis
set xtics 12
set xrange [4:64]

plot 'integration.runge.txt' using 1:2 with linespoints lw 1 lt rgb 'blue' title 'Runge inegration speedup', \
 'integration.montecarlo1.txt' using 1:2 with linespoints lw 1 lt rgb 'red' title 'Montecarlo integration speedup(N=10^7)', \
 'integration.montecarlo2.txt' using 1:2 with linespoints lw 1 lt rgb 'purple' title 'Montecarlo integration speedup(N=10^8)', \
 'linear.txt' using 1:2 with linespoints lw 1 lt rgb 'black' title 'Linear speedup'
