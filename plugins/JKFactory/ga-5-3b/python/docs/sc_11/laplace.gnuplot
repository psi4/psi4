set terminal png font "./FreeSerif.ttf,16"
set output 'laplace.png'
#set terminal epslatex
#set output 'laplace.eps'
set key inside left bottom horizontal
#set title "laplace.py Strong Scaling for 10K x 10K Matrix"
set logscale xy
set xlabel "Cores"
set ylabel "Time (s)"
set pointsize 1.5
plot \
"laplace.dat" u 1:2 t "Solver Time" w linespoints lt 1 pt 3 lc rgb "#D57500", \
"laplace.dat" u 1:3 t "Wall Time"   w linespoints lt 2 pt 5 lc rgb "#707276"
