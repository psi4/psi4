set terminal png font "./FreeSerif.ttf,16"
set output 'python.png'
#set terminal epslatex
#set output 'laplace.eps'
set key inside left top vertical
#set title "'import numpy' Strong Scaling"
set logscale xy
set xlabel "Cores"
set ylabel "Time (s)"
set pointsize 1.5
plot \
"python.dat" u 1:2 t "shared python" w linespoints lt 1 pt 3 lc rgb "#D57500", \
"python.dat" u 1:3 t "bcastf"        w linespoints lt 2 pt 5 lc rgb "#707276", \
"python.dat" u 1:4 t "node python"   w linespoints lt 3 pt 7 lc rgb "#242424", \
"python.dat" u 1:5 t "bcastf + node python" w linespoints lt 4 pt 9 lc rgb "#003698"
