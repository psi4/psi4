#set terminal png giant
set terminal png font "/usr/X11/lib/X11/fonts/TTF/Vera.ttf,16"
set output 'laplace.png'
set key out vert bottom center
set title "laplace.py Strong Scaling"
set logscale xy
set xlabel "Cores"
set ylabel "Time (s)"
set pointsize 2
plot \
"laplace.dat" u 1:2 t "N=10,000" w linespoints, \
"laplace.dat" u 1:4 t "N=100,000" w linespoints
#"laplace.dat" u 1:2 t "Solver Time, N=10,000" w linespoints, \
#"laplace.dat" u 1:4 t "Solver Time, N=100,000" w linespoints
#"laplace.dat" u 1:3 t "Wall Time, N=10,000" w linespoints, \
#"laplace.dat" u 1:5 t "Wall Time, N=100,000" w linespoints

