#Ploting data from a pulse
set term X11 1
set title "Energy(t)"
plot "energy.dat"

set term X11 2
set title "Norm"
set yrange [0.9999:1.0002]
plot "norm.dat"

set term X11 3
set title "z-dipole"
set yrange [-0.2:0.25]
plot "zdip.dat"

set term X11 4
set yrange [-0.25:0.25]
set title "accel"
plot "accel.dat"
