unset log
unset label
set xtic auto
set ytic auto
set title "Weak Scaling"
set ylabel "Zeit [s]"
set xlabel "Problemgroeße"

plot "WEAK_SCALING_JA.dat" using 1:4 title 'Jacobi' with linespoints,\
     "WEAK_SCALING_GS.dat" using 1:4 title 'Gauss' with linespoints
