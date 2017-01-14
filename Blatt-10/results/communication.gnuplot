unset log
unset label
set xtic auto
set ytic auto
set title "Communication"
set ylabel "Zeit [s]"
set xlabel "Anzahl Knoten"

plot "COMMUNICATION_A_JA.dat" using 2:4 title 'Jacobi' with linespoints,\
     "COMMUNICATION_A_GS.dat" using 2:4 title 'Gauss' with linespoints
