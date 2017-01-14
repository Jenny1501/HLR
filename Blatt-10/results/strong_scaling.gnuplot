unset log
unset label
set xtic auto
set ytic auto
set title "Strong Scaling"
set ylabel "Zeit [s]"
set label "Anzahl Prozesse"

plot "STRONG_SCALING_JA.dat" using 1:4 title 'Jacobi' with linespoints,\
     "STRONG_SCALING_GS.dat" using 1:4 title 'Gauss' with linespoints
