set xlabel "Anzahl paralleler Threads"
set ylabel "Programmlaufzeit in Sekunden"
set title "Laufzeit partdiff_openmp Jacobi mit 512 Interlines"

set datafile separator ","

set terminal png
set output 'messung1.png'

set grid
set xrange [0:13]
plot 'messung1.csv' using 1:2 with line

set xlabel "Anzahl Interlines"
set title "Laufzeit partdiff_openmp Jacobi mit 12 Threads"

set output 'messung2.png'
set xrange [0:1024]
plot 'messung2.csv' using 1:2 with line

set terminal win