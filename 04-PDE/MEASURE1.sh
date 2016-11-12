#!/bin/bash

#SBATCH --partition=west
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1
#SBATCH --output=Messung1.out

touch Laufzeiten1.out
run="Laufzeiten1.out"
TIMEFORMAT='%E'
rm $run

for (( i = 1; i <= 12; i++ )); do
	echo "$i. Thread:" >> $run
	echo -n "$i. Thread: "
	runtime=$(time (./partdiff-openmp $i 2 512 2 2 1024 >> $run ))	
	echo $runtime
done
