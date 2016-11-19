#!/bin/bash

#SBATCH --partition=west
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1
#SBATCH --output=Messung2.out

touch Laufzeiten2.out
run="Laufzeiten2.out"
TIMEFORMAT='%E'

rm $run

for k in {1,2,4,8,16,32,64,128,256,512,1024}; do
	echo "$k. Interline: " >> $run
	echo -n "$k. Interline: "

	for (( l = 1; l <= 3 ; l++ )); do
		echo "$l Lauf: " 
		runtime=$(time (./partdiff-openmp 12 2 $k 2 2 1024 >> $run))
		echo $runtime
	done
done
