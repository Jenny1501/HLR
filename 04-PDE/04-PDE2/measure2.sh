#!/bin/bash

#SBATCH --partition=west
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1
#SBATCH --output=messung2.out

run_output="run2.out"
TIMEFORMAT='%E'

rm $run_output

for i in {1,2,4,8,16,32,64,128,256,512,1024}; do
	echo "$i Interlines:" >> $run_output
	echo "$i Interlines: "

	for (( j = 1; j <= 3 ; j++ )); do
		time_for_run=$(time (./partdiff-openmp 12 2 $i 2 2 1024 >> $run_output 2>&1) 2>&1)
		echo Lauf $j: $time_for_run seconds
	done
done

