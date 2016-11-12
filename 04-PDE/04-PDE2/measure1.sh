#!/bin/bash

#SBATCH --partition=west
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1
#SBATCH --output=measure1.out

run_output="run1.out"
TIMEFORMAT='%E'

rm $run_output

for (( i = 1; i <= 12; i++ )); do
	echo "$i Thread(s):" > $run_output
	echo -n "$i Thread(s): "

	time_for_run=$(time (./partdiff-openmp $i 2 512 2 2 1024 >> $run_output 2>&1) 2>&1)
	echo $time_for_run seconds
done

