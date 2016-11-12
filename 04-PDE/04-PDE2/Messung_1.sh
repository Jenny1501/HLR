#!/bin/bash

#SBATCH -N 1
#SBATCH -c 12
#SBATCH --ntasks-per-node 1
#SBATCH --output=Messung1.out

touch "Lauf_1.out"
run="Lauf_1.out"
TIMEFORMAT='%E'
rm $run

for (( i = 7; i <= 8; i++ ));
	do
		echo -e "$i. Thread: \n"  >> $run
		echo -n "$i Thread: "

		needed_time=$(time (./partdiff-seq $i 2 512 1 2 1024 >> $run 2>&1) 2>&1)

		echo $needed_time sekunden
done
		
