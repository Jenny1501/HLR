#!/bin/bash

#SBATCH --time=1
#SBATCH --nodes=3
#SBATCH --ntasks=15
#SBATCH --partition=west
#SBATCH --output=job.out
#SBATCH --error=job.err

. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" != 'abu' ]
then
	export MPICH_NEMESIS_NETMOD=tcp
fi

mpiexec ./partdiff-par 1 1 0 1 2 47
