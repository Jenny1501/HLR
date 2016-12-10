#!/bin/bash
	
#SBATCH --time=1
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --partition=west
#SBATCH --output=job.out
#SBATCH --error=job.err

. /etc/profile.d/wr-spack.sh
spack load --dependencies mpi

if [ "${SLURM_PARTITION}" != 'abu' ]
then
	export MPICH_NEMESIS_NETMOD=tcp
fi

mpiexec ./partdiff-seq 1 2 0 2 2 85
