#!/bin/bash

#SBATCH --time=1
#SBATCH --nodes=3
#SBATCH --ntasks=6
#SBATCH --partition=west
#SBATCH --output=job.out
#SBATCH --error=job.err

. /etc/profile.d/wr-spack.sh
spack load --dependencies mp
if [ "${SLURM_PARTITION}" != 'abu' ]
then
    export MPICH_NEMESIS_NETMOD=tcp
fi

mpiexec ./timempi2
