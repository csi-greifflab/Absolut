#!/bin/bash
#SBATCH --account=nn9603k
#SBATCH --job-name=2FD6_AU
#SBATCH --time=165:00:0
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=16

set -o errexit # Make bash exit on any error
set -o nounset # Treat unset variables as errors

module restore system   # Restore loaded modules to the default
module load intel/2017a

mpiexec -n 64 ./AbominationMPI repertoire 2FD6_AU UniqueCDR3s.txt 48 /cluster/work/users/pprobert/2FD6_AUAll/ &> BatchOutput2FD6_AUAll.txt
