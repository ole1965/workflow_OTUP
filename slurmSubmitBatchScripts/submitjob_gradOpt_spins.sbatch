#!/bin/bash

#SBATCH --job-name=job_gradOpt_spins
#SBATCH --output=job_gradOpt_spins.out
#SBATCH --error=job_gradOpt_spins.err
# Queue (Partition):
#SBATCH --partition=general
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=1
#SBATCH --cpus-per-task=72
#
#SBATCH --mail-type=end
#SBATCH --mail-user=ole.geldschlaeger@tuebingen.mpg.de
#
# Wall clock limit:
#SBATCH --time=23:59:00

module load matlab
matlab -nodisplay < main1_gradOpt_spins.m