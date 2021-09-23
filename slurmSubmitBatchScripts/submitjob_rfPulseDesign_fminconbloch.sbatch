#!/bin/bash

#SBATCH --job-name=job_rfPulseDesign_fminconbloch
#SBATCH --output=job_rfPulseDesign_fminconbloch.out
#SBATCH --error=job_rfPulseDesign_fminconbloch.err
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
matlab -nodisplay < main3_rfPulseDesign_fminconBloch_cluster.m