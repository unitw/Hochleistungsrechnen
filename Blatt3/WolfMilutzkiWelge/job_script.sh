#!/bin/sh

#SBATCH -N 4
#SBATCH -n 16
#SBATCH --ntasks-per-node 4
#SBATCH --output timescript.out

srun timescript.sh >> timescript.out
echo fertig > job_script.out