#!/bin/bash
#SBATCH --output=log.o
#SBATCH --error=log.e
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH -t 72:00:00 

source /bsuhome/mpaul/.bashrc
module load slurm
conda activate ff2
python sims/simulation.py
