#!/bin/bash
#SBATCH --output=log.o
#SBATCH --error=log.e
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH -t 72:00:00 

source /home/gwenwhite/.bashrc
module load slurm
conda activate forcefield
python sim/simulation.py
