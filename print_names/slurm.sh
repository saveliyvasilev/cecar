#!/bin/bash
#
#SBATCH --job-name=pnn_sample_noarray
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --nodes=4
#SBATCH --tasks-per-node=10
#SBATCH --time=00:05
#SBATCH --mem-per-cpu=15MB

srun python pnn.py "$SLURM_JOB_ID $SLURM_JOB_NAME"
srun pwd
