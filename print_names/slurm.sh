#!/bin/bash
#
#SBATCH --job-name=pnn_sample
##SBATCH --partition=freeriders
#SBATCH --output=%j.out
#SBATCH --error=%j.err
##SBATCH --ntasks=17
#SBATCH --nodes=5
#SBATCH --tasks-per-node=16
#SBATCH --time=00:05
#SBATCH --mem-per-cpu=15
##SBATCH --array=1-10

#srun python pnn.py $SLURM_ARRAY_JOB_ID

srun python pnn.py "$SLURM_JOB_ID $SLURM_JOB_NAME"

