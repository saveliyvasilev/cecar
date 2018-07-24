#!/bin/bash
#
#SBATCH --job-name=pnn_sample_array
#SBATCH --output=%A-%a.out
#SBATCH --error=%A-%a.err
#SBATCH --nodes=5
#SBATCH --tasks-per-node=16
#SBATCH --time=00:05
#SBATCH --mem-per-cpu=15
#SBATCH --array=1-10

srun python pnn.py "$SLURM_JOB_ID $SLURM_JOB_NAME $SLURM_ARRAY_JOB_ID"


