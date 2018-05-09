#!/bin/bash
#
#SBATCH --job-name=test1
#SBATCH --output=res1.%j.txt
#SBATCH --ntasks=2
#SBATCH --time=00:05
#SBATCH --mem-per-cpu=1

srun ./a.out
