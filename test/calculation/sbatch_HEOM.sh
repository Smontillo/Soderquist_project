#!/bin/bash
#SBATCH -p action
#SBATCH -o output.log
#SBATCH --mem=60GB
#SBATCH --time=10-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --job-name=sec
#SBATCH --open-mode=append

time /scratch/smontill/HEOM/bin/rhot ./input.json

