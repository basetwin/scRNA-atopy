#!/bin/bash
#SBATCH --nodes=1                         # Total number of nodes required
#SBATCH --ntasks-per-node=1               # Number of tasks (processes) per node
#SBATCH --cpus-per-task=56                # Number of CPU cores per task
#SBATCH --partition=gpu4                  # Partition name
#SBATCH --job-name=RunRScript             # Job name
#SBATCH -o SLURM.%N.%j.out                # STDOUT file name (Standard output)
#SBATCH -e SLURM.%N.%j.err                # STDERR file name (standard error)

hostname
date

module load ANACONDA/2021.05             # Load the required module (e.g., Anaconda)


export CONDA_SHLVL=0
source ~/miniconda3/etc/profile.d/conda.sh

conda activate monocle3

Rscript atopyskinanalyse.R 
