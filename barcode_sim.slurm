#!/bin/bash
  
#SBATCH -J barcode_simulation
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4                            # one CPU core per task
#SBATCH -t 00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssrikan2@jhu.edu

# Load software
ml anaconda
conda activate r4-base

# Run R script with a command line argument
Rscript simulation_barcodes.R