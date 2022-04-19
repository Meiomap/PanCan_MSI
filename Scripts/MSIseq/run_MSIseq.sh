#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 70G
#SBATCH --account pcawg
#SBATCH -c 1
#SBATCH -t 24:00:00

source activate R35

Rscript Code/MSIseq/MSIseq_jan_2022.R
#59364332