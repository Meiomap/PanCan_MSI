#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 60G
#SBATCH -c 1
#SBATCH --account pcawg
#SBATCH -t 24:00:00

source activate R35
Rscript /home/simong/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Code/Count_Mononucleotide_deletions.R
Rscript /home/simong/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Code/Count_Mononucleotide_insertions.R
Rscript /home/simong/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Code/Count_Dinucleotide_deletions.R
Rscript /home/simong/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Code/Count_Dinucleotide_insertions.R

#59359416