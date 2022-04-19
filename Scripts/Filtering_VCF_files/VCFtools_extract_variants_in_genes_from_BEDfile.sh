#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 15G
#SBATCH -c 1
#SBATCH -t 10:00:00

source /com/extra/vcftools/LATEST/load.sh

for F in *.vcf.gz; do \
vcftools 
--gzvcf $F #CADD score vcf file
--bed /home/simong/PCAWG/simong_SV/subprojects/MMR_MSI/Data/yeast_genes.bed \ #BED-file with locations or intervals you want to extract
--recode-INFO-all \
--recode \
--out CADD_variants.out.tsv


