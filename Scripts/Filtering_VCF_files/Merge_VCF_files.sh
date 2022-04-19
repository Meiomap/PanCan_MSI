#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 40G
#SBATCH --account pcawg
#SBATCH -c 1
#SBATCH -t 60:00:00

####
# 
# Merges a folder of .vcf files into one combined .vcf file
# 
###


#sbatch /home/simong/PCAWG/simong_SV/subprojects/MMR_MSI/Filter_vcfs/Code/Merge_VCF_files.sh

mkdir merged

header="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tVALUES"
  #  echo -e "$header" > test/combined.vcf; \
  touch merged/combined.vcf;
    for F in *vcf; do \
        s=$F; \
        s=${s##*HMFregCPCT_}; \
        s=${s%.annotated*}; \
        s=${s%._post_processed*}; \
        s=${s%.consensus*}; \
        s=${s%.vcf*}; \
    #   echo  "$header" > barcoded/${s}.vcf; \
        grep -v "#" ${F} | awk -F "\t" -v var="$s" '{print $1"\t"$2"\t"$3"|"var"|"$8"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10}' >> merged/combined.vcf; \
   done

cd merged/

source /com/extra/bcftools/LATEST/load.sh

sort -k1,1V -k2,2n combined.vcf > sorted.combined.tsv

bgzip sorted.combined.tsv
tabix -p vcf sorted.combined.tsv.gz