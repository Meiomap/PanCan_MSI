library(data.table)
library(tidyverse)

#setwd("/Volumes/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI//")
#Rscript Code/R_analysis/annotate_v16_CADD.R around 30gb RAM

####First take the Assemble variants and run them through CADD v.1.6 on the cluster to get the CADD values cleeeeaaaan

db = fread("Data/Variants/Old/HMF_v16_anno.tsv", select = c(1:5, 115,116))
db2 = fread("Data/Variants/Old/PCAWG_v16_anno.tsv", select = c(1:5, 115,116))
db = bind_rows(db, db2)
db = dplyr::select(db, CHROM ='#Chrom', POS = Pos, REF = Ref, ALT = Alt, rawScore=RawScore, CADD_phred=PHRED)

db$CHROM = as.character(db$CHROM)
db$POS = as.integer(db$POS)
#db$rawScore = NULL

HMF = fread("Data/Variants/HMF_assembled_DDR_variants.tsv")
HMF$CADD_phred = NULL
HMF = left_join(HMF, db)%>%
  distinct()
HMF$"Donor_ID" = str_remove(string = HMF$Sample_ID, pattern = "T$|TI$|TII$|TIII$|TIIII$|TIV$") #Fix donor and sample id HMF
HMF = filter(HMF, !is.na(CADD_phred))

#Add VAF to HMF
VAF = fread("../../DDR_project_Jan_2021/Data/Variants/VAF_MAF.tsv")
VAF = dplyr::select(VAF, CHROM, POS, REF, ALT, Sample_ID = Tumor_sample_barcode, VAF)
VAF$CHROM = as.character(VAF$CHROM)
HMF = left_join(HMF, VAF)

#Add clinvar annotation
clinvar = fread("Data/ClinVar_variant_summary.txt")
clinvar=dplyr::select(clinvar, gene = GeneSymbol, CHROM = Chromosome, POS = Start,
                      REF = ReferenceAlleleVCF, ALT = AlternateAlleleVCF,
                      ClinVar = ClinicalSignificance )
HMF = left_join(HMF, clinvar, by = c("CHROM", "POS", "REF","ALT"))
write.table(HMF, "Data/Variants/HMF_assembled_DDR_variants2.tsv", sep = "\t", col.names = T, row.names = F)

print("Done HMF")
rm(HMF)

PCAWG = fread("Data/Variants/PCAWG_assembled_DDR_variants.tsv")
PCAWG$CADD_phred = NULL
PCAWG = left_join(PCAWG, db)%>%
  distinct()
PCAWG = filter(PCAWG, !is.na(CADD_phred))
PCAWG = left_join(PCAWG, clinvar, by = c("CHROM", "POS", "REF","ALT"))
write.table(PCAWG, "Data/Variants/PCAWG_assembled_DDR_variants2.tsv", sep = "\t", col.names = T, row.names = F)
