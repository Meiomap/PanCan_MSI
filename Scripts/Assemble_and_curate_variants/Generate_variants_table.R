library(data.table)
library(tidyverse)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
#

#PCAWG
PCAWG = fread("Data/Variants/PCAWG_assembled_DDR_variants2.tsv")

#Annotate Donor_ID and Sample_ID
info = fread("Data/repository_1610444206.tsv")
info = info%>%
  dplyr::select(Donor_ID = 'ICGC Donor', Sample_ID = 'Sample ID', 'Object ID', file = "File Name")%>%
  separate(file, into = c("ID", NA), sep = "\\.")%>%
  pivot_longer(values_to = "ID", cols = -c("Donor_ID", "Sample_ID"))%>%
  dplyr::select(-name)%>%
  distinct(ID, .keep_all = T)

PCAWG = left_join(PCAWG, info, by = c("Tumor_sample_barcode" = "ID"))

#HMF
HMF = fread("Data/Variants/HMF_assembled_DDR_variants2.tsv")
print("Done loading PCAWG and HMF")

#Add yeast genes
yeast = fread("Filter_vcfs/Results/Yeast_assembled_variants.tsv")

d = bind_rows(PCAWG, HMF, yeast)
d = dplyr::select(d, Donor_ID, Sample_ID, CHROM, POS, END, REF, ALT, GENE, HGVS_P, EFFECT,
                  Variant_Type, Gnomad_AF, VAF, CADD_phred, ClinVar)
#d = as.data.frame(d)



#Filter by VAF
#d = filter(d, VAF >= 0.20|is.na(VAF))

#Format LOF by CADD and clinvar
d$ClinVar[is.na(d$ClinVar)] = "Uncertain significance"

d2=d%>%
  mutate(
         infered_state = ifelse(ClinVar %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic, risk factor", "Likely pathogenic"), "Clinvar Pathogenic", "Uncertain significance"),
         infered_state = ifelse(ClinVar %in% c("Benign", "Benign/Likely benign", "Likely benign"), "Clinvar Benign", infered_state),
         infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred>=25, "Pathogenic VUS (CADD > 25)", infered_state),
         infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred >= 20, "VUS (CADD > 20)", infered_state),
         infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred >= 15, "VUS (CADD > 15)", infered_state),
         infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred >= 10, "VUS (CADD > 10)", infered_state),
         infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred <= 10, "Benign VUS (CADD < 10)", infered_state)
         )%>%
  mutate(Mutation_level = ifelse(infered_state == "Clinvar Pathogenic", 7, NA),
         Mutation_level = ifelse(infered_state == "Pathogenic VUS (CADD > 25)", 6, Mutation_level),
         Mutation_level = ifelse(infered_state == "VUS (CADD > 20)",  5, Mutation_level),
         Mutation_level = ifelse(infered_state == "VUS (CADD > 15)",  4, Mutation_level),
         Mutation_level = ifelse(infered_state == "VUS (CADD > 10)",  3, Mutation_level),
         Mutation_level = ifelse(infered_state == "Benign VUS (CADD < 10)",  2, Mutation_level),
         Mutation_level = ifelse(infered_state == "Clinvar Benign",  1, Mutation_level)
         
  )

table(d2$infered_state)
table(d2$Mutation_level)

#whitelist_pcawg
pcawg = d2[grep("SA", d2$Sample_ID),]
white = fread("Data/whitelist_PCAWG_15_2_2021.tsv")
pcawg = filter(pcawg, Sample_ID %in% white$Sample_ID)

hmf = d2[grep("SA", d2$Sample_ID, invert = T),]
hmf = hmf[grep("T$", hmf$Sample_ID),]
d5 = bind_rows(hmf, pcawg)


write.table(d5, "Results/Variants.tsv", sep ="\t", col.names = T, row.names = F, quote = F)

soma = filter(d5, Variant_Type %in% c("Somatic SNV", "Somatic insertion", "Somatic deletion"))
write.table(soma, "Data/Somatic_Variants.tsv", sep ="\t", col.names = T, row.names = F, quote = F)


###Annotate lynch
