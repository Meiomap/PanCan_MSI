############
# Assemble variants of HMF
# 
############

library(tidyverse)
library(data.table)

#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")

  print("Loading somatic")
  ####SNV####
  SNV = fread("/home/simong/HartwigMedical/faststorage/Simons_analysis/Extract_DDR_regions/snpSift_filtered/barcoded/DDR_variants.table")

  SNV = SNV%>%
    distinct_all(.keep_all = T)%>%
    filter(FILTER == "PASS")%>%
    mutate(
      CHROM = as.factor(CHROM),
      Study = "HMF",
      Variant_Type = ifelse(nchar(ALT) > nchar(REF), "Somatic insertion", "Somatic deletion"),
      Variant_Type = ifelse(nchar(ALT) ==  nchar(REF), "Somatic SNV", Variant_Type)
   #   CADD_phred = ifelse(EFFECT == "frameshift_variant" & is.na(CADD_phred), 30, CADD_phred)
    )%>%
    separate(col = Tumor_sample_barcode, into = c(NA, "Sample_ID"), sep = "\\_")%>%
    dplyr::select(CHROM, POS, REF, ALT, Sample_ID, GENE, CADD_phred, HGVS_P, EFFECT, Variant_Type, Study, FILTER, Gnomad_AF)
  
  
  ####Germline_SNV####
  print("Loading germline")
  G_SNV = fread("/home/simong/HartwigMedical/faststorage/Simons_analysis/Extract_DDR_regions_germline/snpSift_filtered/barcoded/DDR_variants.table")
 
  G_SNV = G_SNV%>%
    distinct_all(.keep_all = T)%>%
    filter(FILTER == "PASS")%>%
    mutate(
      CHROM = as.factor(CHROM),
      Variant_Type = ifelse(nchar(ALT) > nchar(REF), "Germline insertion", "Germline deletion"),
      Variant_Type = ifelse(nchar(ALT) ==  nchar(REF), "Germline SNV", Variant_Type),
     # CADD_phred = ifelse(EFFECT == "frameshift_variant" & is.na(CADD_phred), 30, CADD_phred),
      Study = "HMF"
    )%>%
    separate(col = Tumor_sample_barcode, into = c(NA, NA, "Sample_ID"), sep = "\\_")%>%
    dplyr::select(CHROM, POS, REF, ALT, Sample_ID, GENE, CADD_phred, HGVS_P, EFFECT, Variant_Type, Study, FILTER, Gnomad_AF)
  
  G_SNV$Sample_ID = paste(G_SNV$Sample_ID,"T", sep ="")

  d = bind_rows(SNV, G_SNV)
  print("d is done")
  table(d$Variant_Type)
  
  
  d = d%>%
    group_by(CHROM, POS ,ALT, REF )%>%
    mutate(n = n(),
           END = POS + nchar(ALT)
           )%>%
    ungroup()%>%
    filter(FILTER == "PASS")
  
 # table(duplicated(d[,c("Donor_ID","CHROM","POS","END", "REF", "ALT")]))
  
  d = distinct(d, Sample_ID, CHROM, POS ,END , REF, ALT, .keep_all = T)

  d = d%>%
    group_by(CHROM, POS, END, Variant_Type)%>%
    mutate(n = n())%>%
    ungroup()%>%
    filter(n < 50)

  
  write.table(d, "Data/Variants/HMF_assembled_DDR_variants.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
  