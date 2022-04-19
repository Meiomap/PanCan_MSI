library(tidyverse)
library(data.table)

#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")

  ####SNV####
  SNV = fread("/home/simong/2019_VCF_downloads/SNV/snpSift_filtered/barcoded/DDR_variants.table")
  SNV = SNV%>%
    mutate(FILTER = replace_na(FILTER,"PASS"),
           CHROM = as.factor(CHROM),
           Variant_Type = "Somatic SNV",
           Study = "PCAWG"
           )%>%
    filter(FILTER == "PASS", NumCallers > 1)%>%
    dplyr::select(CHROM, POS, REF, ALT, Tumor_sample_barcode, GENE, CADD_phred, HGVS_P, EFFECT, Variant_Type, Study, FILTER, Gnomad_AF, VAF)%>%
    separate(Tumor_sample_barcode, sep = "\\.", into = c("Tumor_sample_barcode", NA), remove = T)
  
  # caf = read.table("/home/simong/PCAWG/simong_SV/Midtterm/Data/caf_values.tsv", sep = "\t", header= T) #We have CAF values in PCAWG
  # SNV = left_join(SNV, caf, by = c("CHROM", "POS", "ALT"))
  # SNV = filter(SNV, is.na(CAF)| CAF > 0.25)
  # 
  ####INDELS####
  INDELS =  fread("/home/simong/2019_VCF_downloads/INDELS/snpSift_filtered/barcoded/DDR_variants.table")
  
  INDELS = INDELS %>%
    mutate(FILTER = replace_na(FILTER,"PASS"),
          as.factor(CHROM),
          Variant_Type = ifelse(nchar(ALT) > nchar(REF), "Somatic insertion", "Somatic deletion"),
          Study = "PCAWG",
    ) %>%
    filter(FILTER == "PASS", NumCallers > 1)%>%
    dplyr::select(CHROM, POS, REF, ALT, Tumor_sample_barcode, GENE, CADD_phred, HGVS_P, EFFECT, Variant_Type, Study, FILTER, Gnomad_AF, VAF)%>%
    separate(Tumor_sample_barcode, sep = "\\.", into = c("Tumor_sample_barcode", NA), remove = T)
  
  # INDELS = left_join(INDELS, caf, by = c("CHROM", "POS", "ALT"))
  # INDELS = filter(INDELS, is.na(CAF)| CAF > 0.25)
  
  
  ####Germline_SNV####
  G_SNV =  fread("/home/simong/2019_VCF_downloads/GERMLINE_SNV_old/snpSift_filtered/barcoded/DDR_variants.table")

  G_SNV = G_SNV%>%
    mutate(FILTER = replace_na(FILTER,"PASS"),
           CHROM = as.factor(CHROM),
           Variant_Type = "Germline SNV",
           Study = "PCAWG"
    )%>%
    filter(FILTER == "PASS")%>%
    dplyr::select(CHROM, POS, REF, ALT, Tumor_sample_barcode, GENE, CADD_phred, HGVS_P, EFFECT, Variant_Type, Study, FILTER, Gnomad_AF)%>%
    separate(Tumor_sample_barcode, sep = "\\.", into = c("Tumor_sample_barcode", NA), remove = T)
  
  ####Germline_INDEL####
  G_INDELS = fread("/home/simong/2019_VCF_downloads/GERMLINE_INDELS_old/snpSift_filtered/barcoded/DDR_variants.table")
  
  G_INDELS = G_INDELS %>%
    mutate(as.factor(CHROM),
           Variant_Type = ifelse(nchar(ALT) > nchar(REF), "Germline insertion", "Germline deletion"),
           Study = "PCAWG",
          # CADD_phred = ifelse(EFFECT == "frameshift_variant" & is.na(CADD_phred), 30, CADD_phred)
    ) %>%
    filter(FILTER == "PASS")%>%
    dplyr::select(CHROM, POS, REF, ALT, Tumor_sample_barcode, GENE, CADD_phred, HGVS_P, EFFECT, Variant_Type, Study, FILTER, Gnomad_AF)%>%
    separate(Tumor_sample_barcode, sep = "\\.", into = c("Tumor_sample_barcode", NA), remove = T)
 
  d = bind_rows(SNV, INDELS, G_INDELS, G_SNV)
  
  
  d = d%>%
    group_by(CHROM, POS ,ALT, REF )%>%
    mutate(n = n(),
           END = POS + nchar(ALT)
           )%>%
    ungroup()%>%
    filter(FILTER == "PASS")

  table(duplicated(d[,c("Tumor_sample_barcode","CHROM","POS","END", "REF", "ALT")]))
  
  d = distinct(d, Tumor_sample_barcode, CHROM, POS ,END , REF, ALT, .keep_all = T)
  #d = filter(d, is.na(Gnomad_AF) | Gnomad_AF < 0.005)
  
  d = d%>%
    group_by(CHROM, POS, END, Variant_Type)%>%
    mutate(n = n())%>%
    ungroup()%>%
    filter(n < 50)
  
  write.table(d, "Data/Variants/PCAWG_assembled_DDR_variants.tsv", sep = "\t", col.names = T, row.names = F)
  
  
  # ggplot(distinct(d, CHROM, POS,END, .keep_all = T), aes(x = n, fill = Variant_Type))+
  #   geom_histogram(stat = "count")+
  #  facet_wrap(~Variant_Type, scales = "free")
  
  #Add project_code
  # primaryTumorLocation = fread("Data/PrimaryTumorLocations_allSamples.tsv", select = c("Sample_ID","primaryTumorLocation"))
  # d = left_join(d, primaryTumorLocation)
  # 
  #Filter if present in whitelist from sigToolsLib article (cause 6 bone samples are missing)
  # features  = fread("Data/COMBINED_featureCounts.tsv")
  # d = filter(d, Sample_ID %in% features$Sample_ID)
  
  #Filter away secondary samples
 # d = filter(d, !is.na(d$primaryTumorLocation))
  
  # d = mutate(d,
  #            CADD_phred = ifelse(is.na(CADD_phred), 0, CADD_phred), 
  #            Gnomad_AF = ifelse(is.na(Gnomad_AF), 0, Gnomad_AF),
  #            Impact = cut(CADD_phred, breaks = c(-10, 10, 25, 1e5), labels = c("LOW", "MODERATE", "HIGH"), include.lowest = T)
  # )
#  t = as.data.frame(table(d$EFFECT[d$CADD_phred == 0]))
  
  
  #Filter if more than 15 high impact variants across the DDR genes
  # d2 = d%>%
  #   filter(Impact == "HIGH")%>%
  #   group_by(Sample_ID)%>%
  #   filter(n() > 25)%>%
  #   distinct(Sample_ID)
  # 
#  d = filter(d, !Sample_ID %in% d2$Sample_ID) 
  
#  write.table(d, "Data/PCAWG_assembled_DDR_variants.tsv", sep = "\t", col.names = T, row.names = F)
