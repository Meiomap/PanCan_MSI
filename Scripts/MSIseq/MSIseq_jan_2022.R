library(data.table)
library(tidyverse)
library(IRanges)
library(doParallel)
library(foreach)

#PCAWG SNV: 31082785
#51069789

#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")

#Necessary for MSIseq
source("Code/MSIseq/Compute_input_variables.R")
repeats = fread("Data/MSIseq/rep_DNA.tsv")
repeats = filter(repeats, End_Position >= Start_Position)

#Load indels
pcawg = fread("../../DDR_project_May_2021/Combined_somatic_files/tmp/PCAWG_indels.tsv")
hmf = fread("../../DDR_project_May_2021/Combined_somatic_files/tmp/HMF_indels_NO_PON.tsv")
colnames(pcawg) = c("CHROM", "POS", "INFO", "REF", "ALT", "QUAL", "FILTER", "FORMAT")
colnames(hmf) = c("CHROM", "POS", "INFO", "REF", "ALT", "QUAL", "FILTER", "FORMAT")

#Filter
pcawg = pcawg[grep("LOWSUPPORT|NORMALPANEL",pcawg$FILTER, invert = T),]
hmf = filter(hmf, FILTER == "PASS")

#Add Sample ID
pcawg = separate(pcawg, col = INFO, into = c(NA, "Tumor_Sample_Barcode", NA), sep = "\\||\\.INDELS|\\.consensus", remove = F)
hmf = separate(hmf, col = INFO, into = c(NA, "Tumor_Sample_Barcode", NA), sep = "\\||\\_p", remove = F)
hmf = separate(hmf, col = Tumor_Sample_Barcode, into = c(NA, "Tumor_Sample_Barcode"), sep = "_", remove = T)

all = bind_rows(pcawg, hmf)

vcf = all%>%
  mutate(
    Variant_Type = ifelse(nchar(REF) > nchar(ALT), "DEL", "INS"),
    Variant_Type = ifelse(nchar(REF) == nchar(ALT), "SNP", Variant_Type),
    width = abs(nchar(ALT) - nchar(REF)) ,
    Start_Position = POS,
    End_Position = POS,
  )%>%
  dplyr::select(Chrom = "CHROM", Start_Position, End_Position, Variant_Type, Tumor_Sample_Barcode, width)
#
tmp = Compute.input.variables(
  data = vcf,
  repeats = repeats, 
  uniform.seq.len = 1
)

out = tmp
out$"Tumor_Sample_Barcode" = row.names(out)


write.table(out, "Data/Del_ins_counts/INDELs_MSIseq.tsv", sep ="\t", col.names = T, row.names = F)
