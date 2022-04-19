library(tidyverse)
library(data.table)
library('BSgenome.Hsapiens.UCSC.hg19')
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")
setwd("/home/simong/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")

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

pcawg = bind_rows(pcawg, hmf)
pcawg = filter(pcawg, CHROM != "MT")

pcawg = pcawg%>%
  rowwise()%>%
  mutate(
    length_ref = nchar(REF),
    length_alt = nchar(ALT),
    TYPE = ifelse(length_ref>length_alt, "DEL", "INS"),
    Ins_seq = sub('.', '', ALT), #The way its annotated in the vcf, we have the first base of the reference which is the same as the alt. We remove this to see the deleted seq
    ins_length = nchar(Ins_seq)
  )%>%
  filter(TYPE == "INS", ins_length == 2)

#Add flanks to categorize the variants
chr <- paste0("chr",pcawg$CHROM)
position <- pcawg$POS
#alleles <- paste0("[", pcawg$HGVS_P, "]")
offset <- 10
ins_length = pcawg$ins_length

pcawg$m5_seq = paste(getSeq(Hsapiens,chr,position-offset,position-1), sep = "")
pcawg$m3_seq = paste(getSeq(Hsapiens,chr,position+1,position+offset), sep ="")
#pcawg$m3__del_seq = paste(getSeq(Hsapiens,chr,position+1,position+del_length), sep ="")


#Remoove deletions with only one base (might be repeated many times)
pcawg = pcawg%>%
  rowwise()%>%
  mutate(
    first_base = unlist(str_split(Ins_seq,""))[1], #First base in each insertion (could be any of the bases in the insertion)
    n_rep_first_char = str_count(Ins_seq, first_base), #Count how many time this base reoccurs in the event
    min_string = Ins_seq #Copy this first base too 3 times, as this is the minimal string tto search for
  )%>%
  dplyr::filter(n_rep_first_char != ins_length) #Only keep cases where bases are different


#Test if min-string occurs at start of M3
pcawg$m3_rep = str_detect(pcawg$m3_seq, paste0("^",pcawg$min_string))
pcawg$m5_rep = str_detect(pcawg$m5_seq, paste0(pcawg$min_string,"$"))
pcawg$first_base = NULL; pcawg$n_rep_first_char = NULL; pcawg$min_string = NULL; pcawg$FORMAT = NULL; pcawg$ins_length = NULL #Cleanup

pcawg2 = pcawg%>%
  filter(m3_rep == T | m5_rep ==T)%>%
  group_by(Tumor_Sample_Barcode)%>%
  mutate(
    n_m3 = sum(m3_rep),
    n_m5 = sum(m5_rep)
  )%>%
  dplyr::select(Tumor_Sample_Barcode, n_m3, n_m5)%>%
  distinct(Tumor_Sample_Barcode, .keep_all = T)

missing = filter(pcawg, !Tumor_Sample_Barcode %in% pcawg2$Tumor_Sample_Barcode)%>%
  dplyr::select(Tumor_Sample_Barcode)%>%
  mutate(n_m3 = 0, n_m5 = 0)
pcawg = bind_rows(pcawg2, missing)

# hist(pcawg$n_m3)
# plot(pcawg$n_m3, pcawg$n_m5)

write.table(pcawg, "Data/Del_ins_counts/Dinucleotide_repeat_insertions_06_01_22.tsv", sep ="\t", col.names = T, row.names = F)

# # #Quick analysis of the file
# setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")
# pcawg = fread("Data/PCAWG_Mononucleotide_repeat_deletions_21_11_21.tsv") #Update with new pcwag where i fixed the ID spllit
# n_distinct(pcawg$Tumor_Sample_Barcode)
# 
# #d = bind_rows(d, pcawg)
# plot(d$n_m3, d$n_m5)
# 
# #Update IDs on pcawg
# #Annotate Donor_ID and Sample_ID
# info = fread("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_May_2021/Data/repository_1610444206.tsv")
# info = info%>%
#   dplyr::select(Donor_ID = 'ICGC Donor', Sample_ID = 'Sample ID', 'Object ID', file = "File Name")%>%
#   separate(file, into = c("ID", NA), sep = "\\.")%>%
#   pivot_longer(values_to = "ID", cols = -c("Donor_ID", "Sample_ID"))%>%
#   dplyr::select(-name)%>%
#   distinct(ID, .keep_all = T)
# 
# pcawg = left_join(pcawg, info, by = c("Tumor_Sample_Barcode" = "ID"))
# # white = fread("Data/whitelist_PCAWG_15_2_2021.tsv")
# # d = filter(d, Sample_ID %in% white$Sample_ID)
# # PCAWG$Tumor_sample_barcode = NULL
# cat("Distinct PCAWG samples with variants:", n_distinct(pcawg$Sample_ID) )
# table(is.na(pcawg))
# 
# whitelist = fread("Data/whitelist_PCAWG_07_04_2021.tsv")
# pcawg$Sample_ID = NULL
# pcawg = left_join(pcawg, whitelist)
# 
# #Update ID on HMF
# d = fread("Data/Mononucleotide_repeat_deletions_29_10_21.tsv")
# hmf = d[grep("CPCT|DRUP", d$Tumor_Sample_Barcode), ]
# n_distinct(hmf$Tumor_Sample_Barcode)
# hmf = dplyr::rename(hmf, Sample_ID = Tumor_Sample_Barcode)
# hmf$Donor_ID = stringr::str_replace(hmf$Sample_ID, "T$|TI$|TI*$|TIV", "")
# n_distinct(hmf$Donor_ID)
# table(is.na(hmf))
# 
# pcawg$Study = "PCAWG"
# hmf$Study = "HMF"
# d = bind_rows(hmf, pcawg)
# #d$Tumor_Sample_Barcode = NULL
# table(d$Study)
# 
# #MSIseq scores
# msi = fread("Data/MSI_status_all_samples.tsv")
# d = filter(d, Sample_ID %in% msi$Sample_ID)
# missing_pcawg = filter(msi, !Sample_ID %in% d$Sample_ID)
# table(missing_pcawg$Donor_ID %in% pcawg$Donor_ID)
# 
# d = left_join(d, msi)
# table(is.na(d$msiStatus)) #Got info on 5.5 k
# 
# d = d[!is.na(d$Sample_ID),]
# 
# cor(d$n_m3, d$S.ind)
# 
# ggplot(d, aes(x = S.ind, y = n_m3))+
#   geom_point()
# 
# #MMR state
# mmr = readxl::read_excel("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/DDR_project_May_2021/Feedback_from_reviewers/Split_by_dataset/Supp_tables/Sup_tables.xlsx", sheet = "S9")
# #mmr$Sample_ID = stringr::str_replace(mmr$Sample_ID, "T$|TI$|TI*$", "")
# d2 = left_join(d, mmr, by = "Sample_ID")
# 
# table(d2$MMR_LOF)
# ggplot(d2, aes(x = S.ind, y = n_m3))+
#   geom_point(aes(color = MMR_LOF))+
#   ggpubr::theme_classic2(base_size = 9)+
#   xlab("MSIseq\n#Indels across 10 mill. sites of repetitive DNA")+
#   ylab("Deletion with 3'-mononucletide stretch\n>4 bases")
# 
# #For Santi
# d3 = dplyr::select(d2, Sample_ID, MSIseq= S.ind, n_m3_mono = n_m3, msiStatus, msiStatus_orig, BRCA_LOF, MMR_LOF, TP53_mono_loss, primaryTumorLocation)
# write.table(d3, "/Volumes/Macintosh HD/Users/au460892/Documents/SV_HMF_Barcelona/Data/MSI_09_11.tsv", sep ="\t", col.names = T, row.names = F)
# ggsave("Figures/Amru_pres_02_11_21.pdf", device = "pdf", width = 4, height = 4)
# 
# d3 = filter(d2, n_m3 < 100, S.ind > 700)
