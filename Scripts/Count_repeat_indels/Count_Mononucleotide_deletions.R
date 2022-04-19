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
  mutate(
    length_ref = nchar(REF),
    length_alt = nchar(ALT),
    TYPE = ifelse(length_ref>length_alt, "DEL", "INS")
  )%>%
  filter(TYPE == "DEL")
pcawg$Del_seq = sub('.', '', pcawg$REF) #The way its annotated in the vcf, we have the first base of the reference which is the same as the alt. We remove this to see the deleted seq
pcawg$del_length = nchar(pcawg$Del_seq)

#Add flanks to categorize the variants
chr <- paste0("chr",pcawg$CHROM)
position <- pcawg$POS
#alleles <- paste0("[", pcawg$HGVS_P, "]")
offset <- 10
del_length = pcawg$del_length

pcawg$m5_seq = paste(getSeq(Hsapiens,chr,position-offset,position-1), sep = "")
pcawg$m3_seq = paste(getSeq(Hsapiens,chr,position+1,position+offset), sep ="")
#pcawg$m3__del_seq = paste(getSeq(Hsapiens,chr,position+1,position+del_length), sep ="")


#Identify deletions with only one base (might be repeated many times)
pcawg = pcawg%>%
  rowwise()%>%
  mutate(
    first_base = unlist(str_split(Del_seq,""))[1], #First base in each insertion (could be any of the bases in the insertion)
    n_rep_first_char = str_count(Del_seq, first_base), #Count how many time this base reoccurs in the event
    min_string = strrep(first_base, 3) #Copy this first base too 3 times, as this is the minimal string tto search for
  )%>%
  dplyr::filter(n_rep_first_char == del_length) #Only keep cases where all bases are he samee


#Test if min-string occurs at start of M3
pcawg$m3_rep = str_detect(pcawg$m3_seq, paste0("^",pcawg$min_string))
pcawg$m5_rep = str_detect(pcawg$m5_seq, paste0(pcawg$min_string,"$"))
pcawg$first_base = NULL; pcawg$n_rep_first_char = NULL; pcawg$min_string = NULL; pcawg$FORMAT = NULL; pcawg$del_length = NULL #Cleanup

pcawg = pcawg%>%
  filter(m3_rep == T | m5_rep ==T)%>%
  group_by(Tumor_Sample_Barcode)%>%
  mutate(
    n_m3 = sum(m3_rep),
    n_m5 = sum(m5_rep)
  )%>%
  dplyr::select(Tumor_Sample_Barcode, n_m3, n_m5)%>%
  distinct(Tumor_Sample_Barcode, .keep_all = T)

write.table(pcawg, "Data/Del_ins_counts/Mononucleotide_repeat_deletions_06_01_22.tsv", sep ="\t", col.names = T, row.names = F)