library(tidyverse)
library(data.table)
library(GenomicRanges)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")
#setwd("/home/simong/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")

samples = fread("Data/MSI_status_all_samples.tsv")
trans = fread("Data/repository_1610444206.tsv")
pcawg = fread("../../../../HartwigMedical/faststorage/Originial_data_july2019/Somatic_variants/vcf/CPCT02120042R_CPCT02120042T_post_processed.vcf",
              skip = "#CHROM")
pcawg = fread("../../../../2019_VCF_downloads/INDELS/9310c247-1343-5559-a5c9-2e0e6d285644.INDELS_consensus.vcf",
              skip = "#CHROM")


pcawg = fread("../../DDR_project_May_2021/Combined_somatic_files/tmp/PCAWG_indels.tsv", nrow = 1e6)
hmf = fread("../../DDR_project_May_2021/Combined_somatic_files/tmp/HMF_indels_NO_PON.tsv", nro = 1e4)
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
all = filter(all, CHROM != "MT")

all = all%>%
  mutate(
    length_ref = nchar(REF),
    length_alt = nchar(ALT),
    TYPE = ifelse(length_ref>length_alt, "DEL", "INS"),
    width = ifelse(TYPE == "INS", length_alt-1, length_ref-1)
  )

#Bethesda regions
bet = readxl::read_excel("Data/Del_ins_counts/Bethesda_panel_coordinates.xlsx")
gr1 = GRanges(seqnames=bet$Chrom,ranges=IRanges(start=bet$Start,width=bet$End-bet$Start))

ids = unique(all$Tumor_Sample_Barcode)
for(id in ids){
  tmp = filter(all, Tumor_Sample_Barcode == id)
  gr2 = GRanges(seqnames=tmp$CHROM,ranges=IRanges(tmp$POS,width=tmp$width))
  suppressWarnings(expr = {
    gr3 = intersect(gr1,gr2)
  })
  gr3 = as.data.frame(gr3)
  tmp_out = data.frame()%>%
    mutate(Tumor_Sample_Barcode = id, n_bethesda = nrow(gr3))
  
  if(id == ids[1]){
    out = tmp_out
  }else{
    out = bind_rows(out, tmp_out)
  }
}

write.table(out, "Data/Del_ins_counts/Bethesda_events.tsv", sep ="\t", col.names = T, row.names = F)


