library(tidyverse)
library(data.table)
#setwd("/Volumes/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI/Filter_vcfs")

d = fread("HMF_out/Germ/barcoded/combined.ann.filtered.vcf", skip = 0, sep = "\t", header = T)
d = filter(d, !`ANN[0].EFFECT` %in% c("intron_variant","intragenic_variant"))
cat("N_row", nrow(d))
#d = fread("combined.ann.filtered.vcf", skip = 0, sep = "\t", header = T)
d = d[2:nrow(d),]
d = rename(d, HGVS_P = 'ANN[0].HGVS_P', GENE = 'ANN[0].GENE', EFFECT = "ANN[0].EFFECT", AA_POS = "ANN[0].AA_POS", CDS_POS = "ANN[0].CDS_POS")
d = separate(d, col = ID, into = c("ID","Tumor_sample_barcode", "INFO"), sep = "\\|")
print("Loaded combined.ann.filtered.vcf")

KVsep <- fixed(";")  #key-value separator
Vsep <- fixed("=")     #value separator

d2 <-  d %>%
  mutate(KVpairs = str_split(INFO, KVsep)) %>%
  unnest(KVpairs) %>%
  separate(KVpairs, into = c("key", "value"), Vsep) %>% #Split INFO
  spread(key, value) %>%
  dplyr::select(-INFO)%>%
  select_if(~sum(!is.na(.)) > 0)

print("Starting to annotate CADD and GNOMAD")

db_SNV = fread("data/CADD_SNV.tsv", sep ="\t", header = F)
db_indels = fread("data/CADD_INDELS.tsv", sep ="\t", header = F)

# db_SNV = fread("/Volumes/GenomeDK/PCAWG/CADD/DDR_genes/DDR_CADD_SNV.tsv", sep ="\t", header = T)
# db_indels = fread("/Volumes/GenomeDK/PCAWG/CADD/DDR_genes/DDR_CADD_INDELS.tsv", sep ="\t", header = T)

db = bind_rows(db_SNV, db_indels)
#db$INFO = NULL; db$FILTER = NULL
colnames(db) = c("CHROM", "POS", "REF", "ALT", "CADD", "CADD_phred")
db$CHROM = as.character(db$CHROM)
db$POS = as.integer(db$POS)

d3 = left_join(d2, db, by = c("CHROM", "POS", "REF", "ALT"))

#gnomad = fread("/Volumes/GenomeDK/PCAWG/gnomad/DDR_variants_gnomad.tsv")
#gnomad = fread("/home/simong/PCAWG/gnomad/DDR_variants_gnomad.tsv")

#gnomad$ID = NULL
#gnomad = rename(gnomad, Gnomad_AF = AF)
#d4 = left_join(d3, gnomad, by = c("CHROM", "POS", "REF", "ALT"))

print("Time to save")
write.table(d3, "data/HMF_germ_variants.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
#write.table(d4, "/Volumes/GenomeDK/2019_VCF_downloads/GERMLINE_SNV/snpSift_filtered/barcoded/DDR_variants.table", sep = "\t", col.names = T, row.names = F, quote = F)
