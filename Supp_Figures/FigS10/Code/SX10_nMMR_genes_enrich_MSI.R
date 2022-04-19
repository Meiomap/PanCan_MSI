######
library(data.table)
library(tidyverse)
library(maftools)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

d_all = fread("Results/All_Pathogenic_variants.tsv")
n_distinct(d_all$Sample_ID)
genes = c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "POLE", "EXO1", "POLD1")
#genes = c("MSH6", "MSH3", "MLH1", "PMS1", "MRE11A", "RAD50", "NBN", "POLE", "ATR", "TOPBP1", "TOP3A", "PALB2", "XRCC4")
MMRd = filter(d_all, Gene %in% genes)
d_all2 = filter(d_all, Sample_ID %in% MMRd$Sample_ID, Gene %in% genes)
#d_all2 = filter(d_all, Sample_ID %in% MMRd$Sample_ID)

n_distinct(d_all2$Sample_ID)

msi = fread("Data/MSI_status_all_samples.tsv")
msi = filter(msi, Sample_ID %in% d_all2$Sample_ID)
msi$S.ind = msi$n_rep_ind

d = d_all2%>%
  group_by(Sample_ID)%>%
  mutate(
    n = n_distinct(Gene)
  )%>%
  dplyr::select(Sample_ID, n, S.ind)%>%
  distinct()

d$scaled_ind = scale(d$S.ind)
cor(d$n, d$S.ind)
co = cor(d$n, d$scaled_ind)

ggplot(mutate(d, n = ifelse(n<30, n, 30)), aes(x = as.factor(n), y = S.ind))+
  geom_boxplot(outlier.shape = NA, size = 0.2)+
  geom_jitter(width = 0.2, size = 0.3, col = "grey30")+
  scale_y_log10()+
  theme_bw(base_size = 7)+
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7)
  )+
  ylab("No. of mono- and dinucleotide repeat indels")+
  xlab("No. of MMR genes with pathogenic event")+
  ggpubr::stat_compare_means(comparisons = list(c(1,2), c(2,3), c(3,4)), 
                             size = 2.5)

ggsave(filename = "Supp. Figures/SX3.pdf", device = cairo_pdf, width = 3, height = 3)
