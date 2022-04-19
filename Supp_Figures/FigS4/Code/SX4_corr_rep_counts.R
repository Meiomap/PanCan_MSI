# 1 Correlation between MSIse and mono and di-nucleotide counts (Show that Di-nucleotide outliers are skin cancers)
library(tidyverse)
library(data.table)
library(ggpubr)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Supp. Figures/Del_ins_counts_S4_to_S8/")

msi = fread("INDELs_MSIseq.tsv")
msi = dplyr::select(msi, ID = Tumor_Sample_Barcode, msiseq = S.ind)

#load counts
mo_del = fread("Mononucleotide_repeat_deletions_06_01_22.tsv")
mo_del = mutate(mo_del, group = "mono_del")
mo_ins = fread("Mononucleotide_repeat_insertions_06_01_22.tsv")
mo_ins = mutate(mo_ins, group = "mono_ins")
di_del = fread("Dinucleotide_repeat_deletions_06_01_22.tsv")
di_del = mutate(di_del, group = "di_del")
di_ins = fread("Dinucleotide_repeat_insertions_06_01_22.tsv")
di_ins = mutate(di_ins, group = "di_ins")
d = bind_rows(mo_del, bind_rows(mo_ins, bind_rows(di_del, di_ins)))%>%
  dplyr::rename(ID = Tumor_Sample_Barcode)

#Add missing (cases with zero)
d2 = tidyr::expand_grid(ID = unique(d$ID), group = unique(d$group), n_m3 =0, n_m5 = 0)%>%
  anti_join(d, by = c("ID", "group"))%>%
  distinct()

d3 = bind_rows(d, d2)%>%
  mutate(n = n_m3+n_m5)%>%
  distinct()

d4 = d3%>%
  ungroup()%>%
  dplyr::select(-c(n_m3, n_m5))%>%
  pivot_wider(names_from = group, values_from = n)%>%
  mutate(n_di = di_del + di_ins, 
         n_mo = mono_ins + mono_del)%>%
  dplyr::select(ID, n_di, n_mo)%>%
  left_join(msi)%>%
  distinct()

#Di vs Mo
##Add ID
source("Code/Translate_IDs.R")
d5 = translate_id(d4)
#write.table(d5, "MSI_scores_all_samples_jan_22.tsv", sep = "\t", col.names = T, row.names = F)

d5 = mutate(d5, col = ifelse(primaryTumorLocation == "Skin", "Skin", "Other"))
cor(d5$n_di, d5$n_mo, method = "spearman")
g1 = ggplot(d5, aes(x = n_di, y = n_mo))+
  geom_point(size = 0.5)+
  theme_classic2(base_size=7)+
  xlab("Dinucleotide repeat indels")+
  ylab("Mononucleotide repeat indels")+
  scale_y_continuous(labels = scales::comma_format())+
  scale_x_continuous(labels = scales::comma_format())+
  geom_smooth(method = "lm", col = "grey60", size = 0.5)+
  annotate(geom = "text", x = 15e3, y = 8e3, label = "Spearman Corr. 0.89", size = 2.5)
  #geom_point(data = filter(d5, col == "Skin"),aes(color = col))
  # scale_x_log10()+
  # scale_y_log10()
g1

# #Read highest to see if it can be true — appears it can. It's a huge stomach file (230mb) full of indels!
# tt = fread("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/Originial_data_july2019/Somatic_variants/vcf/CPCT02450014R_CPCT02450014T_post_processed.vcf", skip = "#CHROM")
# tt2 = tt%>%
#   dplyr::select(REF, ALT, FILTER)%>%
#   filter(FILTER == "PASS")%>%
#   mutate(
#     n_ref = nchar(REF),
#     n_alt = nchar(ALT)
#   )
# 
# table(tt2$n_ref == tt2$n_alt)

#Per cancertype di/mo
d6 = d5%>%
  mutate(ratio = n_di/(n_mo+n_di), n = n_di+n_mo)%>%
  filter(n > 10)

g2 = ggplot(d6, aes(x = ratio))+
  geom_histogram()+
  theme_classic2(base_size=7)+
  xlab("Ratio of di-nucleotide repeat indels\n(No. of dinucleotide repeat indels / total no. of repeat indels)")+
  ylab("No. of tumours")+
  geom_rug(sides = "b")+
  scale_y_continuous(labels = scales::comma_format())
g2

# #Msiseq
ggplot(d5, aes(x = msiseq+1))+
  geom_histogram()+
  scale_x_log10()+
  geom_vline(xintercept = 700, lty = 2, col = "red")

table(d4$msiseq>700)
g3 = ggplot(d5, aes(x = n_di+n_mo, y = msiseq))+
  geom_point(size = 0.5)+
  geom_smooth(method = "lm", col = "grey60", size = 0.5)+
 # geom_vline(xintercept = 9649, lty = 2, col = "red")+
  #geom_vline(xintercept = 11338, lty = 2, col = "blue")+
  theme_classic2(base_size=7)+
  annotate(geom = "text", x = 2.1e5, y = 200, label = "Spearman Corr. 0.81", size = 2.5)+
  xlab("Mono- and dinucleotide repeat indels")+
  scale_y_continuous(labels = scales::comma_format())+
  scale_x_continuous(labels = scales::comma_format())
g3    
cor(d4$n_di+d4$n_mo, d4$msiseq, method = "spearman")
cor(log(d4$n_di+d4$n_mo+1), log(d4$msiseq+1), method = "pearson")

library(patchwork)
layout = "
AB
CC
"

out =g1 + g3 + g2 +
  patchwork::plot_layout(design = layout, guides = "collect") +
  plot_annotation(
    title = "Supplementary figure X4",
    subtitle = '',
    #caption = ''
    tag_levels = 'a'
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    text = element_text( size = 7),
    title =  element_text( size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7)
  )
out

ggsave(plot = out, "../SX4.pdf", device = cairo_pdf, width = 6, height = 4)




# d4$n = d4$n_di + d4$n_mo
# lm(data = d4, formula = "msiseq ~ n")
# 684 / 0.024
  
  

  
