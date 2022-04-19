# 2 Proportion of mono-nucleotide and di-nucleotide stretches per gene (~30% mono, 8%di in average)
library(tidyverse)
library(data.table)
library(ggpubr)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Supp. Figures/Del_ins_counts_S4_to_S8/")

#Get per-sample ploidy


cnv_hmf = fread("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/DDR_project_May_2021/Combined_somatic_files/HMF_CNV_COI.cnv")
cnv_hmf$total_cn = cnv_hmf$minorAllelePloidy+cnv_hmf$majorAllelePloidy
cnv_hmf$length = cnv_hmf$end-cnv_hmf$start
cnv_pcawg = fread("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/DDR_project_May_2021/Combined_somatic_files/PCAWG_CNV_Consensus.tsv")
cnv_pcawg$length = cnv_pcawg$end-cnv_pcawg$start

d = bind_rows(
  dplyr::select(cnv_hmf, ID = Sample_ID, total_cn, length),
  dplyr::select(cnv_pcawg, ID, total_cn, length),
)%>%
  mutate(
    weighted_cnv = length*total_cn,
  )%>%
  group_by(ID)%>%
  mutate(
    full_length = sum(length),
    full_cnv = sum(weighted_cnv, na.rm = T), 
    weighted_cn = full_cnv / full_length
  )%>%
  summarise(median = median(weighted_cn, na.rm = T), 
            q10 = quantile(weighted_cn, 0.1, na.rm = T), 
            q90 = quantile(weighted_cn, 0.9, na.rm = T))

  
hist(d$median)
summary(d$median)

#translate IDs
source("Code/Translate_IDs.R")
d = translate_id(d)

g1 = ggplot(d, aes(x = (median)))+
  geom_histogram()+
  scale_x_continuous(breaks = 0:8, labels = 0:8, expand = c(0.1, 0))+
  theme_classic2(base_size = 7)+
  scale_y_continuous(labels = scales::comma_format())+
  ylab("No. of tumours")+
  xlab("Median ploidy across the genome")
g1


#Add msi
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Supp. Figures/Del_ins_counts_S4_to_S8/")
msi = fread("../../Data/MSI_status_all_samples.tsv")
d = left_join(d, msi)
d$n = d$n_di + d$n_mo
d = na.omit(d)
d$MSI = d$msiStatus

cor(d$median, d$n, method = "spearman")
g2 = ggplot(d, aes(y = n, x = (median)))+
  geom_point(size = 0.5, alpha = 1, aes(col = MSI))+
  scale_x_continuous(breaks = 0:8, labels = 0:8, expand = c(0.1, 0))+
  theme_classic2(base_size = 7)+
  scale_y_continuous(labels = scales::comma_format())+
  xlab("Median ploidy across the genome")+
  ylab("Mono- and dinucleotide repeat indels")+
  scale_color_manual(values = c("brown", "darkblue"))
g2

library(patchwork)
layout = "
AA
BB
"

out =g1 + g2 +
  patchwork::plot_layout(design = layout, guides = "collect") +
  plot_annotation(
    title = "Supplementary figure S5",
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

#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Supp. Figures/Del_ins_counts_S4_to_S8/")
#ggsave(plot = out, "../SX6.pdf", device = cairo_pdf, width = 4, height = 3.5)

# library('grid')
# library(gridExtra)
# grid.newpage()
# grid.draw(cbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
# g = arrangeGrob(g1,g2, nrow=2)

