library(data.table)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(cutpointr)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")

#Add MSI at PCAWG — from MSIseq/decide_cutoff_and_msiStatus.R
msi = fread("Data/MSI_status_all_samples.tsv")
msi = mutate(msi, msiStatus_orig = ifelse(msiStatus_orig == "MSI", msiStatus_orig, NA))

dev.off()
gB <<- ggplot(msi, aes(x = n_rep_ind))+
  geom_histogram(color = "black",  size = 0.2, bins = 1e2,fill = "grey90")+
  theme_classic2()+
  theme(
    legend.position = c(0.7, 0.9),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = NA)
  )+
 # facet_grid(rows = vars(Study))+
  xlab("No. of mono- and dinucleotide repeat indels")+
  ylab("No. of tumours")+
  geom_rug(aes(color = msiStatus_orig))+
  scale_color_manual(na.value = NA, values = c("brown"))+
  geom_vline(aes(xintercept = cutoff), lty = 2, col = "black", size = 0.2)+
  geom_text(data = distinct(msi, Study, .keep_all = T), aes(x = cutoff+1000, y = 100, label = paste("Cutoff = ", scales::comma( round(cutoff,1) ), sep ="")),
            angle =0, hjust = 0, size = 2.5)+
  guides(color = guide_legend(title ="| MSI status from other studies*", label = F, direction = "horizontal"))+
  facet_grid(rows = vars(Study))+
  scale_x_log10(labels = scales::comma_format()) +
  theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        #strip.text =  element_text(family = "Arial", size = 7),
        axis.title = element_text(size = 7),
        #panel.border = element_rect(color = "black", fill = NA),
        panel.spacing = unit(1, "lines"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(color = NA))+
  geom_text(aes(label = Study), data = distinct(msi, Study), x = log(1, base = 10), y = 180, hjust = 0, size = 2.5)+
  scale_y_continuous(expand = c(0,0))
  
gB
#ggsave(plot = gB, filename = "Figures/Figure_1/fig1B.pdf", device = "pdf", width = 4, height = 2.4)

  
