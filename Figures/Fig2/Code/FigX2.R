#Fig 1
library(data.table)
library(tidyverse)
library(patchwork)
library(cutpointr)
library("ggplotify")
#library(extrafont)
#loadfonts(device = "pdf")
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

#Intro of the data figure
source("Figures/FigX2/Figure_1A.R")
gA = d1
gA
#ggsave(plot = gA, filename = "Figures/Figure_1/fig1A.pdf", device = "pdf", width = 6.8, height = 2.4)

#B
source("Figures/FigX2/Figure_1_B.R")
gB
#ggsave(plot = gB, filename = "Figures/Figure_1/fig1B.pdf", device = "pdf", width = 3.5, height = 2.4)

#C
source("Figures/FigX2/Fig1_C.R")
gC
#ggsave(plot = gC, filename = "Figures/Figure_1/fig1C.pdf", device = "pdf", width = 2.4, height = 2.4)

layout = "
 AAAAAAAAAA
 BBBBBCCCC#
"
full_plot <<- gA + gB + gC +
  plot_annotation(
    title = "Figure 2X",
    #subtitle = '',
    #caption = ''
    tag_levels = 'a'
  ) +
  patchwork::plot_layout(design = layout) & 
  theme(text = element_text( size = 7),
        title =  element_text( size = 7),
        legend.text = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title = element_text(size = 7)
        )

full_plot
#out
ggsave(plot = full_plot, "Figures/FigX2/FigX2.pdf",
       height = 5.5, width = 6.5, device = "pdf", useDingbats = T)


 
