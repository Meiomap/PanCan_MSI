###
## Association Between CT and MSI
##
###

library(data.table)
library(tidyverse)
library(extrafont)
#library(cutpointr)
#library(patchwork)

#Figure A
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

source("Figures/FigX3/Figure_3A.R")
gA
#ggsave(plot = gA, "Figures/Figure_2/Fig2_A.pdf", device = cairo_pdf, width = 7, height = 3)

#Figure C
source("Figures/FigX3//Figure_3C.R")
source("Figures/FigX3/Figure_3C_PanCancer.R")
gC_pan

# gC = gC + theme(
#   #axis.text.y = element_blank(),
#   axis.title.y = element_blank(),
#   axis.ticks.y = element_blank()
# )
# gC

library(patchwork)
layout = "
aaaaaaaaaa
aaaaaaaaaa
aaaaaaaaaa
aaaaaaaaaa
bbcccccccc
bbcccccccc
"

out =gA + gC_pan + gC +
  patchwork::plot_layout(design = layout) +
  plot_annotation(
    title = "Figure X3",
    subtitle = '',
    #caption = ''
    tag_levels = 'a'
  ) &
    theme(
      text = element_text( size = 7),
      title =  element_text( size = 7),
      axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 7),
      axis.title = element_text(size = 7)
    )
out

ggsave(plot = out, "Figures/FigX3/FigX3.pdf", device = cairo_pdf, width = 8, height = 7)
