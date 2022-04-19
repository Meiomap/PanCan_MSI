library(data.table)
library(tidyverse)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MMR_MSI_06_21//")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

#Regenerae germlien and somatic analysis output

# source("Figures/FigX4//Fig4a_Somatic.R")
# source("Figures/FigX4/Fig34_Germline.R")

#3A
source("Figures/FigX4/Fig4a.R")
gA = gA + theme(axis.text = element_text(size = 7))+
  guides(colour = guide_legend(override.aes = list(size=5)))
gA
#3b
source("Figures/FigX4/Fig4b.R")
gB = gB+ theme(axis.text.y = element_text(size = 7),
               axis.text.x = element_text(size = 7, face = "italic"))+
  scale_y_continuous(limits = c(0,0.3))
gB

layout = "
 AA
 BB
"
full_plot <<- gA + gB +
  plot_annotation(
    title = "Figure 4",
    #subtitle = '',
    #caption = ''
    tag_levels = 'a'
  ) +
  patchwork::plot_layout(design = layout) &
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    text = element_text( size = 7),
    title =  element_text( size = 7),
 #   axis.text.x = element_text(size = 7, angle = 40, hjust = 1, vjust = 1),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7)
  )

full_plot
#out
ggsave(plot = full_plot, "Figures/FigX4/FigX4.pdf",height = 5, width = 6, device = cairo_pdf)
  


