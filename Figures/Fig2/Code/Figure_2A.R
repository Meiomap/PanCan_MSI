####Overview of the data
library(data.table)
library(tidyverse)
library(patchwork)
library(cutpointr)
#library("ggplotify")
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

d = fread("Data/MSI_status_all_samples.tsv")
table(duplicated(d$Donor_ID))
white = fread("Data/210127_extra_hmfIds.tsv")

hmf = d[grep("SA", d$Sample_ID, invert = T),]
hmf = filter(hmf, Sample_ID %in% white$sampleId)
PCAWG = d[grep("SA", d$Sample_ID, invert = F),]
d = bind_rows(hmf, PCAWG)


d = d%>%
  group_by(cancertype)%>%
  mutate(n = n())%>%
  ungroup()

d1 <<- ggplot(d, aes(x = reorder(cancertype, -n)))+
  geom_histogram(stat = "count", aes(fill = Study), color = "darkslategrey", 
                 show.legend = T, width = 0.6, size = 0.2)+
  theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        legend.position = c(0.76, 0.75),
       # panel.border = element_rect(color = "black", fill = NA),
        legend.background = element_blank(),
        axis.line = element_line(color = NA))+
  xlab("")+
  ylab("No. of tumours")+
  scale_fill_manual(values = c("grey30", "grey90"))+
  guides(fill = guide_legend(title = "", direction = "horizontal"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,1e3), labels = scales::comma_format())

d1
# toPlot <-
#   d %>%
#   group_by(cancertype)%>%
#   mutate(bin = cut(S.ind, pretty(S.ind, 10))) %$%
#   table(cancertype, bin, Study) %>%
#   as.data.frame() %>%
#   mutate(plotVal = ifelse(Study == "PCAWG"
#                           , -1*Freq
#                           , Freq))
# 
# ggplot(toPlot, aes(x = bin, y = plotVal, fill = Study)) +
#   geom_col()

# layout = "
# AAA
# BBB
# BBB
# BBB
# "
# 
# cowplot::plot_grid(d2 + theme(axis.text.x = element_blank()), d1, 
#                    nrow = 2, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 3))
# 
# d2 + d1 + 
#   plot_annotation(tag_levels = 'A',
#                   title = 'Figure 1',
#                   subtitle = 'MSI annotation across PCAWG and HMF cancer samples',
#                   caption = 'janaury 13 2021'
#   ) +
#   plot_layout(design = layout, guides = "collect") & 
#   theme(legend.position = 'bottom')
