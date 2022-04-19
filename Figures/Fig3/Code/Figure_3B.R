library(data.table)
library(tidyverse)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

#Load MSI stats
msi = fread("Data/MSI_status_all_samples.tsv")
msi$msiStatus_orig = NULL
#d = msi
#d = filter(d, Study == "PCAWG")

all = fread("Results/All_Pathogenic_variants.tsv")

d = filter(all, VAF>0.20 | is.na(VAF))
d = filter(d, Gnomad_AF < 0.005|is.na(Gnomad_AF))


d = dplyr::filter(d, Gene %in% c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "POLE", "EXO1", "POLD1"))%>%
  dplyr::select(Sample_ID, Gene, Mutation_level, msiStatus, state)%>%
  distinct(Sample_ID, Gene, state, .keep_all = T)%>%
  group_by(Sample_ID, Gene)%>%
  mutate(
    state = ifelse(sum(state == "Germline")>0, "Germline", "Somatic"),
    Double = ifelse(n()>1, T, F),
    LOF = T
  )%>%
  distinct(Sample_ID, Gene, .keep_all = T)%>%
  dplyr::select(Sample_ID, Gene, state, Double, msiStatus)

d2 = d%>%
  filter(msiStatus == "MSI")%>%
  group_by(Gene, state)%>%
  mutate(n = n_distinct(Sample_ID))%>%
  ungroup()%>%
  mutate(
    prop = n/108
  )%>%
  distinct(Gene, state, n, prop)%>%
  filter(Gene != "None")

#Add missing
d3 = expand_grid(Gene = d2$Gene, state = d2$state)%>%
  mutate(n = 0)%>%
  anti_join(d2, by = c("Gene", "state"))%>%
  distinct()

d3 = bind_rows(d2,d3)
d3 = distinct(d3)
d3$prop[is.na(d3$prop)] = 0
f2pan <<-d3

gC_pan <<- ggplot(d3, aes(x = Gene, y = n, fill = state))+
  geom_histogram(stat = "identity", width = 0.6, position = "dodge", color = "black")+
  scale_color_brewer()+
  facet_grid(cols =vars("Pan-Cancer\nn = 108 MSI cancers"))+
  ylab("No. of tumours")+
  guides(fill = guide_legend(title = ""))+
  theme_bw(base_size = 9) +
  theme(axis.title = element_text(size = 7),
        axis.text.x = element_text(angle= 45, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA),
        strip.text = element_text(size  = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.line = element_line(color = NA))+
  xlab("")+
  scale_fill_manual(values = c("darkslategrey", "goldenrod"))

gC_pan
#ggsave(plot = gB, "Figures/Figure_2/Fig2C_PanCancer.pdf", device = "pdf", width = 1.5, height = 2.5)  
  


