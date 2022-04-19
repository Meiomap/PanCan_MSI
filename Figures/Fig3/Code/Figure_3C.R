library(data.table)
library(tidyverse)
library(extrafont)
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

d = dplyr::filter(d, 
                  Gene %in% c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "POLE", "EXO1", "POLD1"),
                  cancertype %in% c("Small intestine cancer", "Endometrial cancer", "Colon cancer", "Prostate cancer"))%>%
  dplyr::select(Sample_ID, cancertype, Gene, Mutation_level, msiStatus, state)%>%
  distinct(Sample_ID, Gene, state, .keep_all = T)%>%
  group_by(Sample_ID, Gene)%>%
  mutate(
    state = ifelse(sum(state == "Germline")>0, "Germline", "Somatic"),
    Double = ifelse(n()>1, T, F),
    LOF = T
  )%>%
  distinct(Sample_ID, Gene, .keep_all = T)%>%
  dplyr::select(Sample_ID, cancertype, Gene, state, Double, msiStatus)

msi2 = filter(msi, cancertype %in% c("Small intestine cancer", "Endometrial cancer", "Prostate cancer", "Colon cancer"))%>%
  dplyr::select(cancertype, Sample_ID, msiStatus)%>%
 # anti_join(d, by = "Sample_ID")%>%
  mutate(Gene = "None")

d2 = bind_rows(d, msi2)%>%
  filter(msiStatus == "MSI")%>%
  group_by(cancertype)%>%
  mutate(n_ct = n_distinct(Sample_ID))%>%
  group_by(cancertype, Gene, state)%>%
  mutate(n = n_distinct(Sample_ID))%>%
  ungroup()%>%
  mutate(
    prop = n/n_ct
  )%>%
  distinct(cancertype, Gene, state, n, n_ct, prop)%>%
  filter(Gene != "None")

d3 = expand_grid(cancertype = d2$cancertype, Gene = d2$Gene, state = d2$state)%>%
  mutate(n = 0)%>%
  anti_join(d2)%>%
  distinct()

d3 = left_join(d3, dplyr::select(d2, cancertype, n_ct))
d3 = bind_rows(d2,d3)
d3 = distinct(d3)
d3$prop[is.na(d3$prop)] = 0
f2c <<-d3
ll = c( "Small intestine cancer","Prostate cancer","Endometrial cancer", "Colon cancer")
d3$cancertype = factor(d3$cancertype, levels = (ll), labels = c("Small intestine cancer\n3 MSI cancers",
                                                                "Prostate cancer\n21 MSI cancers",
                                                                "Endometrial cancer\n11 MSI cancers", 
                                                                "Colon cancer\n27 MSI cancers"))

gC = ggplot(d3, aes(x = Gene, y = n, fill = state))+
  geom_histogram(stat = "identity", width = 0.6, position = "dodge", color = "black")+
  scale_color_brewer()+
  facet_grid(cols =vars(cancertype), scales = "free")+
  ylab("No. of tumours")+
  guides(fill = guide_legend(title = ""))+
  theme_bw(base_size = 9) +
  theme(axis.title = element_text(size = 7),
        axis.text.x = element_text(angle= 45, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        panel.border = element_rect(color = "black", fill = NA),
        strip.text = element_text(size  = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        axis.line = element_line(color = NA))+
  xlab("")+
  scale_fill_manual(values = c("darkslategrey", "goldenrod"))

gC
#ggsave(plot = gC, "Figures/Figure_2/Fig2C.pdf", device = "pdf", width = 7, height = 2.5)  
  


