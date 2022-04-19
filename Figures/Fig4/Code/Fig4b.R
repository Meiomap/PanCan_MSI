###
## Association between MSI and germline/somatic mutations
## 
###
library(data.table)
library(tidyverse)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")

#Somatic genes
gs = fread("Figures/FigX4/somatic_fig3.tsv")
genes_soma = filter(gs, q_val < 0.05)%>%mutate(g = "Somatic")

#Germline genes
gs = fread("Figures/FigX4/germline_fig3.tsv")
genes_germ = filter(gs, q_val < 0.05)%>%mutate(g = "Germline")

#Both
both = filter(genes_germ, Gene %in% genes_soma$Gene)%>% mutate(state  ="Both")
# genes_soma = filter(genes_soma, !Gene %in% both$Gene)%>%mutate(state  ="Somatic")
# genes_germ = filter(genes_germ, !Gene %in% both$Gene)%>% mutate(state  ="Germline")
genes = bind_rows(genes_soma, genes_germ)%>%
  mutate(state = ifelse(Gene %in% both$Gene, "Both", ""))

g = genes%>%
  mutate(
    Gene = Gene,
    MSI = n_MSI_HMF + n_MSI_PCAWG,
    n_mutated = n_mutated_HMF + n_mutated_PCAWG,
    MSS = n_mutated - MSI
  )%>%
  distinct(Gene, MSI, MSS, state, g)%>%
  pivot_longer(cols = c(MSI,MSS))%>%
  dplyr::rename(msiStatus = name, n = value)%>%
  mutate(
    n_tumours = ifelse(msiStatus == "MSS", 5949, 108),
    prop = n/n_tumours
  )

g$state = factor(g$state, levels = c("Both", ""))
gB1 <<- ggplot(g, aes(x = reorder(Gene, -prop), y = prop, fill = msiStatus))+
  geom_histogram(stat = "identity", position = "dodge", color = "black", width = 0.5, size = 0.2)+
  geom_text(aes(label = paste(" n=", n,sep =""), y = prop+0.001),size = 2.5,
            angle = 90, hjust = 0, vjust = 0.5 ,position = position_dodge(width = 0.5))+
  ggh4x::facet_nested(~state+g,
                      scales = "free",space = "free_x", switch = "both")+
  theme_bw(base_size = 7)+
  xlab("")+
  ylab("Proportion of tumours")+
  scale_fill_manual(values = c("grey30", "grey90"))+
  guides(fill = guide_legend(title =""))+
  scale_y_continuous(limits = c(0, 0.25))

gB1

 gB <<- gB1+theme_bw(base_size = 7) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        axis.title = element_text(size = 7),
        title = element_text(size = 7),
        plot.title = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      #  legend.position = c(0.9, 0.8),
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = NA))
 # ggtitle("Somatic and germline events in MSI-associated genes")+ylab("Tumours mutated in gene\n(Proportion)")

#gA
gB
#ggsave(plot = gC, "Figures/Fig3_new/Fig3C.pdf", device = "pdf", width = 7.5, height = 2.5)



#ggsave(plot = gA, "Figures/Figure_3/Fig3A.pdf", device = "pdf", width = 7, height = 3)

#Ratio
# d3 = d2%>%
#   dplyr::select(-prop)%>%
#   tidyr::pivot_wider(names_from = c(state), values_from = n_mutated)
# d3[is.na(d3)] = 0
# d3$"sum" = d3$Somatic + d3$Germline
# d3$"prop" = d3$Somatic/d3$sum
# 
# ggplot(d3, aes(x = reorder(Gene, -prop), y = prop, fill = msiStatus))+
#   geom_histogram(stat = "identity", position = "dodge")+
#   geom_text(data = d4, aes(label = paste("n=", value,sep =""), y = prop*0.5), position = position_dodge(width = 1))+
#   facet_grid(rows = vars(paste(cancertype, ", n=", n_patients, sep = "")), 
#             # cols = vars(msiStatus),
#              scales = "free")+
#   ylab("Proportion somatic hits\nin patients with pathogenic hits")
