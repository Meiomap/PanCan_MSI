library(data.table)
library(tidyverse)
library(scales)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

#Load MSI stats
msi = fread("Data/MSI_status_all_samples.tsv")
table(duplicated(msi$Sample_ID))
table(duplicated(msi$Donor_ID))

d = msi
d$S.ind = d$n_rep_ind
table(duplicated(d$Donor_ID))
#d = filter(d, Study == "PCAWG")

d2 = d%>%
  distinct(Sample_ID, .keep_all = T)%>%
  dplyr::select(Sample_ID, Donor_ID, cancertype, S.ind, msiStatus)%>%
  mutate(
 #   msiStatus = ifelse(S.ind >= cutoff, "HIGH", "LOW"),
    n_HIGH_global = sum(msiStatus == "MSI"),
    n_LOW_global = sum(msiStatus == "MSS")
  )%>%
  group_by(cancertype)%>%
  mutate(
    n_ct = n_distinct(Sample_ID)
  )%>%
  ungroup()%>%
  distinct(Sample_ID, .keep_all = T)%>%
  group_by(cancertype)%>%
  mutate(
    n = sum(msiStatus == "MSI")
  )%>%
  ungroup()%>%
  distinct(cancertype, n_ct, n, n_HIGH_global, n_LOW_global)%>%
  mutate(prop = n/n_ct,
         n_global = (n_HIGH_global+ n_LOW_global),
         global_prop = 108/6047)%>%
#         global_prop = ifelse(Study == "PCAWG", 94/2681, 90/3517))%>%
  distinct(cancertype, .keep_all = T)%>%
  #filter(msiStatus == "MSI")%>%
  rowwise()%>%
  mutate(
    p_val = binom.test(x = n, n = n_ct, p = global_prop, alternative = "greater")$p.value,
  )

# source("Code/Helper_functions/sumlog.R")
# d3 = dplyr::select(d2, cancertype, Study, n, n_ct, p_val)%>%
#   pivot_wider(names_from= Study, values_from = c(p_val, n, n_ct))%>%
#   rowwise()%>%
#   mutate(
#     p_val = ifelse(!is.na(sum(p_val_HMF,p_val_PCAWG)) ,
#                    sumlog(c(p_val_HMF, p_val_PCAWG))$p, 
#                    max(p_val_HMF, p_val_PCAWG, na.rm = T))
#   )

d3 = filter(d2, n >1)
d3$q_val = p.adjust(d3$p_val, method = "fdr")

d3 = mutate(d3, 
            sign = ifelse(q_val < 0.05, "q \U2264 0.05", "NS"),
            prop = n/n_ct
)

dout <<-d3

gA = ggplot(d3, aes(x = reorder(paste(cancertype," (MSI: ", n, " / ", n_ct, ")", sep =""), -prop), y = prop))+
  geom_histogram(stat = "identity",aes(fill = sign), color = "black", width = 0.6, show.legend = T, size = 0.2)+
  scale_fill_manual(values = c("brown", "grey60"),
                    breaks = c("q \U2264 0.05", "NS"))+
  geom_text(data = filter(d3, q_val < 0.01), aes(label = scales::scientific(q_val, digits = 2,prefix = " q \U2264 ")),
            vjust = 0.5, hjust = 0, angle = 90, nudge_y = 0, size = 2.5)+
  geom_text(data = filter(d3, q_val > 0.01), aes(label = paste(" q \U2264 ", round(q_val, 2), sep = "")),
            vjust = 0.5, hjust = 0, angle = 90, nudge_y = 0, size = 2.5)+
  theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7),
        axis.title = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.background = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.direction = "horizontal",
        plot.margin = margin(5,5,5,50),
        axis.line = element_line(color = NA))+
  xlab("")+
  ylab("Tumours with MSI\n(Proportion)")+
  scale_y_continuous(limits =c(0,0.42), expand = c(0,0))+
  guides(fill = guide_legend(title = ""))
gA
#ggsave(plot = gA, "Figures/Figure_2/Fig2_A.pdf", device = cairo_pdf, width = 7, height = 3)
