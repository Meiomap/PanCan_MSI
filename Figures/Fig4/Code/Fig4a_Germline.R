###
## Association between MSI and germline/somatic mutations
## 
###
library(data.table)
library(tidyverse)
library(scales)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MMR_MSI_06_21/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")


all = fread("Data/LOF_25_all.tsv")
all$Study = NULL
all$cancertype = NULL
all$primaryTumorLocation = NULL
n_distinct(all$Sample_ID)
msiStatus = fread("Data/MSI_status_all_samples.tsv")
msiStatus$Donor_ID = NULL
n_distinct(msiStatus$Sample_ID)
all = left_join(all, msiStatus, by = c("Sample_ID"))
table(is.na(all$msiStatus))
d = filter(all, Germline == T, Mutation_level >= 2)
#d = filter(d, !is.na(msiStatus))
table(duplicated(d[,c("Sample_ID", "Gene", "msiStatus")]))
d = distinct(d, Sample_ID, Gene, msiStatus, Study, cancertype)
n_distinct(d$Sample_ID) 
#HOW MANY MSI VS HOW ANY MSS ARE MUTATED IN ANY ONE GENE


#Sum over cancertypes
d2 = d%>%
  dplyr::select(-cancertype)%>%
  group_by(Study, Gene, msiStatus)%>%
  mutate(
    n = n_distinct(Sample_ID)
  )%>%
  distinct(Study, Gene, msiStatus, .keep_all = T)%>%
  pivot_wider(names_from = msiStatus, values_from = n)%>%
  group_by(Study, Gene)%>%
  mutate(
    n_MSI = sum(MSI, na.rm = T), #To summarize within study, sometimes two rows for same gene, but not always
    n_MSS = sum(MSS, na.rm = T)
  )%>%
  distinct(Study, Gene, .keep_all = T)%>%
  mutate(
    n_mutated = sum(n_MSI, n_MSS),
    prop_MSI = n_MSI/n_mutated,
    max_hit = max(n_MSI, n_MSS)
  )%>%
  ungroup()%>%
  distinct(Study, Gene, .keep_all = T)%>%
  group_by(Study)%>%
  filter(n_mutated >= 3)%>%
  distinct()%>%
  mutate(
    median_prop = ifelse(Study == "PCAWG", 22/2560, 86/3497)
  )

#Select genes for foreground..
goi = fread("Data/GOI_jan26.tsv")
goi = goi[-c(1:29),]
goi = bind_rows(goi, data.frame(Gene = "POLD1"))
#goi = filter(goi, Pathway %in% c("Mismatch Repair (MMR)", "Non-homologous End Joining (NHEJ)", "Fanconi Anemia (FA)"))

d3 = filter(d2, Gene %in% goi$Gene)%>%
  filter(n_MSI>=0, n_mutated >= 0)

#d3 = d2

# ggplot(data = d3, aes(x = prop_MSI))+
#   geom_histogram(fill = "white", color = "black")+
#   geom_vline(xintercept = overall_prop, lty = 2, col = "red")

#hist(d3$n_MSI)

#Binom-test per sample
for(i in 1:nrow(d3)){
  ti = binom.test(x = d3$n_MSI[i], n = d3$n_mutated[i], p = d3$median_prop[i], alternative = "greater")
  d3$"p_val"[i] = ti$p.value
}

#Add missing gene/study combinations to fix p-val evaluation
ex = tidyr::expand_grid(Gene = d3$Gene, Study = d3$Study)%>%
  distinct(Gene, Study)%>%
  mutate(ID = paste(Gene, Study))%>%
  filter(!ID %in% paste(d3$Gene, d3$Study))%>%
  dplyr::select(-ID)%>%
  mutate(p_val = 1)
d3 = bind_rows(d3, ex)

source("Code/Helper_functions/sumlog.R")
d4 = dplyr::select(d3, Gene, Study, n_MSI, n_mutated, prop_MSI, p_val)%>%
 # pivot_wider(names_from= Study, values_from = c(p_val, n_MSI, n_mutated))%>%
  group_by(Gene)%>%
  mutate(
    p_val2 = p_val,
    p_val = ifelse(n_distinct(p_val)>1,
                   sumlog(c(p_val))$p, 
                   sum(p_val, na.rm = T)),
    p_val = ifelse(p_val == 2, 1, p_val)
  )%>%
  pivot_wider(names_from = Study, values_from = c(n_MSI, n_mutated, prop_MSI, p_val2))
#d4[is.na(d4)]=0
d4 = d4%>%
  mutate(
    n_mutated_PCAWG = ifelse(is.na(n_mutated_PCAWG), 0 , n_mutated_PCAWG),
    n_mutated_HMF = ifelse(is.na(n_mutated_HMF), 0 , n_mutated_HMF),
    n_MSI_PCAWG = ifelse(is.na(n_MSI_PCAWG), 0 , n_MSI_PCAWG),
    n_MSI_HMF = ifelse(is.na(n_MSI_HMF), 0 , n_MSI_HMF),
    prop_MSI_PCAWG = ifelse(is.na(prop_MSI_PCAWG), 0 , prop_MSI_PCAWG),
    prop_MSI_HMF = ifelse(is.na(prop_MSI_HMF), 0 , prop_MSI_HMF),
    sum_n_MSI = n_MSI_PCAWG + n_MSI_HMF,
    sum_n_mutated = n_mutated_PCAWG + n_mutated_HMF,
    sum_prop = sum_n_MSI/sum_n_mutated
  )

  
d4$q_val = p.adjust(d4$p_val, method = "fdr", n = nrow(d4)) #Cause we also test for the germline variatns

write.table(d4, "Figures/FigX4/germline_fig3.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

d4 = mutate(d4, 
            sign = ifelse(p_val < 0.1, "p \U2264 0.10", "NS"),
            sign = ifelse(p_val < 0.05, "p \U2264 0.05", sign)
)

table(d4$sign)
dout <<-d4
d4 = filter(d4, !is.na(sum_n_mutated) & sum_n_MSI > 0)

write.table(d4, "Data/Tmp/Germline_goi.tsv", sep = "\t", col.names = T, row.names = F, quote  = F)

#Barplot
gB_Germ <<- ggplot(data = d4, aes(x = reorder(paste(Gene,": ", sum_n_MSI, " / ", sum_n_mutated, sep = ""), -sum_prop), y = sum_prop))+
  geom_histogram(stat = "identity", aes(fill = sign), position = "stack",color = "black",
                 width = 0.5, show.legend =T, size = 0.2)+
  geom_text(data = filter(d4, p_val < 0.01), aes(label = scientific(p_val, digits = 2,prefix = " p \U2264 ")),
            vjust = 0.5, hjust = 0, angle = 90, size = 2.5, family = "Arial")+
  geom_text(data = filter(d4, p_val > 0.01), aes(label = paste(" p \U2264 ", round(p_val, 2), sep = "")),
            vjust = 0.5, hjust = 0, angle = 90, size = 2.5, family = "Arial")+
  ggtitle("Pathogenic germline variants in DDR genes")+
  ylab("Prop. samples with MSI")+theme_bw(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        title = element_text(size = 7),
        axis.title = element_text(size = 7),
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = NA))+
 scale_fill_manual(values = c("grey", "lightblue", "darkblue"),
                  breaks = c("NS", "p \U2264 0.10", "p \U2264 0.05"))+
  xlab("")+
  ylab("Prop. Mutated samples with MSI")+
  guides(fill = guide_legend(title = ""))+
  scale_y_continuous(limits = c(0,0.35))
  # annotate(geom = "label", x = nrow(d3), y = overall_prop, hjust = 1,
  #          label = "Prop. patients with MSI across all DDR genes", col = "brown")


gB_Germ

# goi = d3%>%
#   dplyr::select(state, Gene, q_val, p_val, sign, prop_MSI)%>%
#   mutate(
#     prop_MSI_deviance = prop_MSI - overall_prop
#   )
#write.table(goi, "Results/somatic_goi.tsv", sep ="\t", col.names = T, row.names = F)
