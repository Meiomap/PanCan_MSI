######
# Figure 1 — Summary figure 
# - 1. Overview of samples per study, cancertype and metastasis
# - 2. Overview of DDR deficiency across samples
# - 3. Overview of types of mutational signals
######
library(data.table)
library(tidyverse)
library(maftools)
library(extrafont)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")

#####
# Make maftools plot
#####

d_all = fread("Results/All_Pathogenic_variants.tsv")
d = filter(d_all, VAF>0.20| is.na(VAF))
d = filter(d, Gnomad_AF < 0.005|is.na(Gnomad_AF))

#d$Gene[d$Gene %in% c( "MSH3", "MSH6")] = "MSH3+MSH6"
#########
d = dplyr::select(d, Hugo_Symbol = Gene, Chromosome = CHROM, Start_Position = POS, End_Position = POS, Reference_Allele = REF, 
                  Alternative_Allele = ALT, Tumor_Sample_Barcode = Sample_ID, Variant_Classification = EFFECT, 
                  Variant_Type = Variant_Type, CADD_phred, Study, cancertype, state, S.ind, msiStatus, VAF)

d$Tumor_Seq_Allele1 = d$Reference_Allele
d$Tumor_Seq_Allele2 = d$Alternative_Allele

#Filter by MSI
d2 = filter(d, msiStatus == "MSI")

genes = c("PMS1","PMS2","MSH2", "MSH3", "MSH6", "MLH1", "MLH3","POLE", "POLD1", "EXO1", "RAD50",  "TOP3A")
#genes = distinct(d2, Hugo_Symbol)$Hugo_Symbol
for(i in 1:length(genes)){
  g1 = genes[i]
  tmp1 = filter(d2, Hugo_Symbol == g1)%>%distinct(Tumor_Sample_Barcode, .keep_all = T)
  tmp_out = data.frame(Gene1 = rep(g1, length(genes)), 
                       n_Gene1 = rep(n_distinct(tmp1$Tumor_Sample_Barcode), length(genes)),
                       Gene2 = rep(NA, length(genes)),
                       n_Gene2 = rep(NA, length(genes)),
                       n_Comutated = rep(NA, length(genes))
  )
    for(j in 1:length(genes)){
      g2 = genes[j]
      tmp2 = filter(d2, Hugo_Symbol == g2)%>%distinct(Tumor_Sample_Barcode, .keep_all = T)
      tmp_out$Gene2[j] = g2 
      tmp_out$n_Gene2[j] = n_distinct(tmp2$Tumor_Sample_Barcode) 
      tmp_out$n_Comutated[j] = sum(tmp2$Tumor_Sample_Barcode %in% tmp1$Tumor_Sample_Barcode)
    }
  if(i == 1){
    out = tmp_out
  }else{
    out = bind_rows(out, tmp_out)
  }
}

#AAAAAALLL this code just to get the damn order right : P
out2 = out%>%
  mutate(
    n_Comutated = ifelse(Gene1 == Gene2, NA, n_Comutated),
    n_samples = n_Gene1+n_Gene2-n_Comutated,
    prop = n_Comutated/n_samples)

out2 <<- out2

out2$Gene1 = factor(out2$Gene1, levels =genes)

out2 = out2%>%
  arrange(by = Gene1)%>%
  mutate(Lab_x = as.factor(paste(Gene1, " (",n_Gene1,")", sep ="")))

out2$Gene1 = factor(out2$Gene1, labels = unique(out2$Lab_x),
                    levels = genes)

out2$Gene2 = factor(out2$Gene2, levels = )

out2 = out2%>%
  arrange(by = (Gene2))%>%
  mutate(Lab_y = as.factor(paste(Gene2, " (",n_Gene2,")", sep ="")))

out2$Gene2 = factor(out2$Gene2, labels = unique(out2$Lab_y),
                    levels =  genes)

ggplot(out2, aes(x = Gene1, 
                 y= Gene2, 
                 fill = n_Comutated))+
  geom_tile(color = "grey60")+
  geom_text(data = filter(out2, !is.na(prop)),aes(label = paste(round(prop, 2)*100,"%", sep = "")), size = 1.7, family = "Arial")+
  scale_fill_gradient2(low = "white", mid = "white", high = "brown", midpoint = 4, na.value = "grey")+
  theme_bw(base_size = 7, base_family = "Arial")+
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(size = 7),
    title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.key.height = unit(0.1, "cm")
  )+
  guides(fill = guide_colorbar(title = "No. of co-mutated samples"))+
  ggtitle("Co-mutation of MMR genes and genes with q<0.10", subtitle = "Across 183 MSI samples")+
  xlab("")+
  ylab("")+
  scale_x_discrete(position = "bottom")+
  scale_y_discrete(limits = rev)

#ggsave("Figures/Fig4/coMut.pdf", device = "pdf", width = 5, height = 3)  


#Filter by dMMR
# #d2 = filter(d, Tumor_Sample_Barcode %in% MMRd$Sample_ID)
# # 
# d2$Variant_Classification = "Pathogenic"
# laml.maf <<- read.maf(
#   maf = d2,
#   clinicalData = d2,
#   vc_nonSyn = levels(as.factor(d2$Variant_Classification)),
# )
# 
# t = somaticInteractions(maf = laml.maf,
#                 genes = c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "MRE11A", "RAD50", "NBN", "POLE", "ATR", "TOPBP1", "TOP3A", "PALB2", "XRCC4"),
#                 #genes = c("dMMR","MRE11A", "RAD50", "NBN", "POLE", "ATR", "TOPBP1", "TOP3A", "PALB2", "XRCC4"),
#                 top = 25,
#                     pvalue = c(0.05, 0.1),
#                     showCounts = T,
#                     showSigSymbols = F)


