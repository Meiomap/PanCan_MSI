######
library(data.table)
library(tidyverse)
library(maftools)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")


d_all = fread("Results/All_Pathogenic_variants.tsv")
n_distinct(d_all$Sample_ID)
#MMRd = filter(d_all, Gene %in% c("POLE"))
#MMRd = filter(d_all, Gene %in% c("MSH6", "MSH3", "MLH1", "MSH2", "PMS1", "PMS2", "MLH3"))
MMRd = filter(d_all, Gene %in% c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "POLE", "EXO1", "POLD1"))

d_all2 = filter(d_all, !Sample_ID %in% MMRd$Sample_ID) #Only search among MMR proficient samples
#d_all2 = d_all
ns = n_distinct(d_all2$Sample_ID)

msi = fread("Data/MSI_status_all_samples.tsv")
msi = filter(msi, Sample_ID %in% d_all2$Sample_ID)
msi$S.ind = msi$n_rep_ind
gp = function(gene = "EXO1"){
  exo = filter(d_all2, Gene %in% gene)
  dfg = filter(msi, Sample_ID %in% exo$Sample_ID)
  dfg = mutate(dfg, status = paste("Mutated, n=", n_distinct(dfg$Sample_ID), sep =""))
  dbg = filter(msi, !Sample_ID %in% exo$Sample_ID)
  dbg = mutate(dbg, status = paste("Not Mutated, n=", n_distinct(dbg$Sample_ID), sep =""))
  d = bind_rows(dfg, dbg)
  
  tmp = data.frame(
    Gene = gene,
    nfg = n_distinct(dfg$Sample_ID),
    nbg = n_distinct(dbg$Sample_ID),
    p =  wilcox.test(dfg$S.ind, dbg$S.ind)$p.value,
    Dif_med = median(dfg$S.ind)- median(dbg$S.ind)
  )
  tmp$n = tmp$nfg + tmp$nbg
  tmp
}
gp(gene = "BRCA2")
#genes1 = c("MSH3", "MSH6", "MLH1", "PMS2", "MRE11A", "RAD50", "NBN", "POLE", "ATR","TOPBP1", "TOP3A", "PALB2", "XRCC4")

goi = fread("Data/GOI_jan26.tsv")
goi = goi[-c(1:29),]
goi = bind_rows(goi, data.frame(Gene = "POLD1"))
#write.table(goi, "Data/KnijnenburgEtAl_CoreDDR.tsv",sep= "\t", col.names = T, row.names = F, quote = F)
genes = goi$Gene

for(g in genes){
  try({
  t = gp(gene = g)
  if(g == genes[[1]]){
    out = t
  }else{
    out = bind_rows(out, t)
  }
  }, next)
}

#out = filter(out, nfg >= 10)
out$q = p.adjust(out$p, method = "fdr", n = 70+78)

out <<- out
#write.table(out, "Figures/FigX5/MSI_agnostic_enrichment_of_del_rep.tsv", sep ="\t", col.names = T, row.names = F)

# reverselog_trans <- function(base = exp(1)) {
#   trans <- function(x) -log(x, base)
#   inv <- function(x) base^(-x)
#   trans_new(paste0("reverselog-", format(base)), trans, inv, 
#             log_breaks(base = base), 
#             domain = c(1e-100, Inf))
# }

ggplot(out, aes(x = Dif_med, y = -log(p, base = 10), label = paste(Gene," (n=", nfg,")",
                                                        #"\n", scientific(p, digits = 2,prefix = " p \U2264 "),
                                                        sep ="")))+
  geom_point(size = 0.05, alpha = 0.5)+
  ggrepel::geom_text_repel(data = filter(out, q<0.05, nfg>= 10), size = 2.5, family = "Arial", segment.size = 0.1,
                           nudge_y = 1, min.segment.length = 0.01, col = "darkblue")+
  geom_text(data = filter(out, q>0.05 | nfg<10), aes(label = Gene), size = 1.5, family = "Arial", col = "grey50",
            nudge_y = -0.3, col = "lightblue")+
  geom_hline(yintercept = -log(max(filter(out, q<0.05, nfg>= 10)$p), base = 10), lty = 2, col = "grey50", size = 0.2)+
  annotate(geom = "text", x = 1, y = -log(max(filter(out, q<0.05, nfg>= 10)$p), base = 10)+0.1, label = "p = 0.01", 
           vjust = 0, size = 2.5, family = "Arial",col = "grey50")+
  theme_bw(base_size = 7, base_family = "Arial")+
  theme(
    panel.grid  = element_blank(),
    axis.title = element_text(size = 7),
    title = element_text(size = 7),
    plot.title = element_text(size = 7)
  )+
#  scale_x_continuous(trans = reverselog_trans, breaks = c(-1e2,0,1e3,5e3, 1e4))+
  xlab("Enrichment in median no. of Indels in rep DNA")+
  ylab("-log(Wilcoxon two-tailed p-value)")+
  ggtitle(paste("Enrichment of MSI across", ns, "samples with no pathogenic mutations in MMR genes, EXO1, POLD1 or POLE"))

#ggsave("Figures/Fig4//volcano2.pdf", width = 6, height = 2, device = "pdf")

