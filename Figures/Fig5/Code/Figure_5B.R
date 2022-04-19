######
library(data.table)
library(tidyverse)
library(maftools)
library(ggallin)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")


d_all = fread("Results/All_Pathogenic_variants.tsv")
n_distinct(d_all$Sample_ID)
#MMRd = filter(d_all, Gene %in% c("POLE"))
#MMRd = filter(d_all, Gene %in% c("MSH6", "MSH3", "MLH1", "MSH2", "PMS1", "PMS2", "MLH3"))
MMRd = filter(d_all, Gene %in% c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "POLE", "EXO1", "POLD1"))
n_distinct(MMRd$Sample_ID)

d_all2 = filter(d_all, Sample_ID %in% MMRd$Sample_ID)
ns = n_distinct(d_all2$Sample_ID)

msi = fread("Data/MSI_status_all_samples.tsv")
msi$S.ind = msi$n_rep_ind
table(msi$msiStatus)
msi = filter(msi, Sample_ID %in% d_all2$Sample_ID)

#Load the MMR agnostic results

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

goi = fread("Data/GOI_jan26.tsv")
goi = goi[-c(1:29),]
#write.table(goi, "Data/KnijnenburgEtAl_CoreDDR.tsv",sep= "\t", col.names = T, row.names = F, quote = F)
genes = c(goi$Gene, "POLD1")

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
out$q = p.adjust(out$p, method = "fdr", n = 70+78) #Cause also testing the MMR-free cohort

out <<- out
out$col = NULL
agno = fread("Figures/FigX5//MSI_agnostic_enrichment_of_del_rep.tsv")
agno = filter(agno, q<0.05)
#  dplyr::select(Gene)%>%
#  mutate(col = "Associated to MSI with and without co-mutation with MMR genes")
# out = left_join(out, agno)

out <<- mutate(out,
             col = NA,
             col = ifelse(q<0.05 & Gene %in% agno$Gene, "Associated to MSI with and without co-mutation with MMR genes",col),
             col = ifelse(q>0.05 & Gene %in% agno$Gene, "Exclusively associated with MSI without co-mutation with MMR genes", col),
             col = ifelse(q<0.05 & !Gene %in% agno$Gene, "Exclusively associated with MSI when co-mutated with MMR", col),
             col = ifelse(Gene %in% c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "POLE", "EXO1", "POLD1") & q<0.05, 
                            yes = "MMR gene", no = col),
             col = ifelse(is.na(col), "NS", col)
)

table(out$col)

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}
cutoff = min(-log(dplyr::filter(out, q<0.05)$p, base = 10))-0.05

ggplot(out, aes(x = Dif_med, y = -log(p, base = 10),col = col, label = paste(Gene," (n=", nfg,"/",nfg+nbg ,")",
                                                        #"\n", scientific(p, digits = 2,prefix = " p \U2264 "),
                                                        sep ="")))+
  
  #Significance line
  geom_hline(yintercept = cutoff, lty = 2, col = "grey50", size = 0.2)+
  
  #Significant points
  geom_point(data = filter(out, q<0.05), size = 0.5, alpha = 0.5)+
  ggrepel::geom_text_repel(data = filter(out, q<0.05), size = 2.5, segment.size = 0.1,
                           nudge_y = 0, min.segment.length = 0.01, show.legend = F, max.overlaps = 1e3)+
  
  # #Significant points in non-MMR butt not in MMR
  # geom_point(data = filter(out, col == "Exclusively associated with MSI without co-mutation with MMR genes"), size = 0.5, alpha = 0.5)+
  # ggrepel::geom_text_repel(data = filter(out, col == "Exclusively associated with MSI without co-mutation with MMR genes"), size = 2.5, family = "Arial", segment.size = 0.1,
  #                          nudge_y = 0, min.segment.length = 0.01, show.legend = F, max.overlaps = 1e3)+
  
  #NS points
  geom_point(data = filter(out, q>0.05), size = 0.5, alpha = 0.5, col = "grey50")+
  ggrepel::geom_text_repel(data = filter(out, q>0.05), size = 1, family = "Arial", segment.size = 0.1,
                           nudge_y = 0, min.segment.length = 0.01, show.legend = F, max.overlaps = 1e3, col = "grey50")+
  
 # geom_point(data = filter(out, q>0.05), size = 0.05, alpha = 0.5, color = "grey60")+
 # geom_text(data = filter(out, q>0.05), aes(label = Gene), size = 2.5, family = "Arial", col = "grey50", show.legend = F )+
  annotate(geom = "text", x = -80, y = cutoff, label = "q = 0.05", 
           vjust = 0, size = 2.5,col = "grey50")+
  theme_bw(base_size = 7)+
  theme(
    panel.grid  = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 7),
    title = element_text(size = 7),
    plot.title = element_text(size = 7),
    legend.position = c(0.23, 0.8),
    legend.background = element_blank()
  )+
  #scale_x_log10()+
  scale_x_continuous(trans = pseudolog10_trans, breaks = c(-1e2,0,1e2, 1e3, 1e4, 1e5), labels = scales::label_comma())+
  xlab("Enrichment in median nr. of Indels in rep DNA")+
  ylab("p-value")+
  ggtitle("Enrichment of MSI across 764 samples with mutations in MMR genes, EXO1, POLD1 or POLE")+
  scale_color_manual(values = c("#39A8E0", "brown", "#662681", "goldenrod4"))+
  labs(color = "")+
  scale_y_continuous(breaks = c(0,1,2,3,4,5, 6), labels = c(1, 0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6))

#ggsave("Figures/Fig5/FigX5B.pdf", width = 7, height = 5, device = cairo_pdf)

