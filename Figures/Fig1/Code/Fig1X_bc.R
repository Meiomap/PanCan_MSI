library(tidyverse)
library(data.table)
library(extrafont)
  
 # setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

d = readxl::read_excel("Figures/FigX1_Amrutas_figure//DAVID_F1_15_03_22.xls", sheet = 1) 
d$`List Total` = round(d$Count/(d$'%'/100),0) #DAVID filtered away genes not recognized, but they were still considered in the statistical tests. This way I can get the true number of inputted genes
d = d%>%
  tidyr::separate(Term, into = c("GOterm", "Process"), sep = "~")%>%
  dplyr::rename(Observed = Count)%>%
  mutate(Expected = (`Pop Hits`/`Pop Total`) * `List Total`) #159 is the number of genes i
d$Process = firstup(d$Process)
#  d$PValue = NULL
 # d$FDR = NULL
  for(i in 1:nrow(d)){
    tmp = d[i,]
    
    observed = tmp$Observed
    genelist = tmp$`List Total`
    Genes_in_process = tmp$`Pop Hits`
    total_genes = tmp$`Pop Total`
    
    dat = data.frame(
      "Experiment" = c(observed-1, genelist-observed+1),
      "Population" = c(Genes_in_process, total_genes-Genes_in_process),
      row.names = c("in_process", "Not_in_process"),
      stringsAsFactors = FALSE
    )
    dat
    mosaicplot(dat,
               main = "Mosaic plot",
               color = TRUE
    )  
    d$p[i] = fisher.test(dat, alternative = "greater")$p.value
  }

d=dplyr::select(d, GOterm, Process, Expected, Observed, `Fold Enrichment`, pValue = p, Genes)
hist(d$pValue)
d$FDR = p.adjust(d$pValue, method = "fdr")
increased = d

d1 = dplyr::filter(d, pValue<.05)

d2 = d1%>%
    mutate(Obs2 = Observed,
           p2 = pValue)%>%
    pivot_longer(cols = c(Observed, Expected))%>%
    mutate(pValue = ifelse(name == "Expected", NA, pValue))
d2 = mutate(d2, name = ifelse(name == "Observed", "Observed in trial", "Expected at random"))
  
g1 = ggplot(d2, aes(y = reorder(Process, -p2), x = value, fill = name, 
                     label = paste("  p \u2264" ,scales::scientific(pValue, digits = 2)), sep = ""))+
              #  label = paste(" ", p2, sep ="")))+
    geom_histogram(stat = "identity", width = 0.5, color = "black", size = 0.2,position = "dodge")+
    geom_text(data = na.omit(d2), hjust = 0, size = 2)+
    theme_bw(base_size = 7)+
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(size = 7),
      plot.title = element_text(size = 7),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.2, "cm"),
      legend.direction = "vertical",
      axis.text = element_text(size = 7)
    )+
    scale_x_continuous(limits = c(0,25), expand = c(0,0))+
    ylab("")+
    scale_fill_manual(values = rev(c("black","white")))+
    ggtitle("Increased")+
    xlab("No. of genes")
  
g1
  
#Decreased
d = readxl::read_excel("Figures/FigX1_Amrutas_figure//DAVID_F1_15_03_22.xls", sheet = 2) 
d$`List Total` = round(d$Count/(d$'%'/100),0) #DAVID filtered away genes not recognized, but they were still considered in the statistical tests. This way I can get the true number of inputted genes
d = d%>%
  tidyr::separate(Term, into = c("GOterm", "Process"), sep = "~")%>%
  dplyr::rename(Observed = Count)%>%
  mutate(Expected = (`Pop Hits`/`Pop Total`) * `List Total`) #159 is the number of genes i
d$Process = firstup(d$Process)

#  d$PValue = NULL
# d$FDR = NULL
for(i in 1:nrow(d)){
  tmp = d[i,]
  
  observed = tmp$Observed
  genelist = tmp$`List Total`
  Genes_in_process = tmp$`Pop Hits`
  total_genes = tmp$`Pop Total`
  
  dat = data.frame(
    "Experiment" = c(observed-1, genelist-observed+1),
    "Population" = c(Genes_in_process, total_genes-Genes_in_process),
    row.names = c("in_process", "Not_in_process"),
    stringsAsFactors = FALSE
  )
  dat
  mosaicplot(dat,
             main = "Mosaic plot",
             color = TRUE
  )  
  d$p[i] = fisher.test(dat, alternative = "greater")$p.value
}

d=dplyr::select(d, GOterm, Process, Expected, Observed, `Fold Enrichment`, pValue = p, Genes)
hist(d$pValue)
d$FDR = p.adjust(d$pValue, method = "fdr")
decreased = d


d1 = dplyr::filter(d, pValue<.05)

d2 = d1%>%
  mutate(Obs2 = Observed,
         p2 = pValue)%>%
  pivot_longer(cols = c(Observed, Expected))%>%
  mutate(pValue = ifelse(name == "Expected", NA, pValue))
d2 = mutate(d2, name = ifelse(name == "Observed", "Observed in trial", "Expected at random"))

g2 = ggplot(d2, aes(y = reorder(Process, -p2), x = value, fill = name, 
                    label = paste("  p \u2264" ,scales::scientific(pValue, digits = 2)), sep = ""))+
  #  label = paste(" ", p2, sep ="")))+
  geom_histogram(stat = "identity", width = 0.5, color = "black", size = 0.2,position = "dodge", show.legend = F)+
  geom_text(data = na.omit(d2), hjust = 0, size = 2)+
  theme_bw(base_size = 7)+
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 7),
    axis.text = element_text(size = 7)
  )+
  scale_x_continuous(limits = c(0,25), expand = c(0,0))+
  ylab("")+
  scale_fill_manual(values = rev(c("black","white")))+
  ggtitle("Decreased")+
  xlab("No. of genes")

g2
library(patchwork)
out = g1 + g2 + plot_layout(nrow = 2, heights = c(1,1.9)) & theme(
 # axis.text.y = element_text(angle = 45, hjust = 1, vjust = 1)
)
out
increased = as.data.frame(increased)
increased = arrange(increased, pValue)

decreased = as.data.frame(decreased)
decreased = arrange(decreased, pValue)

#Write to excel
library(xlsx)
wb = loadWorkbook("Figures/FigX1_Amrutas_figure/David_forrmatted.xlsx")
xlsx::removeSheet(wb, sheetName = "Decreased")
xlsx::removeSheet(wb, sheetName = "Increased")
Decreased = xlsx::createSheet(wb, sheetName = "Decreased")
Increased = xlsx::createSheet(wb, sheetName = "Increased")
xlsx::addDataFrame(x = decreased, sheet = Decreased, row.names = FALSE, byrow = F)
xlsx::addDataFrame(x = increased, sheet = Increased, row.names = FALSE, byrow = F)
xlsx::saveWorkbook(wb, "Figures/FigX1_Amrutas_figure/David_forrmatted.xlsx")

ggsave(plot = out, filename = "Figures/FigX1_Amrutas_figure//bothPlots.pdf", device = cairo_pdf,height = 4.5, width = 4.5)           
