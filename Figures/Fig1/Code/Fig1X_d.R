library(tidyverse)
library(scales)

setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")

d = readxl::read_excel("Figures/FigX1_Amrutas_figure//Screen_F1d.xlsx", sheet = 1) 
d$Plot[d$strain == "rad27∆"] = 1
d = dplyr::filter(d, Plot == 1)
d2 = d%>%
  pivot_longer(cols = c(EXP1, EXP2, EXP3))%>%
  group_by(strain, Group)%>%
  mutate(
    value = value*1e5,
    mean = mean(value),
    std  = sd(value),
    Group = ifelse(Group == "Control", "gene∆ + vector only", "gene∆ + h*MLH1*"))

#Sig test
# for(s in unique(filter(d, Plot == 1)$strain)){
#   tmp_test = dplyr::filter(d2, strain == s, Group == "gene∆ + h*MLH1*", Plot == 1)
#   tmp_control = dplyr::filter(d2, strain == "wt", Group == "gene∆ + h*MLH1*", Plot == 1)
#   t = t.test(x = tmp_test$value, y = tmp_control$value, alternative = "greater")
#   d2$pVal[d2$strain == s] = t$p.value
# }
# 
# for(s in unique(filter(d, Plot == 2)$strain)){
#   tmp_test = dplyr::filter(d2, strain == s, Group == "gene∆ + h*MLH1*", Plot == 2)
#   tmp_control = dplyr::filter(d2, strain == s, Group == "gene∆ + vector only", Plot == 2)
#   t = t.test(x = tmp_test$value, y = tmp_control$value, alternative = "greater")
#   d2$pVal[d2$strain == s] = t$p.value
# }
#d2$q = p.adjust(d2$pVal)

d2$strain = factor(d2$strain, levels = rev(unique(d$strain)))
d2$Group = factor(d2$Group, levels = (c("gene∆ + vector only", "gene∆ + h*MLH1*")))

p1 = ggplot(dplyr::filter(d2, Plot == 1), aes(x = mean, y = strain, fill = Group, color = Group))+
  geom_col(stat = "identity", position=position_dodge(0.5), width = 0.5, 
           size = 0.2, color = "black")+
  geom_errorbarh(aes(x = mean, xmin = mean-std, xmax = mean+std), 
                 height = 0.5, color = "black", position=position_dodge(width=0.5), size = 0.1)+
  geom_vline(xintercept = d2$mean[d2$strain == "wt" & d2$Group == "gene∆ + h*MLH1*"][1], lty = 2, size = 0.2, col = "grey20")+
  
  geom_point(aes(x = value), color = "black",
             position=position_dodge(width=0.5), show.legend = F, size = 0.3, shape = 4)+
  xlab(bquote("Mutation rate per cell/per generation"~(x10^-5)))+
  theme_bw(base_size = 7)+
  theme(
    panel.grid = element_blank(),
      legend.text = ggtext::element_markdown(size = 7),
      axis.text.y = element_text(face = "italic",color = "black", size = 7),
     # legend.position = c(0.7,0.92),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.2, "cm"),
      axis.text.x = element_text(size = 7)
    )+
  scale_fill_manual(values = rev(c("red", "grey85")))+
  scale_y_discrete(breaks = rev(unique(d$strain)))+
  ylab("")+
  guides(fill = guide_legend(reverse=T))+
 # scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8), position = "bottom")+
  ggbreak::scale_x_break(c(0.8,1), scales = c(0.3,0.8), ticklabels = c(1,2,3))
p1 



#library(patchwork)
#out = p1 + p2 + plot_layout(nrow = 2, heights = c(5,1))
#out
ggsave(plot = p1, filename = "Figures/FigX1_Amrutas_figure/Fig1X_d.pdf", device = cairo_pdf,height = 4, width = 4)           

