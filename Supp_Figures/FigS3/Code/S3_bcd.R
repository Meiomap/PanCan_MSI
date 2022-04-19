library(tidyverse)
library(scales)

setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")

d = readxl::read_excel("Supp. Figures/Amrutas/S3tables.xlsx", sheet = 1) 
d = pivot_longer(d, cols = c(ex1, ex2))%>%
  mutate(Strain = paste(Strain1, Strain2, sep = " + "))%>%
  group_by(Strain, Figure)%>%
  mutate(
    value = value * 1e5,
    mean = mean(value), 
    std = sd(value))%>%
  ungroup()
d$Strain1 = factor(d$Strain1, levels = unique(d$Strain1))
d$Strain2 = factor(d$Strain2, levels = unique(d$Strain2))
d1 = filter(d, Figure == "B")

p1 = ggplot(d1, aes(x = Strain1, y = mean, fill = Strain2))+
  geom_col(stat = "identity", position=position_dodge(0.5), width = 0.5, 
           size = 0.2, color = "black")+
  geom_errorbar(aes(y = mean, ymin = mean-std, ymax = mean+std), 
                 height = 0.5, color = "black", position=position_dodge(width=0.5), size = 0.1)+
  geom_point(aes(y = value), color = "black",
             position=position_dodge(width=0.5), show.legend = F, size = 0.3, shape = 4)+
  ylab(bquote("Mutation rate"~(x10^-5)))+
  theme_bw(base_size = 7)+
  theme(
    panel.grid = element_blank(),
      legend.text = ggtext::element_markdown(size = 7),
      axis.text.x = ggtext::element_markdown(size = 7),
      legend.position = c(0.8,0.92),
      legend.background = element_blank(),
      legend.title = element_blank(),
      legend.key.height = unit(0.2, "cm"),
      legend.key.width = unit(0.2, "cm"),
      axis.text.y = element_text(size = 7)
    )+
  scale_fill_manual(values = rev(c("red", "grey85")))+
  xlab("")
p1 

d2 = filter(d, Figure == "C")
d2$Strain1 = factor(d2$Strain1, levels = unique(d2$Strain1))
p2 = ggplot(d2, aes(x = Strain1, y = mean, fill = Strain2))+
  geom_col(stat = "identity", position=position_dodge(0.5), width = 0.5, 
           size = 0.2, color = "black", show.legend = F)+
  geom_errorbar(aes(y = mean, ymin = mean-std, ymax = mean+std), 
                height = 0.5, color = "black", position=position_dodge(width=0.5), size = 0.1)+
  geom_point(aes(y = value), color = "black",
             position=position_dodge(width=0.5), show.legend = F, size = 0.3, shape = 4)+
  ylab(bquote("Mutation rate"~(x10^-5)))+
  theme_bw(base_size = 7)+
  theme(
    panel.grid = element_blank(),
    legend.text = ggtext::element_markdown(size = 7),
    axis.text.x = ggtext::element_markdown(size = 7),
    axis.text.y = element_text(size = 7)
  )+
  scale_fill_manual(values = rev(c("blue", "grey85")))+
  xlab("")
p2 


layout = "
AA
B#
"
library(patchwork)
out = p1 + p2 + plot_layout(design = layout)
out

ggsave(plot = out, filename = "Supp. Figures/Amrutas/S3_bc.pdf", device = cairo_pdf,height = 3.5, width = 3.3)           

