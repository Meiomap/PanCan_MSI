library(data.table)
library(tidyverse)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MMR_MSI_06_21//")
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21//")


ds = fread("Figures/FigX4/somatic_fig3.tsv")
dg = fread("Figures/FigX4/germline_fig3.tsv")

ds$"state" = "Somatic"
dg$"state" = "Germline"
d = bind_rows(ds,dg)

d = mutate(d, sign = ifelse(q_val < 0.05, "q \U2264 0.05", "NS"))

gA <<- ggplot(d, aes(x = sum_prop, y = -log(p_val, base = 10), 
              label = paste(Gene,": ", sum_n_MSI, " / ", sum_n_mutated, sep = ""),
       col = sign))+
  geom_point(size = 0.5+
 ggrepel::geom_text_repel(data = filter(d, q_val<0.05), size = 2.5,segment.size = 0.1,
                          nudge_y = 0.11, min.segment.length = 0.001, col = "brown", show.legend = F)+
  # geom_text(data = filter(d, q_val>0.1),
  #                          aes(label = Gene),
  #           #segment.size = 0.1,
  #                          size = 2.5, family = "Arial", col = "grey20",
  #          nudge_y = -0.1)+
  ggrepel::geom_text_repel(data = filter(d, q_val>0.05), aes(label = Gene), size = 2, segment.size = 0.1,
                           nudge_y = 0.1, min.segment.length = 0.001, max.time = 2, show.legend = F)+
  theme_bw(base_size = 7)+
  theme(
    axis.text.x = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.title = element_blank(),
    legend.text = element_text(size = 7),
    plot.margin = margin(5,5,5,50),
    axis.line = element_line(color = NA),
   strip.text = element_text(size = 7),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank()
  )+
  xlab("Mutated samples with MSI\n(Proportion)")+
  ylab("-log10(p-value)")+
  facet_wrap(~state, ncol = 2, scales = "free_x")+
#  facet_grid(rows = vars(state), scales = "free_x")+
  scale_color_manual(values = c("grey20", "brown"),
                    breaks = c("NS", "q \U2264 0.05"))+
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6), labels = c(1, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6))+
  ylab("p-value")

gA

#ggsave("Figures/Fig3_new/Fig3a.pdf", device = cairo_pdf, width = 7.4, height = 2.4)
