# 2 Proportion of mono-nucleotide and di-nucleotide stretches per gene (~30% mono, 8%di in average)
library(tidyverse)
library(data.table)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")

#Note, pretty sure i only count mono-nuclleotide events! But might be a fair proxy for di-nucleotide events as well
#Length of indel flanks
d = fread("Data/Del_ins_counts//Length_of_indels.tsv")
d = mutate(d, length = ifelse(group == "Dinucleotide", length * 2, length))
#d$length = d$length+1 #Cause we need to count the single base of the indel as well
hist(d$length)
table(d$length)
d2 = d%>%
  group_by(length)%>%
  mutate(n = n())%>%
  distinct(length, n)%>%
  filter(length >= 3)

bet = data.frame(length = c(24, 26, 38, 40, 44))
g1 = ggplot(d2, aes(x = length, y = n))+
  geom_histogram(stat = "identity")+
  theme_classic2(base_size = 7)+
  scale_y_continuous(labels = scales::comma_format())+
  xlab("Repeat sequence length\n(No. of bases)")+
  ylab("Mono- and dinucleotide indels\n(735 DDR genes, ~6,000 tumours)")+
  scale_x_continuous(limits = c(0,70), breaks = c(3,10,20,40,60))+
  geom_rug(sides = "b", col = "red", data = bet, aes(x = length), inherit.aes = F, show.legend = T)+
  annotate(label = "Bethesda panel\nmarker lengths", col = "red", geom = "text", size = 2.5, x = 44, 
           y = 5e3, angle = 0, hjust = 0)
g1
  
#Load all mono stretches
mo = fread("Data/Del_ins_counts/mononucleotide_stretches_all_genes_with_introns.tsv")
mo = mutate(mo, group = "Mo")
mo = filter(mo, Base != "N")
mo2 = mo%>%
  group_by(Length)%>%
  mutate(n = n())%>%
  distinct(Length, n)

g2 = ggplot(mo2, aes(x = Length, y = n))+
  geom_histogram(stat = "identity")+
  theme_classic2(base_size = 7)+
  scale_y_continuous(labels = scales::comma_format())+
  xlab("Repeat sequence length\n(No. of bases)")+
  ylab("Mono- and dinucleotide repeat sequences\n(735 DDR genes, hg19)")+
  scale_x_continuous(limits = c(0,70), breaks = c(3,10,20,40,60))
g2

# #Rate for each length
# d3 = left_join(dplyr::rename(mo2, n_rep_seq = n, length = Length), dplyr::rename(d2, n_indels = n))%>%
#   mutate(rate = n_indels/n_rep_seq)
# 
# d3[is.na(d3)] = 0
# bet = data.frame(length = c(24, 26, 38, 40, 44))
# g3 = ggplot(d3, aes(x = length, y = rate))+
#   geom_histogram(stat = "identity")+
#   theme_classic2(base_size = 7)+
#   scale_y_continuous(labels = scales::comma_format())+
#   xlab("Length of repeat sequence\n(No. of bases)")+
#   ylab("Rate of indels per repeat sequence")+
#   scale_x_continuous(limits = c(0,70), breaks = c(3,10,20,40,60))+
#   geom_rug(sides = "b", col = "red", data = bet, aes(x = length), inherit.aes = F, show.legend = T)+
#   annotate(label = "Bethesda panel\nmarker lengths", col = "red", geom = "text", size = 2.5, x = 44, 
#            y = 1.3, angle = 0, hjust = 0)
# 
# g3


#Bet rates
#24: 0.07 (We observe an event in 7%)
#26: 0.02 
#38: 0
#40: 0
#44: 0


library(patchwork)

layout = "
AA
BB
"s

out =g2 + g1 + 
  patchwork::plot_layout(design = layout, guides = "collect") +
  plot_annotation(
    title = "Supplementary figure X7",
    subtitle = '',
    #caption = ''
    tag_levels = 'a'
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    text = element_text( size = 7),
    title =  element_text( size = 7),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7)
  )
out

setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/Supp. Figures/Del_ins_counts_S4_to_S8/")
ggsave(plot = out, "../SX7.pdf", device = cairo_pdf, width = 4, height = 4)


