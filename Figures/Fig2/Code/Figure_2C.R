library(data.table)
library(tidyverse)
library(scales)


msi = fread("Data/MSI_status_all_samples.tsv")

t  = as.data.frame(table(msi[,c("msiStatus","Study")]))

gC <<- ggplot(t, aes(x = Study, y = Freq, fill = msiStatus))+
  geom_histogram(stat ="identity", position = "stack", col = "black", width=  0.2, size = 0.2)+
  xlab("")+
  ylab("No. of tumours")+
  theme_bw(base_size = 7)+
  annotate(geom = "text", x = c(1.14, 1.14, 2.14, 2.14), y = c(1500, 3450, 1500, 2550),
           label = paste("n = ",
                         c(
                           format(t$Freq[2], scientific = FALSE, big.mark = ','),
                           format(t$Freq[1], scientific = FALSE, big.mark = ','),
                           format(t$Freq[4], scientific = FALSE, big.mark = ','),
                           format(t$Freq[3], scientific = FALSE, big.mark = ',')),
                         sep =""),
           hjust = 0, size = 2.5, angle = 0)+
  scale_fill_manual(values = c("grey30", "grey90"))+
  scale_y_continuous(labels = scales::comma_format(), limits = c(0,5000),expand = c(0,0)) +
  guides(fill = guide_legend(title = "", direction = "horizontal"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.position = c(0.7,0.9),
        legend.background = element_blank(),
        axis.title = element_text(size = 7),
      #  panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = NA))
gC  

# 
# 
# 
# 
# table(msi$msiStatus)
# t  = as.data.frame(table(msi[,c("msiStatus","Study")]))
# 
# par(mfrow=c(2,1))
# 
# t1 = filter(t, Study == "PCAWG")
# slices <- t1$Freq
# lbls <- paste(t1$msiStatus, "\nn=", t1$Freq, sep ="")
# pie(slices, labels = lbls, main="#MSI patients in PCAWG")
# 
# 
# t1 = filter(t, Study == "HMF")
# slices <- t1$Freq
# lbls <- paste(t1$msiStatus, "\nn=", t1$Freq, sep ="")
# pie(slices, labels = lbls, main="#MSI patients in HMF")
# 
# 
# 
