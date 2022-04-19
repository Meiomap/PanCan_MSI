library(data.table)
library(tidyverse)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/MSI_MMR/MMR_MSI/")

d = fread("Results/Variants.tsv")
d = filter(d, !is.na(Sample_ID))

#Add cancertype
ct = fread("Data/Cancertypes_jan2021.tsv")
d = left_join(d, ct)%>%
  filter(!is.na(cancertype))%>%
  distinct()

#Rename CT
source("Figures/Figure_1/Cancertypes_rename.R")
d = rename_ct(d)

#Filter by impact
d = filter(d, Mutation_level>=6)
d = filter(d, VAF>0.20 | is.na(VAF))
d = filter(d, Gnomad_AF < 0.005|is.na(Gnomad_AF))

#Filter MSI status of samples
msi = fread("Data/MSI_status_all_samples.tsv")
d = left_join(d, dplyr::select(msi, Sample_ID, S.ind, msiStatus))
d = filter(d, !is.na(msiStatus))
#d = filter(d, msiStatus == "HIGH")

#Annotate state
d$"state"[grep("Germ", d$Variant_Type)] = "Germline"
d$"state"[grep("Germ", d$Variant_Type, invert = T)] = "Somatic"
table(is.na(d$state))

#Remove hypermutated variants — already removed all 50+ per study actually
d = d%>%
  group_by(CHROM, POS, REF, ALT)%>%
  mutate(n = n())%>%
  filter(n < 200)%>%
  dplyr::select(-n)

table(duplicated(d$Sample_ID)) #5972— meaning that some samples have no variants

d = rename(d, Gene = GENE)

#Remove the two non-consent HMF samples
white = fread("Data/210127_extra_hmfIds.tsv")
hmf = d[grep("SA", d$Sample_ID, invert = T),]
hmf = filter(hmf, Sample_ID %in% white$sampleId)
PCAWG = d[grep("SA", d$Sample_ID, invert = F),]
d = bind_rows(hmf, PCAWG)

d = filter(d, Gene != "TP53")

write.table(d, "Results/All_Pathogenic_variants.tsv", sep = "\t", col.names = T, row.names = F)
