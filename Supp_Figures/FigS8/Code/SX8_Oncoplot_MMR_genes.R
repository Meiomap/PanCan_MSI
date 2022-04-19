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

# source("Supp. Figures/Oncoplot//Helpers/custom_onco.R")
# source("Supp. Figures/Oncoplot/Helpers/maftoolsHelperFuncs.R")

#####
# Make maftools plot
#####

d_all = fread("Results/All_Pathogenic_variants.tsv")

###Load variants
d = d_all%>%
  separate(col = EFFECT, sep = "&", into = c("EFFECT", NA))
d$msiStatus = NULL; d$S.ind = NULL
msi = fread("Data/MSI_status_all_samples.tsv")
msi = dplyr::select(msi, Sample_ID, msiStatus, n_rep_ind)
d = left_join(d, msi)

d = filter(d, VAF>0.20| is.na(VAF))
d = filter(d, Gnomad_AF < 0.005|is.na(Gnomad_AF))

#Load genes of interest
genes = data.frame(Gene = c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "POLE", "EXO1", "POLD1"))
# goi = fread("Data/GOI_jan26.tsv")
# goi = goi[-c(1:29),]
# goi = bind_rows(goi, data.frame(Gene = "POLD1"))
# genes = goi
# # #Somatic genes
# gs = fread("Figures/Fig3_new/somatic_fig3.tsv")
# genes_soma = gs%>%dplyr::select(Gene)%>%unique()
# 
# #Germline genes
# gs = fread("Figures/Fig3_new/germline_fig3.tsv")
# genes_germ =gs%>%dplyr::select(Gene)%>%unique()
# 
# genes = bind_rows(genes_soma, genes_germ)
d = filter(d, Gene %in% genes$Gene)


#d$EFFECT[d$EFFECT %in% c("protein_protein_contact", "structural_interaction_variant", "sequence_feature","synonymous_variant", "intron_variant", "intragenic_variant")] = "Other"
d$EFFECT[d$EFFECT %in% c("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant")] = "splice site variant"
d$EFFECT[d$EFFECT %in% c("frameshift_variant", "disruptive_inframe_insertion")] = "Frameshift variant"
d$EFFECT[d$EFFECT %in% c("start_lost", "stop_lost", "stop_gained")] = "Start/stop variant"
d$EFFECT[d$EFFECT %in% c("3_prime_UTR_variant", "5_prime_UTR_variant")] = "3'/5' UTR variant"

table(d$EFFECT)

d$EFFECT = str_replace_all(d$EFFECT, "_", " ")
table(d$EFFECT)

#########Convert to MAF
d = dplyr::select(d, Hugo_Symbol = Gene, Chromosome = CHROM, Start_Position = POS, End_Position = POS, Reference_Allele = REF, 
                  Alternative_Allele = ALT, Tumor_Sample_Barcode = Sample_ID, Variant_Classification = EFFECT, 
                  Variant_Type = Variant_Type, CADD_phred, Study, cancertype, state, n_rep_ind, msiStatus, VAF)

d$Tumor_Seq_Allele1 = d$Reference_Allele
d$Tumor_Seq_Allele2 = d$Alternative_Allele

MAF = filter(d, msiStatus == "MSI" , Hugo_Symbol %in% genes$Gene)
n_distinct(MAF$Tumor_Sample_Barcode)

#Add missing donors
msi = fread("Data/MSI_status_all_samples.tsv")
msi = filter(msi, msiStatus == "MSI", !Sample_ID %in% MAF$Tumor_Sample_Barcode)%>%dplyr::select(Tumor_Sample_Barcode = Sample_ID, everything())
d4 = bind_rows(MAF, msi)
n_distinct(d4$Tumor_Sample_Barcode)

#Add pathway info
pathways = readxl::read_excel("Supp. Figures//Oncoplot/Pathways/Goi50_2021-02-19_Eva.xlsx", sheet = 2)
pathways = filter(pathways, Genes %in% genes$Gene)
pathways$Pathway = factor(pathways$Pathway, levels = c(
  "Mismatch repair (MMR)"  , "DNA replication and damage sensing", "Excision repair"    ,    "Fanconi anemia" ,               
   "Homologous recombination (HR)"     ,       "Non-homologous end joining"      ,  
   "Nucleotide excision repair (NER)" ,  "Translesion synthesis"             
))

d4 = left_join(d4, pathways, by = c("Hugo_Symbol" = "Genes"))

n_distinct(d4$Hugo_Symbol)
n_distinct(d4$Tumor_Sample_Barcode)

table(d4$Pathway)
pathways = filter(pathways, Pathway %in% d4$Pathway)
pathways

laml.maf <<- read.maf(
  maf = d4,
  clinicalData = d4,
  vc_nonSyn = levels(as.factor(d4$Variant_Classification)),
) 

dev.off()
library(maftools)
source("Supp. Figures//Oncoplot/Helpers/Onco_25032021.R")
source("Supp. Figures//Oncoplot/Helpers/maftoolsHelperFuncs.R")
#grDevices::cairo_pdf(filename = "Article/Figures_new/Figure_2/Fig2_D_oncoplot.pdf", width = 12, height = 6,  family = "Liberation Sans")
# #png(file = "Article/Figures/Fig1/Oncoplot.png",width = 210, height = 120, units = "mm", res = 180)
#pdf("Supp. Figures/Oncoplot/S2_oncoplot_MMR_genes.pdf", width = 7.5, height = 3.4,  useDingbats=FALSE,  pointsize = 9)

pal = table(d4$Variant_Classification)
pal[1] = "#05B102FF" #Frameshift = green
pal[2] = "#026CCBFF" #Missense = blue
pal[3] = "#F51E02FF" #Splice site = red
pal[4] = "#FB9F53FF" #Start/stop = orange
pal[5] = "#9B9B9BFF" #structural interaction variant = grey
names(pal[6]) = "Multi_hit"
pal[6] = "#BA6222FF"

oncoplot(laml.maf, top = 100,
         #sepwd_samples = -1,
         #clinicalFeatures = c('cancertype', 'Study'),
         sortByAnnotation = F,
         sortByMutation = T,
         keepGeneOrder = F, GeneOrderSort = F,
         titleFontSize = 1,
         topBarData = "n_rep_ind", gene_mar = 20, colors = pal,
         pathways = pathways,
         bgCol = "grey90",
         #sampleOrder = sOrder$Tumor_Sample_Barcode,
         additionalFeature = c("state", "Germline"),
)
#dev.off()

# maftools::plotmafSummary(laml.maf)
# 
# HMF.maf = read.maf(maf = filter(d, Study == "HMF"), vc_nonSyn = levels(as.factor(d$Variant_Classification))) 
# PCAWG.maf = read.maf(maf = filter(d, Study == "HMF"), vc_nonSyn = levels(as.factor(d$Variant_Classification))) 
# vc_nonSyn = levels(as.factor(d$Variant_Classification))
# 
# coBarplot(m1 = HMF.maf, m2 = PCAWG.maf, m1Name = "HMF", m2Name = "PCAWG", genes = 2)
######
# Save figure
####
#ggsave(file="a4_output.pdf", width = 210, height = 297, units = "mm")
