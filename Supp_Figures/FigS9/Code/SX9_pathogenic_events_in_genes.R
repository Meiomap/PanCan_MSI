library(tidyverse)
library(data.table)
library(BSgenome.Hsapiens.UCSC.hg19)
library(ggpubr)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI_06_21/")

#Pre-load data
Patho = fread("Results/All_Pathogenic_variants.tsv")


msi = fread("Data/MSI_status_all_samples.tsv")
#Patho = filter(Patho, !Sample_ID %in% dplyr::filter(msi, msiStatus == "MSI")$Sample_ID)#Filter to MSI high or low
#Find mono and dinucleotide flanks
per_case_summary = function(gene = "BRCA2"){
 # gene = "DCLRE1C"
  
  #Extract Pathogenic events in the gene among the right patients
  Patho_tmp <<- dplyr::filter(Patho, Gene == gene)
  n_events = nrow(Patho_tmp)
  
  ##Identify indels
  rep = Patho_tmp%>%
    mutate(
      length_ref = nchar(REF),
      length_alt = nchar(ALT),
      TYPE = "SNP",
      TYPE = ifelse(length_ref>length_alt, "DEL", TYPE),
      TYPE = ifelse(length_ref<length_alt, "INS", TYPE),
      length = ifelse(TYPE == "INS", length_alt-1, length_ref-1),
      seq = ifelse(TYPE == "DEL", sub('.', '', REF),  sub('.', '', ALT))
    )%>%
    filter(TYPE != "SNP")
  
  if(nrow(rep)>0){
    
    #####MONO EVENTS!
    #Add flanks to categorize the variants
    chr <- paste0("chr",rep$CHROM)
    position <- rep$POS
    #alleles <- paste0("[", pcawg$HGVS_P, "]")
    offset <- 100
    del_length = rep$del_length
    
    rep$m5_seq = paste(getSeq(Hsapiens,chr,position-offset,position-1), sep = "") #Lets get both sides
    rep$m3_seq = paste(getSeq(Hsapiens,chr,position+1,position+offset), sep ="")

    #Identify indels with only one base (might be repeated many times)
    Mo = rep%>%
      rowwise()%>%
      mutate(
        first_base = unlist(str_split(seq,""))[1], #First base in each insertion (could be any of the bases in the insertion)
        n_rep_first_char = str_count(seq, first_base), #Count how many time this base reoccurs in the event
        min_string = strrep(first_base, 3) #Copy this first base too 3 times, as this is the minimal string tto search for
      )%>%
      dplyr::filter(n_rep_first_char == length) #Only keep cases where all bases are he samee
    
    #Test if min-string occurs at start of M3 or M5
    Mo$m3_Mo = str_detect(Mo$m3_seq, paste0("^",Mo$min_string))
    Mo$m5_Mo = str_detect(Mo$m5_seq, paste0(Mo$min_string,"$"))
    Mo$first_base = NULL; Mo$n_Mo_first_char = NULL; Mo$min_string = NULL; Mo$FORMAT = NULL
    Mo = mutate(Mo, m3_Mo = ifelse(m3_Mo == T & m5_Mo == T, F, m3_Mo)) #In cases where both are true to not get douyble count
   # Mo = filter(Mo, m3_Mo == T |m5_Mo == T )
    
    out_tmp = data.frame(n_rep_mo = sum(Mo$m3_Mo)+sum(Mo$m5_Mo))
  
  ####DINUCLEOTIDE EVENTS
      
  #Remoove deletions with only one base (might be repeated many times)
  Di = rep%>%
    rowwise()%>%
    mutate(
      first_base = unlist(str_split(seq,""))[1], #First base in each insertion (could be any of the bases in the insertion)
      n_rep_first_char = str_count(seq, first_base), #Count how many time this base reoccurs in the event
      min_string = strrep(first_base, 3) #Copy this first base too 3 times, as this is the minimal string tto search for
    )%>%
  dplyr::filter(n_rep_first_char != length) #Only keep cases where all bases are he different

    
    #Test if min-string occurs at start of M3
    Di$m3_Di = str_detect(Di$m3_seq, paste0("^",Di$min_string))
    Di$m5_Di = str_detect(Di$m5_seq, paste0(Di$min_string,"$"))
    Di$first_base = NULL; Di$n_Di_first_char = NULL; Di$min_string = NULL; Di$FORMAT = NULL; Di$del_length = NULL #Cleanup
    Di = mutate(Di, m3_Di = ifelse(m3_Di == T & m5_Di == T, F, m3_Di)) #In cases where both are true to not get douyble count
    
    out_tmp$n_rep_di = sum(Di$m3_Di)+sum(Di$m5_Di)
    out_tmp$n_pathogenic_events = n_events; out_tmp$n_indels = nrow(rep); out_tmp$Gene = gene
    
  }else{
    out_tmp = data.frame(
      Gene = gene,
      n_rep_di = 0,
      n_rep_mo = 0,
      n_pathogenic_events = n_events,
      n_indels = 0
    )
  }
  
  out = out_tmp
      
  return(out)
}

tt = per_case_summary(gene = "MSH6")
tt


#Run through all genes
genes = readxl::read_excel("../../DDR_project_May_2021/Feedback_from_reviewers/Split_by_dataset/Supp_tables/Sup_tables.xlsx",
                           sheet = 3)

for(i in 1:nrow(genes)){
  mod = genes[i,]
  tmp = per_case_summary(gene = mod$GENE)
  if(i == 1){
    result = tmp
  }else{
    result = bind_rows(result, tmp)
  }
}

#Export
write.table(result, "Supp. Figures/Del_ins_counts_S4_to_S8/MSS_rate_of_patho_rep_per_gen.tsv", sep ="\t",
            col.names = T, row.names = F)

#Make figure and table
msi_pats = fread("Supp. Figures/Del_ins_counts_S4_to_S8/MSI_108_rate_of_patho_rep_per_gen.tsv")
msi_pats$group = "MSI, 108 tumours"
mss_pats = fread("Supp. Figures/Del_ins_counts_S4_to_S8/MSS_rate_of_patho_rep_per_gen.tsv")
mss_pats$group = "MSS, 5,949 tumours"
d = bind_rows(msi_pats, mss_pats)%>%
  mutate(n_rep = n_rep_di+n_rep_mo, rate = n_rep / n_pathogenic_events)%>%filter(n_pathogenic_events >=0)

d$group = factor(d$group, levels = c("MSS, 5,949 tumours","MSI, 108 tumours"))
d$lim[d$group == "MSS, 5,949 tumours"] = 400
d$lim[d$group != "MSS, 5,949 tumours"] =300

#Add MMR genes?
mmr = c("PMS2", "MSH3", "MSH6", "MLH1", "MSH2", "PMS1", "MLH3", "POLE", "EXO1", "POLD1")
#Rep stress
rep = c("TOPBP1","CHEK1","RAD50","TOP3A", "ATR","PRKDC", "MRE11A", "CUL5", "REV3L", "BRCA2", "TP53BP1", 
        "SHPRH", "ERCC6", "FANCD2", "POLL","GEN1", "ATM" )

d = mutate(d,
           gro = "All DDR genes",
           gro = ifelse(Gene %in% mmr, "MMR genes", gro),
           gro = ifelse(Gene %in% rep, "Genes associated with elevated MSI", gro))

ggplot(filter(d, n_pathogenic_events > 4), aes(x = rate, color = gro))+
  geom_histogram(color = "black", size = 0.2)+
  facet_grid(rows = vars(group), scales = "free_y")+
  xlab("Rate of pathogenic indels\n(No. of repeat indels / no. of pathogenic events)")+
  ylab("No. of genes")+
  theme_bw(base_size = 7)+
  geom_text(data = distinct(d, group, .keep_all = T),size = 2.5, aes(label = group, y = lim), x = 0.1, hjust= 0)+
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 7),
    legend.title = element_blank(),
    plot.title = element_text(size = 7),
    panel.grid = element_blank()
  )+
  geom_rug(sides = "b")+
  ggrepel::geom_text_repel(
    data = filter(d, Gene %in% c(mmr, rep), rate > 0.05), 
    aes(label = paste("italic('",Gene, "')", sep =""), y = 1), 
    nudge_y = 100,
    nudge_x = 0.1,
    angle = 0, 
    hjust = 0, 
    vjust = 0, 
    size = 2.5, 
    #direction = "y", 
    segment.alpha = 0.2,
    parse = TRUE,
    max.overlaps = 1e2)+
  scale_color_manual(values = c("black", "blue", "red"))+
  ggtitle("Supplementary Figure S8")
  

ggsave("Supp. Figures/SX8_2.pdf", device = cairo_pdf, width = 7, height = 5)
#write.table(d,)


  

