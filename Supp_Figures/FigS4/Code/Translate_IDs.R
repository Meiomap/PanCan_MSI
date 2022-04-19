library(tidyverse)
library(data.table)

translate_id = function(d = d){
d = d
info = fread("repository_1610444206.tsv")
info = info%>%
  rowwise()%>%
  separate(`File Name`, sep = "\\.", into = c("File Name", NA))

table(d$ID %in% info$`Object ID`)
table(d$ID %in% info$`File Name`)

pcawg1 = filter(d, ID %in% info$`Object ID`)%>%
  left_join(dplyr::select(info, ID = `Object ID`, Sample_ID = `Sample ID`, Donor_ID = `ICGC Donor`))%>%
  distinct()

pcawg2 = filter(d, ID %in% info$`File Name`)%>%
  left_join(dplyr::select(info, ID = `File Name`, Sample_ID = `Sample ID`, Donor_ID = `ICGC Donor`))%>%
  distinct()
                  
pcawg = bind_rows(pcawg1, pcawg2)
hmf = filter(d, !ID %in% pcawg$ID)
hmf$"Donor_ID" = str_remove_all(hmf$ID, "T*$"); hmf$"Donor_ID" = str_remove_all(hmf$Donor_ID, "TI*$"); hmf$"Donor_ID" = str_remove_all(hmf$Donor_ID, "TII*$")
hmf$"Donor_ID" = str_remove_all(hmf$Donor_ID, "TIV$")
hmf$"Sample_ID" = hmf$ID
d = bind_rows(pcawg, hmf)

white = fread("all_donors.tsv")
d = filter(d, Sample_ID %in% white$Sample_ID)
d$ID = NULL
d = left_join(d, white)
return(d)
}
