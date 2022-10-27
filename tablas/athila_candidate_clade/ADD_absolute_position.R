##library_______________________________________________________________________
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)
library(viridis)
library(gggenes)
library(RColorBrewer)
library(ggdark)


##We set directory______________________________________________________________

setwd("/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/")

##Read our data_________________________________________________________________

##Reverse transcriptase_________________________________________________________

meet_RT <- read.csv("blast_RT_retro_LTR_athila_candidate_1e-5",
                    sep = "\t",
                    stringsAsFactors = TRUE,
                    header = FALSE)

if ( meet_RT[1,1] == 'Warning: [blastx] Examining 5 or more matches is recommended' ) {
  meet_RT <- meet_RT[2:length(meet_RT[,1]),]
}


##RNaseH________________________________________________________________________

meet_RNaseH <- read.csv("blast_RNaseH_retro_LTR_athila_candidate_1e-5",
                        sep = "\t",
                        stringsAsFactors = TRUE,
                        header = FALSE)

if ( meet_RNaseH[1,1] == 'Warning: [blastx] Examining 5 or more matches is recommended' ) {
  meet_RNaseH <- meet_RNaseH[2:length(meet_RNaseH[,1]),]
}

#AP_____________________________________________________________________________
meet_AP <- read.csv("blast_AP_retro_LTR_athila_candidate_1e-5",
                    sep = "\t",
                    stringsAsFactors = TRUE,
                    header = FALSE)

if ( meet_AP[1,1] == 'Warning: [blastx] Examining 5 or more matches is recommended' ) {
  meet_AP <- meet_AP[2:length(meet_AP[,1]),]
}

#Integrase______________________________________________________________________
meet_INT <- read.csv("blast_INT_retro_LTR_athila_candidate_1e-5",
                     sep = "\t",
                     stringsAsFactors = TRUE,
                     header = FALSE)

if ( meet_INT[1,1] == 'Warning: [blastx] Examining 5 or more matches is recommended' ) {
  meet_INT <- meet_INT[2:length(meet_INT[,1]),]
}

#Gag____________________________________________________________________________
meet_GAG <- read.csv("blast_GAG_retro_LTR_athila_candidate_1e-5",
                     sep = "\t",
                     stringsAsFactors = TRUE,
                     header = FALSE)

if ( meet_GAG[1,1] == 'Warning: [blastx] Examining 5 or more matches is recommended' ) {
  meet_GAG <- meet_GAG[2:length(meet_GAG[,1]),]
}

##LTR___________________________________________________________________________

meet_LTR <- read.csv("/home/usuario/Data_Rstudio/tesina_g/tablas/LTR_candidate/caughtUp_LTR",
                     sep = "\t",
                     stringsAsFactors = TRUE,
                     header = FALSE)
meet_LTR <- meet_LTR[,c(1,7,8,9,10)]

names(meet_LTR)[1] = "V1"
names(meet_LTR)[c(2,3,4,5)] = c("LTR2","LTR3","LTR4","LTR5")

meet_LTR$V1 <- str_c(meet_LTR$V1,"_")

##idOf_extractor fasta__________________________________________________________

meet_idOf_extractor <- read.csv("/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/idOf-extractor_need-for-phylogenetics-domain",
                     sep = "\t",
                     stringsAsFactors = TRUE,
                     header = FALSE)
meet_idOf_extractor$V2 <- substr(meet_idOf_extractor$V1, start = 2, stop = 20)

meet_idOf_extractor$V3 <- meet_idOf_extractor[,1]
meet_idOf_extractor$V1 <- meet_idOf_extractor[,2]
meet_idOf_extractor$V2 <- meet_idOf_extractor[,3]
meet_idOf_extractor <- meet_idOf_extractor[,c(1,2)]

##Read our data_________________________________________________________________

length_LTR <- read.csv("/home/usuario/Data_Rstudio/tesina_g/tablas/LTR_candidate/caughtUp_LTR",
                       sep = "\t",
                       stringsAsFactors = TRUE,
                       header = FALSE)

ubeity_expand <- read.csv("/home/usuario/Data_Rstudio/tesina_g/tablas/athila_candidate_clade/Id_expandSideways.txt",
                          sep = "\t",
                          stringsAsFactors = TRUE,
                          header = FALSE)

tab_id_db_gy <- read.csv("/home/usuario/Data_Rstudio/tesina_g/tablas/athila_candidate_clade/Id_length_RT_gy_total.csv",
                         sep = "\t",
                         stringsAsFactors = TRUE)
all_id_Rt.blast <- tab_id_db_gy
all_id_Rt.blast$RT4 <- ifelse(tab_id_db_gy$S > tab_id_db_gy$E, tab_id_db_gy$S - tab_id_db_gy$E, tab_id_db_gy$E - tab_id_db_gy$S)

names(all_id_Rt.blast)[c(1, 2, 3)] = c("V1", "V2", "V3")

all_id_Rt.blast$V1 <- str_c(all_id_Rt.blast$V1,"__")

##Data peak_____________________________________________________________________

summary(meet_RNaseH)
str(meet_RNaseH)

summary(meet_RT)
str(meet_RT)

summary(meet_AP)
str(meet_AP)

summary(meet_GAG)
str(meet_GAG)

summary(meet_INT)
str(meet_INT)

summary(meet_idOf_extractor)
str(meet_idOf_extractor)

summary(length_LTR)
str(length_LTR)

summary(ubeity_expand)
str(ubeity_expand)


##Process data__________________________________________________________________

##absolute length of transposons________________________________________________

##Process data__________________________________________________________________

names(ubeity_expand)[2] = "V15"
names(ubeity_expand)[3]= "V16"

length_LTR$V1 <- substring(length_LTR$V1, 1, 17)
join_df<-full_join(length_LTR, ubeity_expand, by="V1")

LtrAndUbeity <- join_df[,c(1,7,8,9,10,15,16)]

LtrAndUbeity$V17 <-ifelse(LtrAndUbeity$V15 > LtrAndUbeity$V16, LtrAndUbeity$V15 - LtrAndUbeity$V16, LtrAndUbeity$V16 - LtrAndUbeity$V15)

LtrAndUbeity$V18 <- ifelse(LtrAndUbeity$V7 < LtrAndUbeity$V9, LtrAndUbeity$V7, LtrAndUbeity$V9)

LtrAndUbeity$V19 <- ifelse(LtrAndUbeity$V8 > LtrAndUbeity$V10, LtrAndUbeity$V8, LtrAndUbeity$V10)

location_LTR <- LtrAndUbeity[,c(1,9,10)]

location_repEmty_LTR <- location_LTR[!duplicated(location_LTR),]

location_repEmty_LTR[!duplicated(location_repEmty_LTR$V19),]

location_repEmty_LTR$V1 <- str_c(location_repEmty_LTR$V1,"_")

location_repEmty_LTR$V4 <- location_repEmty_LTR$V19 - location_repEmty_LTR$V18

##Domains_______________________________________________________________________


domains <- list(meet_RT, meet_INT, meet_GAG, meet_AP, meet_RNaseH)

RT_domain <- data.frame(domains[[1]][1], domains[[1]][2], domains[[1]][7], domains[[1]][8])
INT_domain <- data.frame(domains[[2]][1], domains[[2]][2], domains[[2]][7], domains[[2]][8])
GAG_domain <- data.frame(domains[[3]][1], domains[[3]][2], domains[[3]][7], domains[[3]][8])
AP_domain <- data.frame(domains[[4]][1], domains[[4]][2], domains[[4]][7], domains[[4]][8])
RNaseH_domain <- data.frame(domains[[5]][1], domains[[5]][2], domains[[5]][7], domains[[5]][8])

RT_domain$V15 <- ifelse(RT_domain$V7 > RT_domain$V8, RT_domain$V7 - RT_domain$V8, RT_domain$V8 - RT_domain$V7)
INT_domain$V15 <- ifelse(INT_domain$V7 > INT_domain$V8, INT_domain$V7 - INT_domain$V8, INT_domain$V8 - INT_domain$V7)
GAG_domain$V15 <- ifelse(GAG_domain$V7 > GAG_domain$V8, GAG_domain$V7 - GAG_domain$V8, GAG_domain$V8 - GAG_domain$V7)
AP_domain$V15 <- ifelse(AP_domain$V7 > AP_domain$V8, AP_domain$V7 - AP_domain$V8, AP_domain$V8 - AP_domain$V7)
RNaseH_domain$V15 <- ifelse(RNaseH_domain$V7 > RNaseH_domain$V8, RNaseH_domain$V7 - RNaseH_domain$V8, RNaseH_domain$V8 - RNaseH_domain$V7)


RH_INT <- inner_join(RNaseH_domain, INT_domain, by="V1")
RT_RH_INT <- inner_join(RT_domain, RH_INT, by="V1")
AP_RT_RH_INT <- inner_join(AP_domain, RT_RH_INT, by="V1")
GAG_AP_RT_RH_INT <- inner_join(GAG_domain, AP_RT_RH_INT, by="V1")
LTR_GAG_AP_RT_RH_INT <- inner_join(meet_LTR, GAG_AP_RT_RH_INT, by="V1")


rmDup_LTR_GAG_AP_RT_RH_INT<- LTR_GAG_AP_RT_RH_INT[!duplicated(LTR_GAG_AP_RT_RH_INT$V1),]

initScale <- rmDup_LTR_GAG_AP_RT_RH_INT[, c(2,3,4,5)]
sdLTR5 <- initScale[,4] - initScale[,3]
sdLTR5 <- sd(sdLTR5)
sdLTR3 <- initScale[,2] - initScale[,1]
sdLTR3 <- sd(sdLTR3)

names(rmDup_LTR_GAG_AP_RT_RH_INT)[c(6,7,8,9)] = c("GAG1","GAG2","GAG3","GAG4")
names(rmDup_LTR_GAG_AP_RT_RH_INT)[c(10,11,12,13)] = c("AP1","AP2","AP3","AP4")
names(rmDup_LTR_GAG_AP_RT_RH_INT)[c(14,15,16,17)] = c("RT1","RT2","RT3","RT4")
names(rmDup_LTR_GAG_AP_RT_RH_INT)[c(18,19,20,21)] = c("RH1","RH2","RH3","RH4")
names(rmDup_LTR_GAG_AP_RT_RH_INT)[c(22,23,24,25)] = c("INT1","INT2","INT3","INT4")

rmDup_LTR_GAG_AP_RT_RH_INT[c(7,8,11,12,15,16,19,20,23,24)] <- c(rmDup_LTR_GAG_AP_RT_RH_INT$GAG2 + initScale[,3], rmDup_LTR_GAG_AP_RT_RH_INT$GAG3 + initScale[,3],
                                                rmDup_LTR_GAG_AP_RT_RH_INT$AP2 + initScale[,3], rmDup_LTR_GAG_AP_RT_RH_INT$AP3 + initScale[,3],
                                                rmDup_LTR_GAG_AP_RT_RH_INT$RT2 + initScale[,3], rmDup_LTR_GAG_AP_RT_RH_INT$RT3 + initScale[,3],
                                                rmDup_LTR_GAG_AP_RT_RH_INT$RH2 + initScale[,3], rmDup_LTR_GAG_AP_RT_RH_INT$RH3 + initScale[,3],
                                                rmDup_LTR_GAG_AP_RT_RH_INT$INT2 + initScale[,3], rmDup_LTR_GAG_AP_RT_RH_INT$INT3 + initScale[,3])

position_absolut_RT <- merge(rmDup_LTR_GAG_AP_RT_RH_INT, all_id_Rt.blast, by = "V1")


