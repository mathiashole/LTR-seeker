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


meet_idOf_extractor$V3 <- substr(meet_idOf_extractor$V2, start = 23, stop = (nchar(as.character(meet_idOf_extractor$V2))-2))

meet_idOf_extractor$V4 <- str_extract(meet_idOf_extractor$V3, pattern = "(.+?)[:space:]")
meet_idOf_extractor$V5 <- str_extract(meet_idOf_extractor$V3, pattern = "[:space:](.+)")

meet_idOf_extractor <- head(meet_idOf_extractor, n= 1004)
meet_idOf_extractor$V4 <- as.numeric(meet_idOf_extractor$V4)
meet_idOf_extractor$V5 <- as.numeric(meet_idOf_extractor$V5)
meet_idOf_extractor$lengthLTR <- meet_idOf_extractor$V5 - meet_idOf_extractor$V4

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

##Process data__________________________________________________________________


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
allDomains <- inner_join(GAG_domain, AP_RT_RH_INT, by="V1")
allDomains_LTR <- inner_join(meet_LTR, allDomains, by="V1")


rmDup_allDomain<- allDomains[!duplicated(allDomains$V1),]
rmDup_All <- allDomains_LTR[!duplicated(allDomains_LTR$V1),]
forOnlySeq <- rmDup_All

initScale <- rmDup_All[, c(2,3,4,5)]
sdLTR5 <- initScale[,4] - initScale[,3]
sdLTR5 <- sd(sdLTR5)
sdLTR3 <- initScale[,2] - initScale[,1]
sdLTR3 <- sd(sdLTR3)

##standard deviation of gene core length
rmDup_All[,7] <- ifelse(((rmDup_All[,8]+sd(rmDup_All[,9])) > rmDup_All[,11]), round(rmDup_All[,7]-sd(rmDup_All[,9])), rmDup_All[,7]*1)
rmDup_All[,8] <- ifelse(((rmDup_All[,8]+sd(rmDup_All[,9])) < rmDup_All[,11]), round(rmDup_All[,8]+sd(rmDup_All[,9])), rmDup_All[,8]*1)

rmDup_All[,12] <- ifelse(((rmDup_All[,12]+sd(rmDup_All[,13])) < rmDup_All[,15]), round(rmDup_All[,12]+sd(rmDup_All[,13])), rmDup_All[,12]*1)

rmDup_All[,15] <- ifelse(((rmDup_All[,16]+sd(rmDup_All[,17])) > rmDup_All[,19]), round(rmDup_All[,15]-sd(rmDup_All[,17])), rmDup_All[,15]*1)
rmDup_All[,16] <- ifelse(((rmDup_All[,16]+sd(rmDup_All[,17])) < rmDup_All[,19]), round(rmDup_All[,16]+sd(rmDup_All[,17])), rmDup_All[,16]*1)

rmDup_All[,19] <- ifelse(((rmDup_All[,19]-sd(rmDup_All[,21])) < rmDup_All[,16]), round(rmDup_All[,19]-sd(rmDup_All[,21])), rmDup_All[,19]*1)
rmDup_All[,20] <- ifelse(((rmDup_All[,19]-sd(rmDup_All[,21])) > rmDup_All[,16]), round(rmDup_All[,20]+sd(rmDup_All[,21])), rmDup_All[,20]*1)

rmDup_All[,23] <- ifelse(((rmDup_All[,23]-sd(rmDup_All[,25])) > rmDup_All[,20]), round(rmDup_All[,23]-sd(rmDup_All[,25])), rmDup_All[,23]*1)


#rmDup_All[,5] <- ifelse((rmDup_All[,5]+sdLTR5) < (rmDup_All[,7]+initScale[,3]), rmDup_All[,5]+(sdLTR5), rmDup_All[,5]*1)

rmDup_All[,4] <- ifelse((rmDup_All[,5]) < (rmDup_All[,7]+initScale[,3]), rmDup_All[,4]-(sdLTR5), rmDup_All[,4]*1)
rmDup_All[,4] <- ifelse((rmDup_All[,5]) > (rmDup_All[,7]+initScale[,3]), rmDup_All[,4]-((sdLTR5)*2), rmDup_All[,4]*1)
#rmDup_All[,4] <- ifelse((rmDup_All[,5]+(sdLTR5/2)) > (rmDup_All[,7]+initScale[,3]), rmDup_All[,4]-((sdLTR5)*1.5), rmDup_All[,4]*1)
#rmDup_All[,4] <- ifelse((rmDup_All[,5]+sdLTR5) > (rmDup_All[,7]+initScale[,3]), rmDup_All[,4]-(sdLTR5/2), rmDup_All[,4]*1)


rmDup_All[,5] <- ifelse((rmDup_All[,5]) > (rmDup_All[,7]+initScale[,3]), rmDup_All[,5]-sdLTR5, rmDup_All[,5]*1)
#rmDup_All[,5] <- ifelse((rmDup_All[,5]+((sdLTR5)/2)) > (rmDup_All[,7]+initScale[,3]), rmDup_All[,5]-((sdLTR5)/2), rmDup_All[,5]*1)
#rmDup_All[,5] <- ifelse((rmDup_All[,5]+sdLTR5) > (rmDup_All[,7]+initScale[,3]), rmDup_All[,5]-((sdLTR5)/2), rmDup_All[,5]+(sdLTR5))


##standard deviation of LTR'3 length

#rmDup_All[,3] <- ifelse((rmDup_All[,2] - sdLTR3) > (rmDup_All[,24]+initScale[,3]), rmDup_All[,2]-(sdLTR3), rmDup_All[,2]*1)
rmDup_All[,3] <- ifelse((rmDup_All[,2]) > (rmDup_All[,24]+initScale[,3]), rmDup_All[,3]+(sdLTR3), rmDup_All[,3]*1)
rmDup_All[,3] <- ifelse(rmDup_All[,2] < (rmDup_All[,24]+initScale[,3]), rmDup_All[,3]+((sdLTR3)*2), rmDup_All[,3]*1)
#rmDup_All[,3] <- ifelse((rmDup_All[,2] - ((sdLTR3)/2)) < (rmDup_All[,24]+initScale[,3]), rmDup_All[,3]+((sdLTR3)*1.5), rmDup_All[,3]*1)
#rmDup_All[,3] <- ifelse((rmDup_All[,2] - sdLTR3) < (rmDup_All[,24]+initScale[,3]), rmDup_All[,3]+((sdLTR3)/2), rmDup_All[,3]*1)

rmDup_All[,2] <- ifelse(rmDup_All[,2] < (rmDup_All[,24]+initScale[,3]), rmDup_All[,2] + sdLTR3, rmDup_All[,2]*1)
#rmDup_All[,2] <- ifelse((rmDup_All[,2] - ((sdLTR3)/2)) < (rmDup_All[,24]+initScale[,3]), rmDup_All[,2]+((sdLTR3)/2), rmDup_All[,2]*1)
#rmDup_All[,2] <- ifelse((rmDup_All[,2] - sdLTR3) < (rmDup_All[,24]+initScale[,3]), rmDup_All[,2]-((sdLTR3)/2), rmDup_All[,2]-(sdLTR3))


names(rmDup_All)[c(6,7,8,9)] = c("GAG1","GAG2","GAG3","GAG4")
names(rmDup_All)[c(10,11,12,13)] = c("AP1","AP2","AP3","AP4")
names(rmDup_All)[c(14,15,16,17)] = c("RT1","RT2","RT3","RT4")
names(rmDup_All)[c(18,19,20,21)] = c("RH1","RH2","RH3","RH4")
names(rmDup_All)[c(22,23,24,25)] = c("INT1","INT2","INT3","INT4")

endLTR <- rmDup_All[,c(1,2,3)]
endLTR[2] <- rep("LTR'3", length(rmDup_All$V1))
endLTR[c(1,3,4)] <- rmDup_All[,c(1,2,3)]
names(endLTR)[c(1,2,3,4)] = c("molecule","gene","start","end")

stLTR <- rmDup_All[,c(1,4,5)]
stLTR[2] <- rep("LTR'5", length(rmDup_All$V1))
stLTR[c(1,3,4)] <- rmDup_All[,c(1,4,5)]
names(stLTR)[c(1,2,3,4)] = c("molecule","gene","start","end")

sGAG <- rmDup_All[,c(1,7,8)]
sGAG[2] <- rep("GAG", length(rmDup_All$V1))
sGAG[c(1,3,4)] <- rmDup_All[,c(1,7,8)]
names(sGAG)[c(1,2,3,4)] = c("molecule","gene","start","end")
sGAG[c(3,4)] <- c(sGAG$start + initScale[,3], sGAG$end + initScale[,3])

sAP <- rmDup_All[,c(1,11,12)]
sAP[2] <- rep("AP", length(rmDup_All$V1))
sAP[c(1,3,4)] <- rmDup_All[,c(1,11,12)]
names(sAP)[c(1,2,3,4)] = c("molecule","gene","start","end")
sAP[c(3,4)] <- c(sAP$start + initScale[,3], sAP$end + initScale[,3])

sRT <- rmDup_All[,c(1,15,16)]
sRT[2] <- rep("RT", length(rmDup_All$V1))
sRT[c(1,3,4)] <- rmDup_All[,c(1,15,16)]
names(sRT)[c(1,2,3,4)] = c("molecule","gene","start","end")
sRT[c(3,4)] <- c(sRT$start + initScale[,3], sRT$end + initScale[,3])

sRH <- rmDup_All[,c(1,19,20)]
sRH[2] <- rep("RH", length(rmDup_All$V1))
sRH[c(1,3,4)] <- rmDup_All[,c(1,19,20)]
names(sRH)[c(1,2,3,4)] = c("molecule","gene","start","end")
sRH[c(3,4)] <- c(sRH$start + initScale[,3], sRH$end + initScale[,3])

sINT <- rmDup_All[,c(1,23,24)]
sINT[2] <- rep("INT", length(rmDup_All$V1))
sINT[c(1,3,4)] <- rmDup_All[,c(1,23,24)]
names(sINT)[c(1,2,3,4)] = c("molecule","gene","start","end")
sINT[c(3,4)] <- c(sINT$start + initScale[,3], sINT$end + initScale[,3])

catSeq <- rbind(stLTR, sGAG, sAP, sRT, sRH, sINT, endLTR)
orCatSeq <- catSeq[order(catSeq$molecule),]
orCatSeq[5] <- ifelse(orCatSeq$start < orCatSeq$end, 'forward', 'reverse')
names(orCatSeq)[5]= 'strand'  

pruebaChe <- head(orCatSeq, n= 70)


dummies <- make_alignment_dummies(
  pruebaChe,
  aes(xmin = start, xmax = end, y = molecule, id = gene),
  on = "RT"
)

ggplot(
  pruebaChe,
  aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)
) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  geom_gene_label(align = "left") +
  geom_blank(data = dummies) +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()

ggplot(
  pruebaChe,
  aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)
) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(2, "mm")) +
  geom_gene_label(align = "left") +
  geom_blank() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_brewer(palette = "Set3") +
  theme_genes()


map_dbl(rmDup_All[,c(2, 5, 9, 13, 17, 21, 25)], mean)
map_dbl(rmDup_All[,c(2, 5, 9, 13, 17, 21, 25)], sd)

map_dbl(rmDup_All[,c(2:5,7:9,11:13,15:17,19:21,23:25)], mean)
map_dbl(rmDup_All[,c(2:5,7:9,11:13,15:17,19:21,23:25)], sd)

##Process other data_____________________________________________________________

names(forOnlySeq)[c(6,7,8,9)] = c("GAG1","GAG2","GAG3","GAG4")
names(forOnlySeq)[c(10,11,12,13)] = c("AP1","AP2","AP3","AP4")
names(forOnlySeq)[c(14,15,16,17)] = c("RT1","RT2","RT3","RT4")
names(forOnlySeq)[c(18,19,20,21)] = c("RH1","RH2","RH3","RH4")
names(forOnlySeq)[c(22,23,24,25)] = c("INT1","INT2","INT3","INT4")

forOnlySeq[c(7,8,11,12,15,16,19,20,23,24)] <- c(forOnlySeq$GAG2 + initScale[,3], forOnlySeq$GAG3 + initScale[,3],
                                                forOnlySeq$AP2 + initScale[,3], forOnlySeq$AP3 + initScale[,3],
                                                forOnlySeq$RT2 + initScale[,3], forOnlySeq$RT3 + initScale[,3],
                                                forOnlySeq$RH2 + initScale[,3], forOnlySeq$RH3 + initScale[,3],
                                                forOnlySeq$INT2 + initScale[,3], forOnlySeq$INT3 + initScale[,3])

forOnlySeq$lengthLTR <- forOnlySeq$LTR3 - forOnlySeq$LTR4

getTrueSeq <- merge(forOnlySeq, meet_idOf_extractor, by = c("V1", "lengthLTR"))

getTrueSeq$V2 <- substr(getTrueSeq$V2, start = 2, stop = 19)

##Histogram length of GAG athila/tat____________________________________________

plot_GAG_length <- ggplot(rmDup_All, aes(x= GAG4, fill= GAG4)) + 
  geom_histogram( fill="#43CD80", alpha=0.4, bins = 30)+
  geom_density(aes(y =..count..*(length(rmDup_All$GAG4))*1.5), adjust = 0.6, col = "black", fill = "#43CD80", alpha= 0.15)+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Nº gene", x = "Length of GAG Athila/tat in Acca sellowiana") 

length(rmDup_All$GAG4)
mean(rmDup_All$GAG4)##There are many trash sequence!!!!!!!
sd(rmDup_All$GAG4)

ks.test(rmDup_All$GAG4, pnorm, mean(rmDup_All$GAG4), sd(rmDup_All$GAG4))


##Histogram length of AP athila/tat_____________________________________________

plot_AP_length <- ggplot(rmDup_All, aes(x= AP4, fill= AP4)) + 
  geom_histogram( fill="#43CD80", alpha=0.4, bins = 30)+
  geom_density(aes(y =..count..*10), adjust = 0.6, col = "black", fill = "#43CD80", alpha= 0.15)+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Nº gene", x = "Length of AP Athila/tat in Acca sellowiana") 

length(rmDup_All$AP4)
mean(rmDup_All$AP4)##There are many trash sequence!!!!!!!
sd(rmDup_All$AP4)

ks.test(rmDup_All$AP4, pnorm, mean(rmDup_All$AP4), sd(rmDup_All$AP4))


##Histogram length of RT athila/tat____________________________________________

plot_RT_length <- ggplot(rmDup_All, aes(x= RT4, fill= RT4)) + 
  geom_histogram( fill="#43CD80", alpha=0.4, bins = 30)+
  geom_density(aes(y =..count..*15), adjust = 2, col = "black", fill = "#43CD80", alpha= 0.15)+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Nº gene", x = "Length of RT Athila/tat in Acca sellowiana") 

length(rmDup_All$RT4)
mean(rmDup_All$RT4)##There are many trash sequence!!!!!!!
sd(rmDup_All$RT4)

ks.test(rmDup_All$RT4, pnorm, mean(rmDup_All$RT4), sd(rmDup_All$RT4))

##Histogram length of RH athila/tat_____________________________________________

plot_RH_length <- ggplot(rmDup_All, aes(x= RH4, fill= RH4)) + 
  geom_histogram( fill="#43CD80", alpha=0.4, bins = 30)+
  geom_density(aes(y =..count..*15), adjust = 2, col = "black", fill = "#43CD80", alpha= 0.15)+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Nº gene", x = "Length of RH Athila/tat in Acca sellowiana") 

length(rmDup_All$RH4)
mean(rmDup_All$RH4)##There are many trash sequence!!!!!!!
sd(rmDup_All$RH4)

ks.test(rmDup_All$RH4, pnorm, mean(rmDup_All$RH4), sd(rmDup_All$RH4))

##Histogram length of INT athila/tat____________________________________________

plot_INT_length <- ggplot(rmDup_All, aes(x= INT4, fill= INT4)) + 
  geom_histogram( fill="#43CD80", alpha=0.4, bins = 30)+
  geom_density(aes(y =..count..*103), adjust = 1, col = "black", fill = "#43CD80", alpha= 0.15)+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Nº gene", x = "Length of INT Athila/tat in Acca sellowiana") 

length(rmDup_All$INT4)
mean(rmDup_All$INT4)##There are many trash sequence!!!!!!!
sd(rmDup_All$INT4)

ks.test(rmDup_All$INT4, pnorm, mean(rmDup_All$INT4), sd(rmDup_All$INT4))


##Histogram length of LTR athila/tat____________________________________________

plot_LTR_length <- ggplot(rmDup_All, aes(x= LTR4, fill= LTR4)) + 
  geom_histogram( fill="#43CD80", alpha=0.4, bins = 30)+
  geom_density(aes(y =..count..*1000), adjust = 2, col = "black", fill = "#43CD80", alpha= 0.15)+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Nº LTR", x = "Length of LTR Athila/tat in Acca sellowiana") 

length(rmDup_All$LTR4)
mean(rmDup_All$LTR4)##There are many trash sequence!!!!!!!
sd(rmDup_All$LTR4)

ks.test(rmDup_All$LTR4, pnorm, mean(rmDup_All$LTR4), sd(rmDup_All$LTR4))


##total histogram of genes______________________________________________________

grid.arrange(plot_GAG_length, plot_AP_length, plot_RT_length, plot_RH_length, plot_INT_length, plot_LTR_length, nrow=3, ncol=2)

##Write start and end of genes__________________________________________________

toExtractor_RT <- getTrueSeq[,c(27,16,17)]
colnames(toExtractor_RT)<-NULL
write.table(toExtractor_RT, "/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/toExtractor_RT.txt", row.names = FALSE, sep = "\t", quote = FALSE, qmethod = "double")

toExtractor_RNaseH <- rmDup_All[,c(1,19,20)]
colnames(toExtractor_RNaseH)<-NULL
write.table(toExtractor_RNaseH, "/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/toExtractor_RNaseH.txt", row.names = FALSE, sep = "\t", quote = FALSE, qmethod = "double")

toExtractor_INT <- rmDup_All[,c(1,23,24)]
colnames(toExtractor_INT)<-NULL
write.table(toExtractor_INT, "/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/toExtractor_INT.txt", row.names = FALSE, sep = "\t", quote = FALSE, qmethod = "double")

toExtractor_GAG <- rmDup_All[,c(1,7,8)]
colnames(toExtractor_GAG)<-NULL
write.table(toExtractor_GAG, "/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/toExtractor_GAG.txt", row.names = FALSE, sep = "\t", quote = FALSE, qmethod = "double")

toExtractor_AP <- rmDup_All[,c(1,19,20)]
colnames(toExtractor_AP)<-NULL
write.table(toExtractor_AP, "/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/toExtractor_AP.txt", row.names = FALSE, sep = "\t", quote = FALSE, qmethod = "double")
