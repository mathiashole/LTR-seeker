##library_______________________________________________________________________
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)
library(viridis)

##We set directory______________________________________________________________

setwd("/home/usuario/Data_Rstudio/tesina_g/tablas/LTR_candidate/")

##Read our data_________________________________________________________________

length_LTR <- read.csv("caughtUp_LTR",
                         sep = "\t",
                         stringsAsFactors = TRUE,
                         header = FALSE)

ubeity_expand <- read.csv("/home/usuario/Data_Rstudio/tesina_g/tablas/athila_candidate_clade/Id_expandSideways.txt",
                         sep = "\t",
                         stringsAsFactors = TRUE,
                         header = FALSE)

##we will need mean and standard deviation of the transposons
#gy_db_Length <- read.csv("/home/usuario/Data_Rstudio/tesina_g/tablas/dataBase_statics/Length_LTR_genome_gy_db.csv",
 #                        sep = ",", 
  #                       stringsAsFactors = T)
##we will need mean and standard deviation of the transposons

##Data peak_____________________________________________________________________

summary(length_LTR)
str(length_LTR)

summary(ubeity_expand)
str(ubeity_expand)

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

colnames(location_repEmty_LTR)<-NULL

write.table(location_repEmty_LTR, "/home/usuario/Data_Rstudio/tesina_g/tablas/LTR_candidate/location_LTR.txt", row.names = FALSE, sep = "\t", quote = FALSE, qmethod = "double")

