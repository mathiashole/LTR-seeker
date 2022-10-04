##library_______________________________________________________________________
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)
library(viridis)

##We set directory______________________________________________________________

setwd("/home/usuario/Data_Rstudio/tesina_g/tablas/athila_candidate_clade/")

##Read our data_________________________________________________________________

tab_clade_at <- read.csv("Id_new.clade.csv",
                         sep = "\t",
                         stringsAsFactors = TRUE)
tab_id_db_gy <- read.csv("Id_length_RT_gy_total.csv",
                         sep = "\t",
                         stringsAsFactors = TRUE)

##we will need mean and standard deviation of the transposons
gy_db_Length <- read.csv("/home/usuario/Data_Rstudio/tesina_g/tablas/dataBase_statics/Length_LTR_genome_gy_db.csv",
                         sep = ",", 
                         stringsAsFactors = T)
##we will need mean and standard deviation of the transposons

##Data peak_____________________________________________________________________

summary(tab_clade_at)
str(tab_clade_at)

summary(tab_id_db_gy)
str(tab_id_db_gy)

summary(gy_db_Length)
str(gy_db_Length)

##Process data__________________________________________________________________


  
  for(i in duplicated(tab_clade_at)){
    if(i != FALSE){
    print(i)
    } else {
    print('NO')
    }
  }

  
  
  for(i in duplicated(tab_id_db_gy)){
    if(i != FALSE){
      print(i)
    } else {
      print('NO')
    } 
  }

  
##we are going to integrate the location to our clade___________________________
  
Id_length_clade <- merge(tab_clade_at, tab_id_db_gy, by = "ID")

##Data peak_____________________________________________________________________

summary(Id_length_clade)  
str(Id_length_clade)

##we add mean and standard deviation to the candidates__________________________

meanOfDb <- mean(gy_db_Length$G.length)
sdOfDb <- sd(gy_db_Length$G.length)

only_length <- Id_length_clade[,2:3]

Id_length_clade$Se <- ifelse(Id_length_clade$S > Id_length_clade$E, round(meanOfDb + sdOfDb), round(-meanOfDb-sdOfDb))
Id_length_clade$Es <- ifelse(Id_length_clade$E < Id_length_clade$S, round(-meanOfDb-sdOfDb), round(meanOfDb + sdOfDb))
Id_length_clade$S <- Id_length_clade$S + Id_length_clade$Se
Id_length_clade$E <- Id_length_clade$E + Id_length_clade$Es
Id_length_clade$S <- ifelse(Id_length_clade$S < 0 , (Id_length_clade$S*0) + 1, Id_length_clade$S*1)
Id_length_clade$E <- ifelse(Id_length_clade$E < 0 , (Id_length_clade$E*0) + 1, Id_length_clade$E*1)
#Id_length_clade$ID <- paste0(rep(">",length(Id_length_clade$ID)), Id_length_clade$ID)

length(Id_length_clade$ID)
str(only_length)

colnames(Id_length_clade)<-NULL

write.table(Id_length_clade[,1:3], "/home/usuario/Data_Rstudio/tesina_g/tablas/athila_candidate_clade/Id_expandSideways.txt", row.names = FALSE, sep = "\t", quote = FALSE, qmethod = "double")

(nth <- paste(1:12, c("st", "nd", "rd", rep("th", 9))))
