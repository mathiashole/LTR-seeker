##library to use________________________________________________________________
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)
library(viridis)

##We set directory______________________________________________________________
setwd("/home/usuario/Data_Rstudio/tesina_g/tablas/dataBase_statics/")

##Read our data_________________________________________________________________

gy_db_Length <- read.csv("Length_LTR_genome_gy_db.csv",
                         sep = ",",
                         stringsAsFactors = T)

ltr_db_Length <- read.csv("LTR_length_gydb.csv",
                          sep = ";",
                          stringsAsFactors = T)

rt_db_Length <- read.csv("RT_length_gydb.csv",
                          sep = ";",
                          stringsAsFactors = T)

##Data peek_____________________________________________________________________

summary(gy_db_Length)
str(gy_db_Length)

summary(ltr_db_Length)
str(ltr_db_Length)

summary(rt_db_Length)
str(rt_db_Length)

##Graphics section______________________________________________________________

##Barplots______________________________________________________________________

##Barplot length genome athila/tat______________________________________________

gen_length <- ggplot(gy_db_Length, aes(x= genom, y= G.length, fill= genom)) + 
  geom_bar(stat="identity",position=position_dodge())+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Length genome Athila/Tat (pb)", x = "")+
  scale_fill_viridis_d(option = 'plasma', alpha = 0.7)

mean(gy_db_Length$G.length)
sd(gy_db_Length$G.length)

##Barplot length LTR athila/tat_________________________________________________

ltr_length <- ggplot(ltr_db_Length, aes(x= LTR, y= L.length, fill= LTR)) + 
  geom_bar(stat="identity",position=position_dodge())+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Length LTR Athila/Tat (pb)", x = "")+
  scale_fill_viridis_d(option = 'plasma', alpha = 0.7)

mean(ltr_db_Length$L.length)
sd(ltr_db_Length$L.length)

##Barplot length RT athila/tat_________________________________________________

ltr_length <- ggplot(rt_db_Length, aes(x= athila, y= rt, fill= athila)) + 
  geom_bar(stat="identity",position=position_dodge())+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Length RT Athila/Tat (AA)", x = "")+
  scale_fill_viridis_d(option = 'plasma', alpha = 0.7)

mean(rt_db_Length$rt)
sd(rt_db_Length$rt)


