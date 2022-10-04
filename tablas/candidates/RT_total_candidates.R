##library to use________________________________________________________________
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)
library(viridis)
library(ggdark)

##We set directory______________________________________________________________
setwd("/home/usuario/Data_Rstudio/tesina_g/tablas/candidates/")

##Read our data_________________________________________________________________

candidate_total_rt <- read.csv("inf_o_ex_RT_candidates.csv",
                         sep = ",",
                         stringsAsFactors = T)

##Data peek_____________________________________________________________________

summary(candidate_total_rt)
str(candidate_total_rt)

##Graphics section______________________________________________________________

##Histogram_____________________________________________________________________

##Histogram length RT candiadtes athila/tat_____________________________________

candidate_RT_length <- ggplot(candidate_total_rt, aes(x= Length, fill= Length)) + 
  geom_histogram(binwidth = 20, fill="#69b3a2", alpha=0.4, position = "identity")+
  geom_density(aes(y = ..count..*437*0.11),adjust = 0.97, col = "black", fill = "#69b3a2", alpha= 0.1)+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "nº = 437", x = "Length of RT candidate") 

length(candidate_total_rt$Length)
mean(candidate_total_rt$Length)##There are many trash sequence!!!!!!!
sd(candidate_total_rt$Length)

ks.test(candidate_total_rt$Length, pnorm, mean(candidate_total_rt$Length), sd(candidate_total_rt$Length))



##Histogram %GC candiadtes athila/tat___________________________________________

candidate_RT_GC <- ggplot(candidate_total_rt, aes(x= X.GC, fill= X.GC)) + 
  geom_histogram( fill="#d646a3", alpha=0.4)+
  geom_density(aes(y =..count..), adjust = 2.5, col = "black", fill = "#d646a3", alpha= 0.15)+
  theme(legend.position="none",
        plot.title = element_text(size=11), 
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "nº = 437", x = "GC% of RT candidate") 

length(candidate_total_rt$X.GC)
mean(candidate_total_rt$X.GC)##There are many trash sequence!!!!!!!
sd(candidate_total_rt$X.GC)

ks.test(candidate_total_rt$X.GC, pnorm, mean(candidate_total_rt$X.GC), sd(candidate_total_rt$X.GC))

