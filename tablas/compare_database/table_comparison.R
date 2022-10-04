##library to use________________________________________________________________
library(tidyverse)
library(ggplot2)
library(viridis)
library(dplyr)
library(gridExtra)
library(RColorBrewer)

##We set directory_______________________________________________________________
setwd("/home/usuario/Data_Rstudio/tesina_g/tablas/compare_df/")

##Read our data
tab_gypsy <- read.csv("Id_s_f_gypsy_tab.csv", sep = "\t", stringsAsFactors = TRUE)

tab_rex <- read.csv("Id_s_f_rex_tab.csv", sep = "\t", stringsAsFactors = TRUE)

tab_two <- read.csv("Id_s_f_two_tab.csv", sep = "\t", stringsAsFactors = TRUE)

##Compare diferent dataframe____________________________________________________

#We highlight differences_______________________________________________________

#anti_join Keeps the elements of the first table for which there is no 
#information in the second
gy_rex <- anti_join(tab_gypsy, tab_rex)
rex_gy <- anti_join(tab_rex, tab_gypsy)

gy_two <- anti_join(tab_gypsy, tab_two)
two_gy <- anti_join(tab_two, tab_gypsy)

rex_two <- anti_join(tab_rex, tab_two)
two_rex <- anti_join(tab_two, tab_rex)

bind_cols(gy_rex, gy_two, rex_gy, rex_two, two_gy, two_rex)

df_diff <- data.frame(
  "diff" = c(rep("gy", "rex", "two", 3)),
  "gy.rex" = c(0, 8, 176),
  "rex" = c(2, 0, 49),
  "two" = c(126, 5, 0)
)

ggplot(df_diff, aes(x= diff, y= , fill= nombre, alpha= .4)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  theme(legend.position="none",
        plot.title = element_text(size=11),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.background = element_rect(fill = "white",
                                        colour = "grey50")) +
  labs(y = "Largo de elementos Athilas", x = "")+
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(colorCount))


##We set directory_______________________________________________________________
setwd("/home/usuario/Data_Rstudio/tesina_g/tablas/athila_candidate_clade/")

tab_id_db_gy <- read.csv("Id_length_RT_gy_total.csv",
                         sep = "\t", stringsAsFactors = TRUE)

tab_clade_at <- read.csv("Id_new.clade.csv",
                         sep = "\t", stringsAsFactors = TRUE)

#I get start and end of our clade_______________________________________________

duplicated(tab_clade_at)
duplicated(tab_id_db_gy)

Id_length_clade <- merge(tab_clade_at, tab_id_db_gy, by = "ID")

write.csv(Id_length_clade, "/home/usuario/Data_Rstudio/tesina_g/tablas/athila_candidate_clade/Id_clade_length_extractor.csv", row.names = FALSE)

SorR_n <- rep(11500, 288) #half the length of our LTR



