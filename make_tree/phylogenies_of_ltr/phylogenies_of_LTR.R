##Library we need_______________________________________________________________
library(ape)
library(BiocManager)
library(tidytree)
library(phylotools)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(TDbook)
#devtools::install_github("YuLab-SMU/ggtree") for install ggtree
library(ggtree)
library(treeio)
library(dplyr)
library(shadowtext)
library(ggnewscale)
library(ggrepel)

##First, we set the directory  where to work____________________________________

setwd("/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/aligment_abs_all/")

##Then, we read tree____________________________________________________________

t_RT_tree <- read.tree("out_aligment_RT_abs.pep.treefile")

##Tidytree package provides as_tibble method to convert the phylo to a tidy data frame

gy_df <-as_tibble(t_RT_tree) #convert phylo to data frame
#And as_phylo() convert data frame to phylo
#_______________________________________________________________________________

#The dplyr verbs can be applied to tbl_tree directly to manipulate tree data.
#In addition, tidytree provides several verbs to filter related nodes, including
#child(), parent(), offspring(), ancestor(), sibling() and MRCA()
#_______________________________________________________________________________

#groupClade() method accepts an internal node or a vector of internal nodes
#to add grouping information of selected clade/clades.

#groupOTU will trace back from input nodes to most recent common ancestor.
#_______________________________________________________________________________

##Outgroup with root() implemented ape package__________________________________

true_RTtree_outgroup <- root(gy_df_tree, outgroup = "RT_CRM")

##phylogenetics tree____________________________________________________________

p1 <- ggtree(true_RTtree_outgroup) +
  geom_nodelab(
    mapping = aes(
      x = branch
    ),
    nudge_y = 0.36
  ) +
  xlim(-.1, 4.5) +
  scale_size_continuous(range = c(3, 10)) +
  geom_tiplab(
    offset = .14, 
  ) + 
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu")) +
  theme(legend.position = "right") #Little visible


##asigno numero de nodos
ggtree(true_RTtree_outgroup) + geom_text(aes(label=node), hjust=-0.2) + 
  geom_tiplab(size=3, color="purple")


select_clade <- groupClade(true_RTtree_outgroup, c(151, 155, 198))


roundrect_plot <- ggtree(select_clade, layout="roundrect", aes(color=group)) +
  theme(legend.position='none')+ 
  labs(fill = "amino acids letters" ) + 
  geom_tiplab(size= 2, color="black") +
  geom_point2(aes(subset=(node==22)), shape=21, size=2, fill='darkorange') +
  scale_color_manual(values=c("black", 
                              "seagreen3", 
                              "goldenrod1", 
                              "deeppink2"),
                     labels=c("Acca sellowiana athila/tat", "Ogre athila/tat db", "Athila/tat db", "Athila db"))+ 
  theme(legend.position = "none") 


msaplot(roundrect_plot, 
        fasta="/home/usuario/Data_Rstudio/tesina_g/tablas/domainLTR/aligment_abs_all/out_aligment_RT_abs.pep_form_MSA_Rstudio", 
        window=NULL,
        offset = 0,
        width = 1,
        color = NULL)
