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

setwd("/home/usuario/acca/tree/gy_tree/")

##Then, we read tree____________________________________________________________

gy_df_tree <- read.tree("alg_acca_gy_CRM.bionj")

##Tidytree package provides as_tibble method to convert the phylo to a tidy data frame

gy_df <-as_tibble(gy_df_tree) #convert phylo to data frame
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

gy_outgroup <- root(gy_df_tree, outgroup = "RT_CRM")

##phylogenetics tree____________________________________________________________

p1 <- ggtree(gy_outgroup) +
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

#plot_list(p1, p2, p3, ncol=3, tag_levels='A') Creat sum of graph, label letter_

#Drop.tip we want to remove three tips (colored red) from the tree______________
'''
f <- system.file("extdata/NHX", "phyldog.nhx", package="treeio")
nhx <- read.nhx(f)
to_drop <- c("Physonect_sp_@2066767",
             "Lychnagalma_utricularia@2253871",
             "Kephyes_ovata@2606431")
p1 <- ggtree(nhx) + geom_tiplab(aes(color = label %in% to_drop)) +
  scale_color_manual(values=c("black", "red")) + xlim(0, 0.8)

nhx_reduced <- drop.tip(nhx, to_drop)
p2 <- ggtree(nhx_reduced) + geom_tiplab() + xlim(0, 0.8)  
plot_list(p1, p2, ncol=2, tag_levels = "A")
'''
##Subsetting tree by tip label__________________________________________________

#This isn't applicable to the outgroup, because it selects all__________________
##tree_subset()function will internally call the groupOTU()   

p1 <- ggtree(gy_outgroup) + 
  geom_tiplab(offset=.5) +  xlim(0, 10) + theme_tree2()

tree2 = tree_subset(gy_outgroup, 465)

p2 <- ggtree(tree2, aes(color=group)) +
  scale_color_manual(values = c("black", "red"), guide = 'none') +
  geom_tiplab(offset=.2) +  xlim(0, 4.5) + theme_tree2() 

tree3 = tree_subset(gy_outgroup,"RT_Grande1-4", levels_back = 6)#Date frame

p3 <- ggtree(tree3, aes(color=group)) +
  scale_color_manual(values = c("black", "red"), guide = 'none') +
  geom_tiplab(offset=.2) +  xlim(0, 4.5) + theme_tree2() 

tree4 = tree_subset(gy_outgroup,"tig00001711_arrow__1", levels_back = 8)#Date frame

p4 <- ggtree(tree4, aes(color=group)) +
  scale_color_manual(values = c("black", "red"), guide = 'none') +
  geom_tiplab(offset=.2) +  xlim(0, 4.5) + theme_tree2() 

p5 <- p4 + geom_point(aes(fill = rate), shape = 21, size = 4)

##______________________________________________________________________________

##Manipulating tree data for vizualization______________________________________

set.seed(1024)
x <- rtree(30)
y <- rtree(30)
p1 <- ggtree(x, layout='roundrect') + 
  geom_hilight(
    mapping=aes(subset = node %in% c(38, 48, 58, 36),
                node = node,
                fill = as.factor(node)
    )
  ) +
  labs(fill = "clades for tree in left" )

p2 <- ggtree(y)

d1 <- p1$data
d2 <- p2$data

## reverse x-axis and 
## set offset to make the tree on the right-hand side of the first tree
d2$x <- max(d2$x) - d2$x + max(d1$x) + 1

pp <- p1 + geom_tree(data=d2, layout='ellipse') +      
  ggnewscale::new_scale_fill() +
  geom_hilight(
    data = d2, 
    mapping = aes( 
      subset = node %in% c(38, 48, 58),
      node=node,
      fill=as.factor(node))
  ) +
  labs(fill = "clades for tree in right" ) 

dd <- bind_rows(d1, d2) %>% 
  filter(!is.na(label))

pp + geom_line(aes(x, y, group=label), data=dd, color='grey') +
  geom_tiplab(geom = 'shadowtext', bg.colour = alpha('firebrick', .5)) +
  geom_tiplab(data = d2, hjust=1, geom = 'shadowtext', 
              bg.colour = alpha('firebrick', .5))

##We are Manipulating our tree data for vizualization___________________________
#_______________________________________________________________________________

##asigno numero de nodos
ggtree(gy_outgroup) + geom_text(aes(label=node), hjust=-0.2) + 
  geom_tiplab(size=3, color="purple")

##marco las ramas de mi interes
p1 <- ggtree(gy_outgroup) +
  geom_hilight(
    mapping=aes(subset = node %in% c(465),
                node = node,
                fill = as.factor(node),
                size = 5
    )
  ) +
  labs(fill = "clades for tree in left" )+ 
  geom_tiplab(size= 2, color="purple")

##extract clade
clade_athila <- extract.clade(gy_outgroup, 465)

##visualizo con numero de nodos de nuevo
ggtree(clade_athila) + geom_text(aes(label=node), hjust=-0.2) + 
  geom_tiplab(size=3, color="purple")

##cladograma circular
cirk_plot <- ggtree(clade_athila, layout="circular") +
  geom_hilight(
    mapping=aes(subset = node %in% c(478, 465, 463),
                node = node,
                fill = as.factor(node),
                size = 5
    )
  ) +
  labs(fill = "clades for tree in left" )+ geom_tippoint(color="#b5e521", alpha=1/4, size=5)

round_plot <- ggtree(clade_athila, layout="roundrect") +
  geom_hilight(
    mapping=aes(subset = node %in% c(478, 465, 463),
                node = node,
                fill = as.factor(node)
    )
  ) +
  labs(fill = "clades for tree in left" )+ geom_tippoint(color="#574953", alpha=1/2, size=3)

round_plot2 <- ggtree(clade_athila, layout="roundrect") + 
  xlim(NA, 4) +
  geom_cladelab(node=481, label="Athila clade from db", align=TRUE, 
                geom='label', fill='lightblue') +
  geom_cladelab(node=c(465,463), label="Athila/tat clade from db", align=TRUE, 
                geom='label', fill='lightblack') +
  geom_tippoint(color="#574953", alpha=1/2, size=3)




##Displaying nodes/tips_________________________________________________________
## We must improve y-axis scale_________________________________________________

ggtree(gy_outgroup) + 
  geom_point(aes(shape=isTip, color=isTip), size=3)
p <- ggtree(gy_outgroup) + geom_tippoint(color="#b5e521", alpha=1/4, size=5)
p + geom_tiplab(size=3, color="purple")
p + geom_tiplab(as_ylab=TRUE, color='firebrick')

#chek this point________________________________________________________________

ggtree(gy_df_tree, continuous = 'colour', yscale = "trait") + 
  scale_color_viridis_c() + theme_minimal()
#Estos son pruebas______________________________________________________________

plot(gy_outgroup, edge.width = 1, label.offset = 2)
nodelabels()

ggtree(gy_outgroup) + geom_text(aes(label=node), hjust=-0.5)

group_athi <- extract.clade(gy_outgroup, 465)

ggtree(group_athi) + geom_tiplab()

##para esta parte debo tener solo mi clado de interes___________________________

Id_clade_athila <- clade_athila$tip.label

write.csv(Id_clade_athila, "/home/usuario/acca/alignment/Id_clade_athila.csv", row.names = FALSE)

p_alg_athtat <- ggtree(clade_athila) +
  geom_hilight(
    mapping=aes(subset = node %in% c(463, 465, 481),
                node = node,
                fill = as.factor(node)
    )
  ) +
  labs(fill = "amino acids letters" )+ 
  geom_tiplab(size= 2, 
              color="black")


msaplot(p_alg_athtat, 
        fasta="/home/usuario/acca/alignment/alg_rt_branch_athi", 
        window=c(1, 170))

msaplot(p=ggtree(clade_athila), 
        fasta="/home/usuario/acca/alignment/alg_rt_branch_athi", 
        window=c(1, 170))


msaplot(cirk_plot,
        fasta ="/home/usuario/acca/alignment/alg_rt_branch_athi", 
        window=c(120, 200))

k_datab <- groupClade(clade_athila, c(463, 465, 481))

ramas_col <- ggtree(k_datab, aes(color=group)) + 
  theme(legend.position='none')+ 
  labs(fill = "amino acids letters" ) + 
  geom_tiplab(size= 2, color="black") +
  scale_color_manual(values=c("black", 
                              "deeppink2", 
                              "goldenrod1", 
                              "seagreen3"))+ 
  theme(legend.position = "right") 

msaplot(ramas_col, 
        fasta="/home/usuario/acca/alignment/alg_rt_branch_athi", 
        window=c(1, 320))



cut_branch <- ggtree(k_datab, aes(color=group)) + 
  theme(legend.position='none')+ 
  labs(fill = "amino acids letters" ) + 
  geom_text(aes(label=node), hjust=-0.2) +
  geom_tiplab(size = 0.3, align = TRUE)+
  theme_tree2()+
  scale_color_manual(values=c("black", 
                              "deeppink2", 
                              "goldenrod1", 
                              "seagreen3"))+ 
  theme(legend.position = "right") 

msaplot(cut_branch, 
        fasta="/home/usuario/acca/alignment/alg_rt_branch_athi", 
        window=NULL,
        offset = 0,
        width = 1,
        color = NULL)
