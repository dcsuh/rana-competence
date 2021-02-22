library(tidyverse)
library(magrittr)
library(ape)
library(rotl)
library(phangorn)
library(geiger)
library(phytools)
library(patchwork)
library(ggtree)
library(here)

#files
#prev_test.csv
#asup_just_tree.txt
source(knitr::purl(here("/scripts/data_format.Rmd"), quiet=TRUE))

names <- c("Notophthalamus viridescens", "Scaphiopus holbrookii", "Anaxyrus terrestris", "Pseudacris crucifer", "Hyla gratiosa", "Hyla avivoca",
           "Pseudacris nigrita", "Pseudacris ornata", "Gastrophyrne carolinensis", "Lithobates clamitans", "Lithobates sphenocephalus",
           "Acris gryllus", "Hyla chrysoscelis", "Hyla cinerea", "Hyla femoralis", "Lithobates catesbeianus", "Siren intermedia")
resolved <- tnrs_match_names(names)
tree_1 <- tol_induced_subtree(resolved$ott_id)
tree_1$tip.label <- strip_ott_ids(tree_1$tip.label,remove_underscores=T)
vl %<>% select(tnrs_name,value) %>% mutate(.,label = tnrs_name) %>% select(-tnrs_name)
tree_1 %<>% full_join(.,vl)
#horizontal tree
p1.h <- ggtree(tree_1,color="gray")+
  geom_tippoint(aes(color=value))+
  scale_colour_viridis_c("Viral load",direction = -1)+
  geom_tiplab(size=3,hjust=0,offset=0.2,fontface="italic")+
  xlim(0, 10)
#vertical tree
p1.v <- ggtree(tree_1,color="gray")+
  geom_tippoint(aes(color=value))+
  scale_colour_viridis_c("Viral load")+
  geom_tiplab(size=2.2,hjust=1,offset=-0.6,fontface="italic",angle=90)+
  scale_x_reverse(limits=c(20,0))+
  coord_flip() 




#trait1 <- rTraitCont(tree, "OU", theta = 0, alpha=.1, sigma=.01)
data <- read_csv(here("data/weighted_prev_competence_111220.csv"))

#Plot of RV status over time - ridges with count of ind. positive for RV at each sitextime
rv <- data %>% dplyr::select(Month,WetAltID,RV.Status)
rv$Month <- factor(data$Month,levels=c("Feb","Mar","Apr","May","Jun","Jul"))

#looking at rv of species and drawing a line to divide "low" and "high" competent
rv_species <- data %>% dplyr::select(vl.4:vl.42) %>% distinct()
rv_species %<>% pivot_longer(.,vl.4:vl.42) %>% filter(value>0)
rv_species %<>% add_column(log_value = log10(rv_species$value))
rv_species$name <- factor(rv_species$name, levels = rv_species$name[order(rv_species$log_value)])
#adjusted names manually to match the resolved names by tnrs_match_names
rv_species$actual_name <- factor(c("Notophthalmus viridescens", "Scaphiopus holbrookii", "Bufo terrestris", 
                                   "Acris gryllus", "Pseudacris crucifer", "Hyla gratiosa", "Hyla avivoca", 
                                   "Pseudacris nigrita", "Pseudacris ornata", "Gastrophryne carolinensis", 
                                   "Rana clamitans", "Rana sphenocephala"))
rv_species$Species <- factor(c("Notophthalmus viridescens", "Scaphiopus holbrookii", "Anaxyrus terrestris", 
                                   "Acris gryllus", "Pseudacris crucifer", "Hyla gratiosa", "Hyla avivoca", 
                                   "Pseudacris nigrita", "Pseudacris ornata", "Gastrophryne carolinensis", 
                                   "Lithobates clamitans", "Lithobates sphenocephalus"))

rv_species$actual_name <- factor(rv_species$actual_name, levels = rv_species$actual_name[order(rv_species$log_value)])
rv_species$actual_name <- as.character(rv_species$actual_name)
rv_species %<>% mutate(name_value = paste(actual_name, value, sep =" "))

names <- c("Notophthalamus viridescens", "Scaphiopus holbrookii", "Anaxyrus terrestris", "Pseudacris crucifer", "Hyla gratiosa", "Hyla avivoca",
           "Pseudacris nigrita", "Pseudacris ornata", "Gastrophyrne carolinensis", "Lithobates clamitans", "Lithobates sphenocephalus",
           "Acris gryllus", "Hyla chrysoscelis", "Hyla cinerea", "Hyla femoralis", "Lithobates catesbeianus", "Siren intermedia")


resolved <- tnrs_match_names(rv_species$actual_name)

tree <- tol_induced_subtree(resolved$ott_id)
#plot(tree, edge.width = "3")

plot.phylo(tree, edge.width = '3', show.tip.label = TRUE)

species <- resolved$unique_name
species <- gsub(" ", "_", species)

amphibians <- read.nexus("data/asup_just_tree.txt")
#plot(amphibians, edge.width = "1")


setdiff(species,amphibians$tip.label)
species <- gsub("Dryophytes","Hyla",species)
species <- gsub("gratiosus","gratiosa",species)
species <- gsub("cinereus","cinerea",species)
species <- gsub("Anaxyrus","Bufo",species)
setdiff(species,amphibians$tip.label)

#prune
pruned.tree<-drop.tip(amphibians,amphibians$tip.label[-match(species, amphibians$tip.label)])
pruned.tree.vl<-drop.tip(amphibians,amphibians$tip.label[-match(species, amphibians$tip.label)])


for (i in 1:17){
  for (j in 1:12)
  if(pruned.tree.vl$tip.label[i]==rv_species$actual_name[j]){
    pruned.tree.vl$tip.label[i] <- rv_species$name_value[j]
  }
}

plot(pruned.tree.vl)

#phylo distances
phy_dist <- cophenetic.phylo(pruned.tree)

library(phylobase)
ott_in_tree <- resolved$ott_id
tr <- tol_induced_subtree(ott_ids = ott_in_tree)
plot(tr)
pruned.tree$tip.label <- str_replace_all(pruned.tree$tip.label, "_", " ")
pruned.tree$tip.label

#using phylo4d from phylobase
values <- rv_species %>% select(value, log_value) #get values
rownames(values) <- rv_species$actual_name
tree_data <- phylo4d(pruned.tree, values)
plot(tree_data)

values$label <- rv_species$actual_name
tree_vl <- full_join(pruned.tree, values, by='label')

ggtree(tree_vl) +geom_tiplab(hjust=0,offset =5)+
  geom_tippoint(aes(color=value,size=10)) +
  scale_color_gradient(low = "#fdc70c",high = "#e93e3a") +
  theme(legend.position = c(0.10,0.75))


names <- c("Notophthalamus viridescens", "Scaphiopus holbrookii", "Anaxyrus terrestris", "Pseudacris crucifer", "Hyla gratiosa", "Hyla avivoca",
           "Pseudacris nigrita", "Pseudacris ornata", "Gastrophyrne carolinensis", "Lithobates clamitans", "Lithobates sphenocephalus",
           "Acris gryllus", "Hyla chrysoscelis", "Hyla cinerea", "Hyla femoralis", "Lithobates catesbeianus", "Siren intermedia")
resolved <- tnrs_match_names(names)
tree <- tol_induced_subtree(resolved$ott_id)
tree$tip.label <- strip_ott_ids(tree$tip.label,remove_underscores=T)
set.seed(123)
vlExample <- tibble(label=tree$tip.label,vl=runif(length(names)))
tree %<>% full_join(.,vlExample)
#horizontal tree
p1.h <- ggtree(tree,color="gray")+
  geom_tippoint(aes(color=vl))+
  scale_colour_viridis_c("Viral load")+
  geom_tiplab(size=3,hjust=0,offset=0.2,fontface="italic")+
  xlim(0, 10)
#vertical tree
p1.v <- ggtree(tree,color="gray")+
  geom_tippoint(aes(color=vl))+
  scale_colour_viridis_c("Viral load")+
  geom_tiplab(size=2.2,hjust=1,offset=-0.6,fontface="italic",angle=90)+
  scale_x_reverse(limits=c(20,0))+
  coord_flip() 
data(cars)
p2 <- cars %>% ggplot(.,aes(x=speed,y=dist))+geom_point()
p3 <- cars %>% ggplot(.,aes(x=speed,y=dist))+geom_bar(stat="identity")
p1.h|(p2/p3)
p1.v/(p2/p3)






#trait data

if (!require("traitdataform")) install.packages("traitdataform") # pkg cont'g AmphiBIO
require("traitdataform")
library(magrittr)
library(tidyverse)
traitdataform::pulldata("amphibio")
mySpp <- c("Rana clamitans","Rana sphenocephala","Gastrophryne carolinensis","Pseudacris nigrita","Pseudacris ornata","Pseudacris crucifer","Acris gryllus","Hyla avivoca","Hyla gratiosa","Bufo terrestris","Scaphiopus holbrookii","Notophthalmus viridescens","Anaxyrus terrestris","Lithobates clamitans","Lithobates sphenocephalus")
a <- amphibio
a %<>% dplyr::filter(Species %in% mySpp)
setdiff(mySpp,a$Species)
a %<>% select_if(~ !any(is.na(.)))


vl_traits <- full_join(rv_species, a, by="Species")

vl_traits %>% ggplot(aes(x=Body_size_mm, y=log_value)) +
  geom_point() + geom_smooth(method="lm")
cor.test(vl_traits$Body_size_mm, vl_traits$log_value, method = "spearman")

vl_traits %>% ggplot(aes(y=log_value, x=Reproductive_output_y)) +
  geom_point() + geom_smooth(method="lm")

vl_traits %>% ggplot(aes(y=log_value, x=Litter_size_min_n)) +
  geom_point() + geom_smooth(method="lm")

vl_traits %>% ggplot(aes(y=log_value, x=Offspring_size_min_mm)) +
  geom_point() + geom_smooth(method="lm")

vl_traits %>% select(log_value, Body_size_mm, Litter_size_min_n, Litter_size_max_n, Offspring_size_min_mm, Offspring_size_max_mm) %>%
  pairs(.)

vl_traits %>% select(log_value, Body_size_mm, Litter_size_min_n, Litter_size_max_n, Offspring_size_min_mm, Offspring_size_max_mm) %>%
  cor(.)
