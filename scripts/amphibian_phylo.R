
library(here)

#make phylogeny with viral load
#make phylogney with host traits
#calculate phylogenetic distance for each community


#files

#asup_just_tree.txt
source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))


resolved <- tnrs_match_names(vl$tnrs_name)
tree <- tol_induced_subtree(resolved$ott_id)
tree$tip.label <- strip_ott_ids(tree$tip.label,remove_underscores=T)
vl %<>% select(tnrs_name,value) %>% mutate(.,label = tnrs_name) %>% select(-tnrs_name)
tree %<>% full_join(.,vl)

#horizontal tree
p1.h <- ggtree(tree,color="gray")+
  geom_tippoint(aes(color=value))+
  scale_colour_viridis_c("Viral load",direction = -1)+
  geom_tiplab(size=3,hjust=0,offset=0.2,fontface="italic")+
  xlim(0, 10)
#vertical tree
p1.v <- ggtree(tree_1,color="gray")+
  geom_tippoint(aes(color=value))+
  scale_colour_viridis_c("Viral load",direction = -1)+
  geom_tiplab(size=2.2,hjust=1,offset=-0.6,fontface="italic",angle=90)+
  scale_x_reverse(limits=c(20,0))+
  coord_flip() 



#trait data

if (!require("traitdataform")) install.packages("traitdataform") # pkg cont'g AmphiBIO
require("traitdataform")

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

