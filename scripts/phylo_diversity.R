
library(tidyverse)
library(magrittr)
library(phytools)
library(phylobase)
library(treeio)
library(phanghorn)
library(ape)
library(ggtree)
library(tidytree)
library(picante)

data <- read_csv(here("data/weighted_prev_competence.csv"))

tree <- treeio::read.nexus("data/asup_just_tree.txt")

species <- c("Notophthalmus_viridescens", "Scaphiopus_holbrookii", "Bufo_terrestris", 
             "Acris_gryllus", "Pseudacris_crucifer", "Hyla_gratiosa", "Hyla_avivoca", 
             "Pseudacris_nigrita", "Pseudacris_ornata", "Gastrophryne_carolinensis", 
             "Rana_clamitans", "Rana_sphenocephala")

tree %<>% drop.tip(.,tree$tip.label[-match(species, tree$tip.label)])

tree %>% ggtree() + geom_tiplab()

community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() 

community_mat <- community_mat[order(community_mat$WetAltID, community_mat$Month.1),]

community_mat %<>% arrange(.,WetAltID) %>% dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, AB2:AB42)

