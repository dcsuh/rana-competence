
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
library(rotl)
library(here)

data <- read_csv(here("data/weighted_prev_competence.csv"))

tree <- read.nexus("data/asup_just_tree.txt")

names <- read_csv("data/species_names_ids.csv")
resolved <- tnrs_match_names(names$species_name)
resolved$unique_name <- gsub(" ", "_",resolved$unique_name)
resolved$search_string <- paste(toupper(substr(resolved$search_string, 1, 1)), 
                                  substr(resolved$search_string,2,nchar(resolved$search_string)),sep = "")

resolved %<>% rename("species_name" = "search_string") %>% select(species_name, unique_name, ott_id)

names %<>% left_join(.,resolved)

species <- unique(c(names$species_name[which(names$species_name %in% tree$tip.label)], 
                    names$unique_name[which(names$unique_name %in% tree$tip.label)]))

tree %<>% drop.tip(.,tree$tip.label[-match(species, tree$tip.label)])

tree %>% ggtree() + geom_tiplab()

subset <- names[unique(c(which(names$species_name %in% species),which(names$unique_name %in% species))),]

subset <- subset[order(subset$species_code),]
print(subset$species_code)

subset %<>% filter(species_code!=c(31,32))
IDs <- c(subset$species_code)
community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB4:AB42) %>% dplyr::select(-c(ABAmb, ABSal, AB31, AB5, AB6, AB8, AB21)) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() 

no_phy_data <- data %>% dplyr::select(AB5, AB6, AB8, AB21, AB31) %>% distinct()
sum(no_phy_data)

phy_comm <- community_mat %>% dplyr::select(AB4:siteID) %>% column_to_rownames(var = "siteID")
sum(phy_comm)
for (i in 1:nrow(subset)){
  col_name <- paste("AB",IDs[i],sep = "")
  phy_comm <- rename(phy_comm, !!subset$species_name[i] := col_name)
}

pd <- pd(phy_comm, tree, include.root = F)

pd %<>% rownames_to_column(., var="siteID") 
pd %<>% mutate(PD = replace_na(pd$PD, 0)) #is this okay to do?

