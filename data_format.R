## ------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(vegan)
library(patchwork)
library(ggnewscale)
library(ape)
library(rotl)
library(picante)
library(here)


## ------------------------------------------------------------------------------------------------------------------------------------
here()
data <- read_csv(here("data/weighted_prev_competence_111220.csv"))
tree <- ape::read.nexus(here("data/asup_just_tree.txt"))
names <- read_csv(here("data/species_names_ids.csv"))


## ------------------------------------------------------------------------------------------------------------------------------------
comm_summ <- data %>% dplyr::select(WetAltID, Month.1, Month, cc, Month, AB2:AB42, Prevalence, MeanWaterTempPredC) %>% 
  dplyr::select(-ABAmb, -ABSal) %>%  
  mutate(., size = rowSums(.[5:25]), 
    ABHigh = (AB9 + AB21 + AB26 + AB42), 
    ABLow = (AB2 + AB3 + AB4 + AB5 + AB6 + AB8 + AB20 + AB24 + AB27 + AB28 + AB29 + AB31 + AB34 + AB35 + AB38 + AB39 + AB41)) %>% #get totals for community size
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>%
  distinct()
comm_summ$siteID <- gsub("Month", "", comm_summ$siteID)


## ------------------------------------------------------------------------------------------------------------------------------------
vl <- data %>% dplyr::select(vl.4:vl.42) %>% distinct() %>% pivot_longer(vl.4:vl.42)
vl %<>% mutate(species_code = as.double(gsub("vl.","",vl$name)))
vl %<>% full_join(., names, by="species_code")
vl %<>% mutate(name = gsub("_"," ",vl$species_name))
vl$value %<>% na_if(.,0)
vl %<>% mutate(., ln_value = log(value)) %>% mutate(., log10_value = log10(value))
vl %<>% replace_na(.,list(value=0,ln_value=0,log10_value=0))


## ------------------------------------------------------------------------------------------------------------------------------------
community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) 
community_mat$Month <- gsub("Month", "", community_mat$Month.1)
community_mat %<>% mutate(., siteID = paste(WetAltID, Month, sep = "_")) %>% distinct() %>% arrange(.,WetAltID) %>% dplyr::select(siteID, AB2:AB42)

#make abundance and pres/abs matrices with rownames as siteIDs
abundance_mat <- community_mat %>% column_to_rownames(., var = "siteID")
presence_mat <- as.data.frame(ifelse(abundance_mat>0, 1, 0))


## ------------------------------------------------------------------------------------------------------------------------------------
#select for environmental variables
env <- data %>% dplyr::select(WetAltID, Month.1, MeanAirT, MeanWaterTempPredC, DryingScore, CanopyCover, Area, Perimeter) %>% distinct() %>% arrange(.,WetAltID) 
env$Month <- gsub("Month", "", env$Month.1)
#make siteID rowname
env %<>% mutate(., siteID = paste(env$WetAltID, env$Month, sep = "_")) %>% distinct() %>% arrange(.,WetAltID) %>% dplyr::select(siteID, MeanAirT:Perimeter)

#need to remove duplicated rows and determine what to do with NA's
#env_mat <- env %>% column_to_rownames(., var = "siteID")


## ------------------------------------------------------------------------------------------------------------------------------------
#Get species richness from pres/abs matrix
richness <- tibble(siteID="", richness=1:nrow(presence_mat))
for (i in 1:nrow(presence_mat)){
  richness[i,1] <- rownames(presence_mat)[i]
  richness[i,2] <- sum(presence_mat[i,])
}


## ------------------------------------------------------------------------------------------------------------------------------------
#pca on community matrix and gather scores for first two principal components
pca <- princomp(abundance_mat, scores = TRUE)
summary(pca)
pca_scores <- as.data.frame(pca$scores)
pca_scores %<>% dplyr::select(Comp.1, Comp.2) %>% rownames_to_column(., var = "siteID")

site_scores  <- full_join(comm_summ, pca_scores, by = "siteID")

#order months
site_scores$Month <- factor(site_scores$Month, levels = c("Feb", "Mar", "Apr", "May", "Jun", "Jul"))

#rank order components and plot in order
site_scores %<>% mutate(pc1Rank=dense_rank(Comp.1))


## ------------------------------------------------------------------------------------------------------------------------------------
#select for environmental variables
env <- data %>% dplyr::select(WetAltID, Month.1, MeanAirT, MeanWaterTempPredC, DryingScore, CanopyCover, Area, Perimeter) %>% distinct() %>% arrange(.,WetAltID) 
env$Month.1 <- gsub("Month", "", env$Month.1)
#make siteID rowname
env$siteID <- paste(env$WetAltID, env$Month.1, sep = "_")
env %<>% dplyr::select(-WetAltID, -Month.1)
#remove duplicated rows (siteID is duplicated but env. variables are not... not sure why this is)
env <- env[-c(53, 90),]
env_mat <- env %>% column_to_rownames(., var = "siteID")


## ------------------------------------------------------------------------------------------------------------------------------------
richness_cc <- full_join(comm_summ, richness, by = "siteID")

#J = H'/ln(S)
#H' = Shannon diversity index
#S = Species Richness

shannon <- as_tibble(vegan::diversity(abundance_mat, index = "shannon"),rownames= "siteID")

shannon %<>% mutate(h = value) %>% dplyr::select(siteID, h)

evenness <- richness_cc %>% mutate(lnS = log(richness)) %>% rowwise()

evenness %<>% full_join(shannon, evenness, by = "siteID")

evenness %<>% mutate(J = h/lnS) %>% rowwise()
evenness <- evenness[order(evenness$WetAltID),]
#create new df for ordered sites to test lagged richness and evenness on prevalence
lag_evenness <- evenness[-c(4,14,47,90),] #brute force remove site_months out of sequence
lag_evenness$lag_richness <- c(0)
lag_evenness$lag_J <- c(0)
lag_evenness$lag_cc <- c(0)
lag_evenness$lag_size <- c(0)
for(n in 2:92){
  lag_evenness$lag_richness[n] <- lag_evenness$richness[n-1]
  lag_evenness$lag_J[n] <- lag_evenness$J[n-1]
  lag_evenness$lag_cc[n] <- lag_evenness$cc[n-1]
  lag_evenness$lag_size[n] <- lag_evenness$size[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
lag_evenness %<>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)


## ------------------------------------------------------------------------------------------------------------------------------------
#barplots comparing community composition
abundances <- comm_summ %>% full_join(.,evenness)


abundances$siteID <- reorder(abundances$siteID, -abundances$size) #order by total community size
#abundances$siteID <- reorder(abundances$siteID, -abundances$cc) #order by community competence
#abundances$siteID <- reorder(abundances$siteID, -abundances$J) #order by Pielou's J



## ------------------------------------------------------------------------------------------------------------------------------------
prev_cc <- comm_summ

order <- prev_cc[order(prev_cc$WetAltID, prev_cc$Month.1),]

#this removes months out of sequence
order <- order[-c(4,14,47,90),]

order %<>% add_column(lag_cc = NA, lag_size = NA, lag_temp = NA)


#create new column that includes previous month's value for cc
for(n in 2:nrow(order)){
  order$lag_cc[n] <- order$cc[n-1]
  order$lag_size[n] <- order$size[n-1]
  order$lag_temp[n] <- order$MeanWaterTempPredC[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
clean <- order %>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)


## ------------------------------------------------------------------------------------------------------------------------------------
rv_sampled <- data %>% select(Species) %>% distinct() #species codes for those sampled for rv
rv_sampled %<>% filter(Species!=32) #remove species 32 because there is no abundance data
keep <- c(rv_sampled$Species)

resolved <- tnrs_match_names(names$tnrs_name)
resolved %<>% rename("tnrs_name" = "unique_name")
#resolved$search_string <- paste(toupper(substr(resolved$search_string, 1, 1)), 
#                                  substr(resolved$search_string,2,nchar(resolved$search_string)),sep = "")

resolved %<>% select(tnrs_name, ott_id)

names %<>% left_join(.,resolved)

names %<>% filter(species_code %in% keep)

names %<>% filter(tnrs_name!="Dryophytes avivoca") #drop species because there is no abundance data
mySpp <- c(union(names$tnrs_name,names$species_name),"Bufo terrestris") #make vector for all species names
mySpp <- gsub(" ","_",mySpp)
getRid <- setdiff(tree$tip.label,mySpp)
pruned_tree <- drop.tip(tree,getRid)

