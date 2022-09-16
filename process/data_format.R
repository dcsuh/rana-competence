library(here)

source(here("base","src.R"))

data <- read_csv(here("data/weighted_prev_competence_111220.csv"))
tree <- ape::read.nexus(here("data/raw_data/asup_just_tree.txt"))
names <- read_csv(here("data/raw_data/species_names_ids.csv"))


#makes community summary where each row is one wetland in one month (what we consider as unique communities)
#includes data on community competence, prevalence, abundance of each species, and abundances of aggregated high and low competent species
comm_summ <- data %>% dplyr::select(WetAltID, Month.1, Month, cc, Month, AB2:AB42, Prevalence, MeanWaterTempPredC) %>% 
  dplyr::select(-ABAmb, -ABSal) %>%  
  mutate(., size = rowSums(.[5:25]), 
    ABHigh = (AB9 + AB21 + AB26 + AB42), 
    ABLow = (AB2 + AB3 + AB4 + AB5 + AB6 + AB8 + AB20 + AB24 + AB27 + AB28 + AB29 + AB31 + AB34 + AB35 + AB38 + AB39 + AB41)) %>% #get totals for community size
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>%
  distinct()
comm_summ$siteID <- gsub("Month", "", comm_summ$siteID)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# makes community matrix where each cell is the abundance of one species at one unique wetland-month combination
community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) 
community_mat$Month <- gsub("Month", "", community_mat$Month.1)
community_mat %<>% mutate(., siteID = paste(WetAltID, Month, sep = "_")) %>% distinct() %>% arrange(.,WetAltID) %>% dplyr::select(siteID, AB2:AB42)

#make abundance and pres/abs matrices with rownames as siteIDs
abundance_mat <- community_mat %>% column_to_rownames(., var = "siteID")
presence_mat <- as.data.frame(ifelse(abundance_mat>0, 1, 0))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
#same as community summary data but with added columns for PC1 and PC2 scores and PC1 rank
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

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#Get species richness from pres/abs matrix
#each row is one unique wetland-month combination and the host species richness
richness <- tibble(siteID="", richness=1:nrow(presence_mat))
for (i in 1:nrow(presence_mat)){
  richness[i,1] <- rownames(presence_mat)[i]
  richness[i,2] <- sum(presence_mat[i,])
}


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#same as community summary data with data on richness and evenness with and without lags included
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

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#make final community data
comm_data <- comm_summ %>% left_join(., evenness)
comm_data %<>% left_join(., site_scores)
comm_data$Month <- factor(comm_data$Month, levels = c("Feb", "Mar", "Apr", "May", "Jun", "Jul"))


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#includes average viral load data for each species
#each row is one host species
# vl <- data %>% dplyr::select(vl.4:vl.42) %>% distinct() %>% pivot_longer(vl.4:vl.42)
# vl %<>% mutate(species_code = as.double(gsub("vl.","",vl$name)))
# vl %<>% full_join(., names, by="species_code")
# vl %<>% mutate(name = gsub("_"," ",vl$species_name))
# vl$value %<>% na_if(.,0)
# vl %<>% mutate(., ln_value = log(value)) %>% mutate(., log10_value = log10(value))
# vl %<>% replace_na(.,list(value=0,ln_value=0,log10_value=0))

vl <- data %>% group_by(Species) %>% summarize(mean = mean(SQMean),
                                               var = var(SQMean),
                                               se = sqrt(var(SQMean)/n()),
                                               n = n())
vl_log <- data %>% 
  group_by(Species) %>% 
  summarize(log10_mean = mean(log10_SQ), 
            log10_var = var(log10_SQ),
            log10_se = sqrt(var(log10_SQ)/n())) %>% 
  mutate(species_code = Species)
vl %<>% left_join(.,vl_log)
vl %<>% full_join(., names, by="species_code")



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#select for environmental variables
env <- data %>% dplyr::select(WetAltID, Month.1, MeanAirT, MeanWaterTempPredC, DryingScore, CanopyCover, Area, Perimeter) %>% 
  distinct() %>% arrange(.,WetAltID) 
env$Month <- gsub("Month", "", env$Month.1)
#make siteID rowname
env %<>% mutate(., siteID = paste(env$WetAltID, env$Month, sep = "_")) %>% 
  distinct() %>% arrange(.,WetAltID) %>% dplyr::select(siteID, MeanAirT:Perimeter)

#duplicated rows because canopy cover is slightly different. This line keeps only the first instance
#13_1 and 19_4 were duplicated and their percent canopy cover were within three percent of their duplicate
env %<>% filter(duplicated(env$siteID)==F)

env_mat <- env %>% column_to_rownames(., var = "siteID")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------
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


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------

# please process data before this line 
if(dir.exists(here("processed_data")) == FALSE) {
  message("Welcome! Let's make some room for the processed data.")
  dir.create(here("processed_data")) 
} else {
  message("/processed_data exists! Proceeeding to save.")
}


saveRDS(comm_data, file = here("processed_data","comm_data.rds"))
saveRDS(vl, file = here("processed_data","vl.rds"))
saveRDS(pruned_tree, file = here("processed_data","tree.rds"))
saveRDS(community_mat, file = here("processed_data", "comm_mat.rds"))
saveRDS(abundance_mat, file = here("processed_data", "abund_mat.rds"))
saveRDS(env_mat, file = here("processed_data", "env_mat.rds"))
saveRDS(pruned_tree, file = here("processed_data", "pruned_tree.rds"))




