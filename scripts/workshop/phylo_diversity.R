

library(here)


source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))



#subset <- names[unique(c(which(names$species_name %in% species),which(names$unique_name %in% species))),]
subset <- names %>% dplyr::select(species_name, species_code, tnrs_name)
subset <- subset[order(subset$species_code),]
#print(subset$species_code)

#subset %<>% filter(species_code!=c(31,32))
IDs <- c(subset$species_code)


phy_comm <- community_mat %>% column_to_rownames(var = "siteID")
phy_comm_1 <- community_mat %>% column_to_rownames(var = "siteID")
sum(phy_comm)
for (i in 1:nrow(subset)){
  col_name <- paste("AB",IDs[i],sep = "")
  phy_comm <- rename(phy_comm, !!subset$species_name[i] := col_name)
}
#same thing but presence-absence matrix instead of abundances
phy_comm_pa <- as.data.frame(ifelse(phy_comm>0,1,0))

#use community matrix and phylogeny to calculate faith's phylogenetic diversity
pd <- pd(phy_comm, tree, include.root = F)

pd %<>% rownames_to_column(., var="siteID") 
#pd %<>% mutate(PD = replace_na(pd$PD, 0)) #is this okay to do?

pd %<>% full_join(abundances,pd,by="siteID")


pd %>% ggplot(.,aes(x=PD,y=cc)) + geom_point() + geom_smooth(method="loess")
pd %>% ggplot(.,aes(x=PD,y=size)) + geom_point() + geom_smooth(method="loess")

pd %>% ggplot(.,aes(x=PD)) + geom_histogram(binwidth=50)


pd %>% filter(is.na(PD)==F) %>% dplyr::select(siteID, PD) %>% distinct() %>%
  ggplot(., aes(x=siteID, y=PD)) +
  geom_point() +
  labs(title = "Values of phylogenetic diversity for sites") +
  theme_minimal()


#calculate mean pairwise distances
pruned_tree$tip.label <- strip_ott_ids(pruned_tree$tip.label,remove_underscores=T)
dist_mat <- as.data.frame(cophenetic(pruned_tree))
#rename a few species for the community data
subset_1 <- subset
subset_1 %<>% mutate(tnrs_name=if_else(tnrs_name=="Anaxyrus terrestris","Bufo terrestris",tnrs_name))
subset_1 %<>% mutate(tnrs_name=if_else(tnrs_name=="Lithobates catesbeianus","Rana catesbeiana",tnrs_name))
subset_1 %<>% mutate(tnrs_name=if_else(tnrs_name=="Lithobates sphenocephalus","Rana sphenocephala",tnrs_name))
subset_1 %<>% mutate(tnrs_name=if_else(tnrs_name=="Lithobates clamitans","Rana clamitans",tnrs_name))
subset_1 %<>% mutate(tnrs_name=if_else(tnrs_name=="Dryophytes cinereus","Hyla cinerea",tnrs_name))
subset_1 %<>% mutate(tnrs_name=if_else(tnrs_name=="Dryophytes gratiosus","Hyla gratiosa",tnrs_name))
subset_1 %<>% mutate(tnrs_name=if_else(tnrs_name=="Dryophytes femoralis","Hyla femoralis",tnrs_name))
subset_1 %<>% mutate(tnrs_name=if_else(tnrs_name=="Dryophytes chrysoscelis","Hyla chrysoscelis",tnrs_name))

for (i in 1:nrow(subset_1)){
  col_name <- paste("AB",IDs[i],sep = "")
  phy_comm_1 <- rename(phy_comm_1, !!subset_1$tnrs_name[i] := col_name)
}

phy_comm_1 %<>% select(-c("AB8","AB2","AB3","AB5","AB6")) 
#remove because we don't know the identity of this species or species is not sampled for ranavirus

mean_pairwise_distances <- mpd(phy_comm_1,dist_mat)

site_IDs <- row.names(phy_comm_1) #make sure these rownames match properly
amphib_mpd <- tibble("siteID" = site_IDs, "mpd" = mean_pairwise_distances)

amphib_mpd <- inner_join(pd,amphib_mpd,by="siteID")

ggplot(amphib_mpd, aes(x=mpd)) + geom_histogram(bins="40") + theme_minimal()+
  labs(x="Mean Pairwise Distances", y = "Frequency")
