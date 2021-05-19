

library(here)


source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))



#subset <- names[unique(c(which(names$species_name %in% species),which(names$unique_name %in% species))),]
subset <- names %>% dplyr::select(species_name, species_code, tnrs_name)
subset <- subset[order(subset$species_code),]
#print(subset$species_code)

#subset %<>% filter(species_code!=c(31,32))
IDs <- c(subset$species_code)


phy_comm <- community_mat %>% column_to_rownames(var = "siteID")
sum(phy_comm)
for (i in 1:nrow(subset)){
  col_name <- paste("AB",IDs[i],sep = "")
  phy_comm <- rename(phy_comm, !!subset$species_name[i] := col_name)
}
#same thing but presence-absence matrix instead of abundances
phy_comm_pa <- as.data.frame(ifelse(phy_comm>0,1,0))

#use community matrix and phylogeny to calculate pairwise distances
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

