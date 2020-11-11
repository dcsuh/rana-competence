#Daniel Suh and Andrew Park
#6/3/20

#Figure 4 in rv_cc manuscript

library(tidyverse)
library(magrittr)
library(ggforce)
library(here)


data <- read_csv("data/weighted_prev_competence.csv")




#select for community variables
community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) 
community_mat$Month.1 <- gsub("Month", "", community_mat$Month.1)
community_mat %<>% mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() %>% arrange(.,WetAltID) %>% dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, AB2:AB42)
#rownames
community_mat %<>% column_to_rownames(., var = "siteID")

#select for environmental variables
env <- data %>% dplyr::select(WetAltID, Month.1, MeanAirT, MeanWaterTempPredC, DryingScore, CanopyCover, Area, Perimeter) %>% distinct() %>% arrange(.,WetAltID) 
env$Month.1 <- gsub("Month", "", env$Month.1)
#make siteID rowname
env$siteID <- paste(env$WetAltID, env$Month.1, sep = "_")
env %<>% dplyr::select(-WetAltID, -Month.1)
#remove duplicated rows (siteID is duplicated but env. variables are not... not sure why this is)
env <- env[-c(53, 90),]
env_mat <- env %>% column_to_rownames(., var = "siteID")

#replace NA's with zeros. Probably not the best option but will allow analysis to run for now
community_mat[is.na(community_mat)] <- 0
env_mat[is.na(env_mat)] <- 0

#pca on community matrix and gather scores for first two principal components
pca <- princomp(community_mat, scores = TRUE)
scores <- as.data.frame(pca$scores)
scores %<>% dplyr::select(Comp.1, Comp.2) %>% rownames_to_column(., var = "siteID")

tmp <- data %>% dplyr::select(WetAltID, Month.1, cc, Month, AB2:AB42) %>% 
  dplyr::select(-ABAmb, -ABSal) %>%  mutate(., Size = rowSums(.[5:25])) %>% #get totals for community size
  dplyr::select(-AB2:-AB42) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>%
  distinct()
tmp$siteID <- gsub("Month", "", tmp$siteID)

site_scores  <- inner_join(tmp, scores, by = "siteID")

#order months
site_scores$Month <- factor(site_scores$Month, levels = c("Feb", "Mar", "Apr", "May", "Jun", "Jul"))

#rank order components and plot in order
site_scores %<>% mutate(pc1Rank=dense_rank(Comp.1))
site_scores %>% ggplot(.,aes(x=pc1Rank,y=cc))+
  geom_point(aes(color = Month.1, shape = factor(WetAltID)), show.legend = F)+
  theme_classic()+labs(x="Principal Component 1 Rank", y = "Community Competence") + 
  scale_shape_manual(values = rep(1:20, len = 20))+
  facet_zoom(xlim= c(5,25))

site_scores %>% ggplot(.,aes(x=pc1Rank,y=cc))+
  geom_point(aes(color = Month, shape = factor(WetAltID), size = Size))+
  labs(x="Principal Component 1 Rank", y = "Community Competence (CC)") + 
  scale_shape_manual(values = rep(1:20, len = 20)) +
  guides(shape=F) +
  theme_classic() +
  theme(aspect.ratio=1/1.67, legend.position = c(0.85,0.4), legend.box = "horizontal") +
  scale_color_brewer(palette = "Paired")
  
#scale_color_manual(values = c( "#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"))


#make same plot with NMDS scores rather than PCA scores
library(vegan)
set.seed(3000)

vegdist <- vegdist(community_mat, method = "bray")

#stress test for nmds
stress <- c()
for (i in 1:10){
  output <- metaMDS(community_mat, distance = "bray", k = i, trace = F)
  stress[i] <- output$stress
}

stress %<>% as.data.frame(.) %>% rename(., stress = .) %>% mutate(., dim = as.numeric(rownames(.)))

ggplot(stress, aes (x = dim, y = stress)) + geom_point()
#the difference between NMDS and PCA is that PCA uses euclidean distances while NMDS rank orders observations for ordination (this is why it is non-metric?)

NMDS1 <- metaMDS(community_mat, distance = "bray", k = 4, trymax = 100, trace = F)
stressplot(NMDS1)
ordiplot(NMDS1, type = "n")
orditorp(NMDS1, display = "species", col = "red")
orditorp(NMDS1, display = "sites", cex = 1.1)


points <- as.data.frame(NMDS1$points)
points %<>% select(MDS1) %>% mutate(., siteID = rownames(.))
points %<>% left_join(tmp, .)
points %<>% mutate(rank = dense_rank(MDS1))

points %>% ggplot(.,aes(x=rank,y=cc))+
  geom_point(aes(color = Month, shape = factor(WetAltID), size = Size))+
  labs(x="NMDS Component 1 Rank", y = "Community Competence (CC)") + 
  scale_shape_manual(values = rep(1:20, len = 20)) +
  guides(shape=F) +
  theme_classic() +
  theme(aspect.ratio=1/1.67, legend.box = "horizontal") +
  scale_color_brewer(palette = "Paired")
