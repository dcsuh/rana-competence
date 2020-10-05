#Daniel Suh and Andrew Park
#6/3/20

#Figure 4 in rv_cc manuscript

library(tidyverse)
library(magrittr)
library(ggforce)



data <- read_csv("prev_test.csv")



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
  
  scale_color_manual(values = c( "#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"))

