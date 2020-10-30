#Daniel Suh and Andrew Park
#6/3/20

#Figure 3 in rv_cc manuscript

library(tidyverse)
library(magrittr)
library(ecodist)
library(vegan)
library(here)

data <- read_csv(here("data/weighted_prev_competence.csv"))


#make community matrix
community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() %>% arrange(.,WetAltID) %>% dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, AB2:AB42)

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

nrow(community_mat) == nrow(env_mat)

site_labels <- rownames(community_mat)

#bray-curtis dissimilarity scores
spec_dist <- bcdist(community_mat)
spec_nmds <- nmds(spec_dist, mindim = 2, maxdim = 2)
spec_scores <- nmds.min(spec_nmds)

#mahalanobis dissimilarity scoring for sites by environmental variables
env_dist <- ecodist::distance(env_mat, method = "mahalanobis")
env_nmds <- nmds(env_dist, mindim = 2, maxdim = 2)
env_scores <- nmds.min(env_nmds)

#pca on community matrix and gather scores for first two principal components
pca <- princomp(community_mat, scores = TRUE)
scores <- as.data.frame(pca$scores)
scores %<>% dplyr::select(Comp.1, Comp.2) %>% rownames_to_column(., var = "siteID")



output=princomp(env_mat, cor=T) # PCA
env_scores <- as.data.frame(output$scores)
env_scores %<>% dplyr::select(Comp.1, Comp.2) %>% rownames_to_column(., var = "siteID")

biplot(output) # Generates a bi-plot with vectors

vec1=envfit(output$score[,1:2], community_mat, permutations=0) 
plot(vec1, col="blue")
scores12=output$score[,1:2] # get the PC1 and PC2 scores 
scores12_spec=cbind(scores12,community_mat) # merge with species 
correlations=cor(scores12_spec) # calculate correlations 
correlations12=correlations[3:15,1:2] # the loadings we want




site_scores <- data %>% dplyr::select(siteID, cc) %>% distinct() %>% rename_with(., ~ gsub("Month","",.x,fixed = "TRUE"))

scores12 %<>% as.data.frame(.)
scores12 %<>% mutate(siteID=rownames(.))
tmp <- site_scores %>% select(siteID,cc)
scores12 %<>% left_join(.,tmp)
scores12 %<>% mutate(CC=as.factor(ifelse(cc>45000,"Hi","Lo")))
spp.scrs <- as.data.frame(scores(vec1, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
vec2=envfit(output$score[,1:2], env_mat, permutations=0)
env.scrs <- as.data.frame(scores(vec2, display = "vectors"))
env.scrs <- cbind(env.scrs, Species = rownames(env.scrs))

library(ggrepel)
myScaleUp1 <- 5
myScaleUp2 <- 2.5
scores12 %>% drop_na(cc) %>% ggplot(.)+geom_point(aes(x=Comp.1,y=Comp.2,fill=CC,color=CC),pch=21)+
  xlab("PC1 (48.2%)")+ylab("PC2 (27.5%)")+
  scale_fill_manual(values=c("darkgreen","lightgreen"))+
  scale_color_manual(values=c("darkgreen","black"))+
  coord_fixed()+
  geom_segment(data=spp.scrs,
               aes(x = 0, xend = myScaleUp1*Comp.1, y = 0, yend = myScaleUp1*Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_text_repel(data = spp.scrs, aes(x = myScaleUp1*Comp.1, y = myScaleUp1*Comp.2, label = Species),
                  size = 2.5,color="darkorange2")+
  geom_segment(data=env.scrs,
               aes(x = 0, xend = myScaleUp2*Comp.1, y = 0, yend = myScaleUp2*Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_text_repel(data = env.scrs, aes(x = myScaleUp2*Comp.1, y = myScaleUp2*Comp.2, label = Species),
                  size = 2.5,color="purple")

#need to look at getting the dots in there and finding a way to make it legible while keeping the original vector lengths
#I need to be able to better explain what the direction and magnitude of each vector represents

scores12 %>% ggplot(.)+geom_point(aes(x=Comp.1,y=Comp.2),pch=21)+
  xlab("PC1 (48.2%)")+ylab("PC2 (27.5%)")+
  scale_fill_manual(values=c("darkgreen","lightgreen"))+
  scale_color_manual(values=c("darkgreen","black"))+
  coord_fixed()+
  geom_segment(data=spp.scrs,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_text_repel(data = spp.scrs, aes(x = Comp.1, y = Comp.2, label = Species),
                  size = 2.5,color="darkorange2")+
  geom_segment(data=env.scrs,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_text_repel(data = env.scrs, aes(x = Comp.1, y = Comp.2, label = Species),
                  size = 2.5,color="purple")
