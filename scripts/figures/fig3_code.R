#Daniel Suh and Andrew Park
#6/3/20

#Figure 3 in rv_cc manuscript

library(tidyverse)
library(magrittr)
library(ecodist)
library(vegan)
library(here)

data <- read_csv(here("data/weighted_prev_competence_111220.csv"))


#make community matrix
community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42)
community_mat$Month.1 <- gsub("Month", "", community_mat$Month.1) 
community_mat %<>% mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% 
  distinct() %>% arrange(.,WetAltID) %>% dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, AB2:AB42)

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

site_labels <- rownames(community_mat)


pca_env_output=princomp(env_mat, cor=T) # PCA
biplot(pca_env_output) # Generates a bi-plot with vectors

spec_vec=envfit(pca_env_output$score[,1:2], community_mat, permutations=0) 
plot(spec_vec, col="blue")

#get loadings for environmental variables
env_var <- as.data.frame(pca_env_output$loadings[,1:2])


env_scores <- pca_env_output$score[,1:2] # get the PC1 and PC2 scores 
env_comm_scores <- cbind(env_scores,community_mat) # merge with species 
correlations <- cor(env_comm_scores)
spec_corr <- correlations[3:23,1:2] # the loadings we want



site_scores <- data %>% dplyr::select(WetAltID, Month.1, cc) 
site_scores$Month.1 <- gsub("Month", "", site_scores$Month.1)
site_scores %<>% mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% 
  distinct() %>% 
  arrange(.,WetAltID) %>% 
  dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, cc)


env_scores %<>% as.data.frame(.)
env_scores %<>% mutate(siteID=rownames(.))
tmp <- site_scores %>% select(siteID,cc)
env_scores %<>% left_join(.,tmp)
env_scores %<>% mutate(CC=as.factor(ifelse(cc>45000,"Hi","Lo")))

spec_scores <- as.data.frame(scores(spec_vec, display = "vectors"))


#what is the difference between spec_scores and spec_corr

env_scores %>% drop_na(cc) %>% ggplot(.) +
  geom_point(aes(x=Comp.1, y=Comp.2, shape = CC, color = cc)) +
  xlab("PC1 (48.2%)")+ylab("PC2 (27.5%)")+
  geom_segment(data=spec_scores,
             aes(x = 0, xend = myScaleUp1*Comp.1, y = 0, yend = myScaleUp1*Comp.2),
             arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_text_repel(data = spec_scores, aes(x = myScaleUp1*Comp.1, y = myScaleUp1*Comp.2, label = row.names(spec_scores)),
                  size = 2.5,color="darkorange2")+
  geom_segment(data=env_var,
               aes(x = 0, xend = myScaleUp2*Comp.1, y = 0, yend = myScaleUp2*Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_text_repel(data = env_var, aes(x = myScaleUp2*Comp.1, y = myScaleUp2*Comp.2, label = rownames(env_var)),
                  size = 2.5,color="purple")+
  coord_fixed()

env_scores %>% drop_na(cc) %>% ggplot(.) +
  xlab("PC1 (48.2%)")+ylab("PC2 (27.5%)")+
  geom_segment(data=spec_scores,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_text_repel(data = spec_scores, aes(x = Comp.1, y = Comp.2, label = row.names(spec_scores)),
                  size = 2.5,color="darkorange2")+
  geom_segment(data=env_var,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_text_repel(data = env_var, aes(x = Comp.1, y = Comp.2, label = rownames(env_var)),
                  size = 2.5,color="purple")+
  coord_fixed()

  
  
  
vec1=envfit(pca_env_output$score[,1:2], community_mat, permutations=0) 
spp.scrs <- as.data.frame(scores(vec1, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
vec2=envfit(pca_env_output$score[,1:2], env_mat, permutations=0)
env.scrs <- as.data.frame(scores(vec2, display = "vectors"))
env.scrs <- cbind(env.scrs, Species = rownames(env.scrs))

library(ggrepel)
myScaleUp1 <- 1
myScaleUp2 <- 1
env_scores %>% drop_na(cc) %>% ggplot(.)+
#  geom_point(aes(x=Comp.1,y=Comp.2,fill=CC,color=CC),pch=21)+
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

# scores12 %>% ggplot(.)+geom_point(aes(x=Comp.1,y=Comp.2),pch=21)+
#   xlab("PC1 (48.2%)")+ylab("PC2 (27.5%)")+
#   scale_fill_manual(values=c("darkgreen","lightgreen"))+
#   scale_color_manual(values=c("darkgreen","black"))+
#   coord_fixed()+
#   geom_segment(data=spp.scrs,
#                aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
#                arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
#   geom_text_repel(data = spp.scrs, aes(x = Comp.1, y = Comp.2, label = Species),
#                   size = 2.5,color="darkorange2")+
#   geom_segment(data=env.scrs,
#                aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
#                arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
#   geom_text_repel(data = env.scrs, aes(x = Comp.1, y = Comp.2, label = Species),
#                   size = 2.5,color="purple")



#try this again for indirect gradient analysis with pca
community_pca <- princomp(community_mat, scores = TRUE)
summary(community_pca)
biplot(community_pca)
comm_scores <- community_pca$scores[,1:2]
comm_env_scores <- cbind(comm_scores,env_mat)
corrs <- cor(comm_env_scores)
env_corrs <- corrs[3:8, 1:2]
env_corrs %<>% as.data.frame

comm_scores <- as.data.frame(community_pca$loadings[,1:2])
env_vec <- envfit(community_pca$score[,1:2], env_mat, permutations=0) 
env_vec <- as.data.frame(scores(env_vec, display = "vectors"))


comm_scores %>% ggplot(.) +
  xlab("PC1")+ylab("PC2")+
  geom_segment(data=comm_scores,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_text_repel(data = comm_scores, aes(x = Comp.1, y = Comp.2, label = row.names(comm_scores)),
                  size = 2.5,color="darkorange2")+
  geom_segment(data=env_vec,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_text_repel(data = env_vec, aes(x = Comp.1, y = Comp.2, label = rownames(env_vec)),
                  size = 2.5,color="purple")+
  coord_fixed()

comm_scores %>% ggplot(.) +
  xlab("PC1")+ylab("PC2")+
  geom_segment(data=comm_scores,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_segment(data=env_vec,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  coord_fixed()


comm_scores %<>% as.data.frame(.)
comm_scores %<>% mutate(siteID=rownames(.))
comm_scores %<>% left_join(.,tmp)

