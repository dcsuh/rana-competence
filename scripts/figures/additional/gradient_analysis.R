#Daniel Suh and Andrew Park
#6/3/20

#Figure 3 in rv_cc manuscript

library(tidyverse)
library(magrittr)
library(ecodist)
library(vegan)
library(here)
library(ggrepel)


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

#make presence-absence matrix
presence_mat <- as.data.frame(ifelse(community_mat>0, 1, 0))

site_labels <- rownames(community_mat)


pca_env_output=princomp(env_mat, cor=T) # PCA
summary(pca_env_output)
biplot(pca_env_output) # Generates a bi-plot with vectors

spec_vec=envfit(pca_env_output$score[,1:2], community_mat) 
plot(spec_vec, col="blue")

#get loadings for environmental variables
env_var <- as.data.frame(pca_env_output$loadings[,1:2])


site_scores <- data %>% dplyr::select(WetAltID, Month.1, cc) 
site_scores$Month.1 <- gsub("Month", "", site_scores$Month.1)
site_scores %<>% mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% 
  distinct() %>% 
  arrange(.,WetAltID) %>% 
  dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, cc)

env_scores <- pca_env_output$score[,1:2] # get the PC1 and PC2 scores 
env_scores %<>% as.data.frame(.)
env_scores %<>% mutate(siteID=rownames(.))
tmp <- site_scores %>% select(siteID,cc)
env_scores %<>% left_join(.,tmp)
env_scores %<>% mutate(CC=as.factor(ifelse(cc>45000,"Hi","Lo")))

spec_scores <- as.data.frame(scores(spec_vec, display = "vectors"))


env_scores %>% drop_na(cc) %>% ggplot(.) +
  geom_point(aes(x=Comp.1, y=Comp.2, shape = CC, color = cc)) +
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



#try this again for indirect gradient analysis with pca
community_pca <- princomp(community_mat, scores = TRUE)
summary(community_pca)

comm_scores <- as.data.frame(community_pca$loadings[,1:2])
env_vec <- envfit(community_pca, env_mat, na.rm = T) 
env_vec <- as.data.frame(scores(env_vec, display = "vectors"))

comm_scores %>% ggplot(.) +
  xlab("PC1")+ylab("PC2")+
  geom_segment(data=comm_scores,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_segment(data=env_vec,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  coord_fixed()

comm_scores %>% ggplot(.) +
  xlab("PC1 (71.2% Var. Explained)")+ylab("PC2 (26% Var. Explained)")+
  geom_segment(data=comm_scores,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_label_repel(data = comm_scores, aes(x = Comp.1, y = Comp.2, label = row.names(comm_scores)),
               size = 2,color="darkorange1", segment.colour = 'black')+
  geom_segment(data=env_vec,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_label_repel(data = env_vec, aes(x = Comp.1, y = Comp.2, label = rownames(env_vec)),
                  size = 2,color="purple", segment.colour = 'black')+
  coord_fixed()

comm_scores %>% ggplot(.) +
  xlab("PC1")+ylab("PC2")+
  geom_segment(data=env_vec,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_label_repel(data = env_vec, aes(x = Comp.1, y = Comp.2, label = rownames(env_vec)),
                   size = 2.5,color="purple", segment.colour = 'black')+
  coord_fixed()


#indirect gradient analysis with pca for presence-absence matrix
presence_pca <- princomp(presence_mat, scores = TRUE)
summary(presence_pca)

pres_scores <- as.data.frame(presence_pca$loadings[,1:2])
env_vec <- envfit(presence_pca, env_mat, na.rm = T) 
env_vec <- as.data.frame(scores(env_vec, display = "vectors"))

pres_scores %>% ggplot(.) +
  xlab("PC1 (21.1%)")+ylab("PC2 (12.6%)")+
  geom_segment(data=pres_scores,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_label_repel(data = pres_scores, aes(x = Comp.1, y = Comp.2, label = row.names(comm_scores)),
                   size = 2,color="darkorange1", segment.colour = 'black')+
  geom_segment(data=env_vec,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_label_repel(data = env_vec, aes(x = Comp.1, y = Comp.2, label = rownames(env_vec)),
                   size = 2,color="purple", segment.colour = 'black')+
  coord_fixed()
