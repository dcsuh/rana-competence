#Daniel Suh and Andrew Park
#6/3/20

#Figure 2 in rv_cc manuscript

library(tidyverse)
library(magrittr)
library(ggforce)
library(here)



site_scores %>% ggplot(.,aes(x=pc1Rank,y=cc))+
  geom_point(aes(color = Month.1, shape = factor(WetAltID)), show.legend = F)+
  theme_classic()+labs(x="Principal Component 1 Rank", y = "Community Competence") + 
  scale_shape_manual(values = rep(1:20, len = 20))+
  facet_zoom(xlim= c(5,25))

site_scores %>% ggplot(.,aes(x=pc1Rank,y=cc))+
  geom_point(aes(color = Month, shape = factor(WetAltID), size = Size))+
  labs(x="Principal Component 1 Rank", y = "Community Competence (CC)", size = "Size") + 
  scale_shape_manual(values = rep(1:20, len = 20)) +
  scale_size_continuous(range = c(2,10)) +
  guides(shape=F) +
  theme_classic() +
  theme(aspect.ratio=1/1.67, legend.position = c(0.85,0.4), legend.box = "horizontal") +
  scale_color_viridis_d(direction = -1)
  
#scale_color_manual(values = c( "#4575b4", "#91bfdb", "#e0f3f8", "#fee090", "#fc8d59", "#d73027"))


# #make same plot with NMDS scores rather than PCA scores
# library(vegan)
# set.seed(3000)
# 
# vegdist <- vegdist(community_mat, method = "bray")
# 
# #stress test for nmds
# stress <- c()
# for (i in 1:10){
#   output <- metaMDS(community_mat, distance = "bray", k = i, trace = F)
#   stress[i] <- output$stress
# }
# 
# stress %<>% as.data.frame(.) %>% rename(., stress = .) %>% mutate(., dim = as.numeric(rownames(.)))
# 
# ggplot(stress, aes (x = dim, y = stress)) + geom_point()
# #the difference between NMDS and PCA is that PCA uses euclidean distances while NMDS rank orders observations for ordination (this is why it is non-metric?)
# 
# NMDS1 <- metaMDS(community_mat, distance = "bray", k = 6, trace = F, autotransform = FALSE)
# stressplot(NMDS1)
# ordiplot(NMDS1, type = "n")
# orditorp(NMDS1, display = "species", col = "red")
# orditorp(NMDS1, display = "sites", cex = 1.1)
# 
# 
# points <- as.data.frame(NMDS1$points)
# points %<>% select(MDS1) %>% mutate(., siteID = rownames(.))
# points %<>% left_join(tmp, .)
# points %<>% mutate(rank = dense_rank(desc(MDS1)))
# 
# points %>% ggplot(.,aes(x=rank,y=cc))+
#   geom_point(aes(color = Month, shape = factor(WetAltID), size = Size))+
#   labs(x="NMDS Component 1 Rank", y = "Community Competence (CC)") + 
#   scale_shape_manual(values = rep(1:20, len = 20)) +
#   guides(shape=F) +
#   theme_classic() +
#   theme(aspect.ratio=1/1.67, legend.box = "horizontal") +
#   scale_color_brewer(palette = "Paired")
