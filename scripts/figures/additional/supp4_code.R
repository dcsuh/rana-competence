## Script originator: Daniel Suh
## Date created: Mar. 29, 2021

#This script creates supplementary figure 3: PCA

library(here)

source(here("base","src.R"))


abundance_mat <- readRDS(here("processed_data","abund_mat.rds"))
env_mat <- readRDS(here("processed_data","env_mat.rds"))


community_pca <- princomp(abundance_mat, scores = TRUE)
summary(community_pca)

comm_scores <- as.data.frame(community_pca$loadings[,1:2])
env_vec <- envfit(community_pca, env_mat, na.rm = T) 
env_vec <- as.data.frame(scores(env_vec, display = "vectors"))

supp4 <- comm_scores %>% ggplot(.) +
  xlab("PC1 (71.2% Var. Explained)")+ylab("PC2 (26% Var. Explained)")+
  geom_segment(data=comm_scores,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "darkorange1") +
  geom_label_repel(data = comm_scores, aes(x = Comp.1, y = Comp.2, label = row.names(comm_scores)),
                   size = 2,color="darkorange1", segment.colour = 'black', max.overlaps = 35)+
  geom_segment(data=env_vec,
               aes(x = 0, xend = Comp.1, y = 0, yend = Comp.2),
               arrow = arrow(length = unit(0.1, "cm")), colour = "purple") +
  geom_label_repel(data = env_vec, aes(x = Comp.1, y = Comp.2, label = rownames(env_vec)),
                   size = 2,color="purple", segment.colour = 'black')+
  coord_fixed()

ggsave("supp4.png",plot=supp4,device="png",path=here("figures"))
