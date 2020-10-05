#Daniel Suh and Andrew Park
#6/3/20

#Figure 2 in rv_cc manuscript

library(tidyverse)
library(magrittr)
library(metacom)


data <- read_csv("weighted_prev_competence.csv")

#make community matrix
community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() %>% arrange(.,WetAltID) %>% dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, AB2:AB42)

community_mat %<>% column_to_rownames(., var = "siteID")

Imagine(t(community_mat), 
        fill = FALSE, 
        order = TRUE, 
        xlab = "", 
        ylab = "",
        xline = 5,
        yline = 4,
        sitenames = ,
        speciesnames = ,
        binary=T)

ordered_matrix <- OrderMatrix(community_mat)
rownames(ordered_matrix)

results <- Metacommunity(community_mat, binary = T)
