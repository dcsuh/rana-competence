#Script Originator: Daniel Suh

#This script creates figure 7


library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(patchwork)
library(ggnewscale)
library(here)

data <- read_csv(here("data/weighted_prev_competence_111220.csv"))

#Get species richness from matrix

community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() 

community_mat <- community_mat[order(community_mat$WetAltID, community_mat$Month.1),]

community_mat %<>% arrange(.,WetAltID) %>% dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, AB2:AB42)

abundances <- community_mat

community_mat %<>% column_to_rownames(., var = "siteID")

presence_mat <- ifelse(community_mat>0, 1, 0)

siteID <- rownames(presence_mat)

richness_mat <- tibble()
for (i in 1:length(siteID)){
  richness_mat[i,1] <- siteID[i]
  richness_mat[i,2] <- sum(presence_mat[i,])
}
richness_mat %<>% rename(siteID = 1, richness = 2)

order <- richness_mat

#Get species richness for each site and community competence values

tmp <- data %>% dplyr::select(WetAltID, Month.1, cc, Prevalence,AB2:AB8, AB9, AB20:AB42) %>% 
  distinct() %>% 
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_"))
tmp %<>% mutate(.,size = rowSums(.[5:25])) %>% select(-c(AB2:AB42))
tmp <- tmp[order(tmp$WetAltID, tmp$Month.1),]
tmp %<>% dplyr::select(-Month.1)

richness_cc <- full_join(tmp, richness_mat, by = "siteID")



#J = H'/ln(S)
#H' = Shannon diversity index
#S = Species Richness

h <- vegan::diversity(community_mat, index = "shannon")
h <- melt(as.matrix(h))
h %<>% rename(siteID = Var1, h = value) %>% dplyr::select(-Var2)

evenness <- richness_cc %>% mutate(lnS = log(richness)) %>% rowwise()

evenness %<>% full_join(h, evenness, by = "siteID")

evenness %<>% mutate(J = h/lnS) %>% rowwise()

#create new df for ordered sites to test lagged richness and evenness on prevalence
lag_evenness <- evenness[-c(4,14,47,90),]
lag_evenness$lag_richness <- c(0)
lag_evenness$lag_J <- c(0)
lag_evenness$lag_cc <- c(0)
lag_evenness$lag_size <- c(0)
for(n in 2:91){
  lag_evenness$lag_richness[n] <- lag_evenness$richness[n-1]
  lag_evenness$lag_J[n] <- lag_evenness$J[n-1]
  lag_evenness$lag_cc[n] <- lag_evenness$cc[n-1]
  lag_evenness$lag_size[n] <- lag_evenness$size[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
lag_evenness %<>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)

#fit two lines to cc_evenness_prev plot
lag_evenness$high <- c(0)
lag_evenness$low <- c(0)
for (i in 1:nrow(lag_evenness)){
  if (lag_evenness$lag_J[i] > 0.6){
    lag_evenness$high[i] <- 1
    lag_evenness$low[i] <- 1
  }
  else if (lag_evenness$lag_cc[i] < 45000){
    lag_evenness$low[i] <- 1
  }
  else if (lag_evenness$lag_cc[i] > 45000){
    lag_evenness$high[i] <- 1
  }
  else {
    lag_evenness$high[i] <- 0
    lag_evenness$low[i] <- 0
  }
}
low_cc <- lag_evenness %>% filter(low==1)
cor.test(x=low_cc$lag_J, y=low_cc$lag_cc,method = "spearman")
high_cc <- lag_evenness %>% filter(high==1)
cor.test(x=high_cc$lag_J, y=high_cc$lag_cc,method = "spearman")

p1 <- evenness %>% remove_missing() %>% ggplot(aes(x=richness, y=cc)) + 
  geom_smooth(method ="lm") +
  geom_point() +
  labs(x = "Richness", y = "CC", title = ) +
  theme_classic()

cor.test(richness_cc$cc, richness_cc$richness, method="spearman")

p2 <- lag_evenness %>% remove_missing() %>% ggplot(aes(x=lag_richness, y=Prevalence)) + 
  geom_smooth(method ="lm") +
  geom_point() +
  labs(x = "Richness", y = "Prevalence", title = ) +
  theme_classic()

cor.test(x = lag_evenness$lag_richness, y = lag_evenness$Prevalence, method = "spearman")

p3 <- lag_evenness %>% ggplot(aes(x = lag_J, y = lag_cc, color = Prevalence, size = lag_size)) +
  geom_point() +
  labs(x="Evenness", y = "CC", size = "Size") +
  scale_color_gradient(low="blue", high="red") + 
  theme_minimal()

mySpan=1.5
myColor="gray80"
p4 <- ggplot() +
  geom_point(data = lag_evenness, aes(x=lag_J, y=lag_cc, color = Prevalence, size = lag_size)) +
  scale_color_viridis_c("Prevalence\n(t+1)", direction = -1) +
  new_scale_color() +
  geom_smooth(data = low_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  new_scale_color() +
  geom_smooth(data = high_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  theme_minimal()+ylim(0,3e+05)+
  labs(x="Evenness", y = "CC", size = "Size") +
  annotate("text", label = "i", family="Times", fontface="italic", x = 0.134, y = 46000, size = 4, colour = "white")+
  annotate("text", label = "ii", family="Times", fontface="italic", x = 0.427, y = 61000, size = 3, colour = "white")
p4

final <- p4
