#Script Originator: Daniel Suh

#This script creates figure 6


library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(patchwork)
library(ggnewscale)
library(here)

data <- read_csv(here("data/weighted_prev_competence_111220.csv")) #read data

#Get species richness from matrix

community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() #select data and make siteID column

community_mat <- community_mat[order(community_mat$WetAltID, community_mat$Month.1),] #order by wetland and month

community_mat %<>% arrange(.,WetAltID) %>% dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, AB2:AB42)

abundances <- community_mat #make abundances df

community_mat %<>% column_to_rownames(., var = "siteID") #make sitexspecies matrix

presence_mat <- ifelse(community_mat>0, 1, 0) #make pres-abs sitexspecies matrix

siteID <- rownames(presence_mat)

richness_mat <- tibble()
for (i in 1:length(siteID)){
  richness_mat[i,1] <- siteID[i]
  richness_mat[i,2] <- sum(presence_mat[i,])
}
richness_mat %<>% rename(siteID = 1, richness = 2) #make df with siteIDs and corresponding species richness

#Get species richness for each site and community competence values

tmp <- data %>% dplyr::select(WetAltID, Month.1, cc, Prevalence,AB2:AB8, AB9, AB20:AB42) %>% 
    distinct() %>% 
    mutate(., siteID = paste(WetAltID, Month.1, sep = "_"))
tmp %<>% mutate(.,size = rowSums(.[5:25])) %>% select(-c(AB2:AB42))
tmp <- tmp[order(tmp$WetAltID, tmp$Month.1),]
tmp %<>% dplyr::select(-Month.1)

richness_cc <- full_join(richness_mat, tmp, by = "siteID") #make df with siteIDs, richness, and cc



#J = H'/ln(S)
#H' = Shannon diversity index
#S = Species Richness

h <- vegan::diversity(community_mat, index = "shannon") #calculate shannons index for each site-month
h <- melt(as.matrix(h)) #melt into df
h %<>% rename(siteID = Var1, h = value) %>% dplyr::select(-Var2) #get df with shannons index for each site-month

evenness <- richness_cc %>% mutate(lnS = log(richness)) %>% rowwise() 
evenness %<>% full_join(h, evenness, by = "siteID")
evenness %<>% mutate(J = h/lnS) %>% rowwise() #calculate Pielou's J

#create new df for ordered sites to test lagged richness and evenness on prevalence
lag_evenness <- evenness[-c(4,14,47,90),] #remove site-months that are out of sequence
lag_evenness$lag_richness <- c(0)
lag_evenness$lag_J <- c(0)
lag_evenness$lag_cc <- c(0)
lag_evenness$lag_size <- c(0)
for(n in 2:91){ #lag values by one month
  lag_evenness$lag_richness[n] <- lag_evenness$richness[n-1]
  lag_evenness$lag_J[n] <- lag_evenness$J[n-1]
  lag_evenness$lag_cc[n] <- lag_evenness$cc[n-1]
  lag_evenness$lag_size[n] <- lag_evenness$size[n-1]
}
#remove the first entry for each wetland to remove the carryover from the last wetland
lag_evenness %<>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)
lag_evenness %<>% select(siteID, lag_richness:lag_size)


#barplots comparing community composition



evenness %<>% select(-WetAltID)
abundances %<>% left_join(.,evenness, by = "siteID")
abundances %<>% left_join(.,lag_evenness, by = "siteID")
abundances <- mutate(abundances, ABHigh = (AB9 + AB21 + AB26 + AB42)) 
abundances <- mutate(abundances, ABLow = (AB2 + AB3 + AB4 + AB5 + AB6 + AB8 + AB20 + AB24 + AB27 + AB28 + AB29 + AB31 + AB34 + AB35 + AB38 + AB39 + AB41))
abundances <- mutate(abundances, total = ABHigh + ABLow)
#final df for plots. includes data for abundance, prevalence, cc, diversity measure, lagged values


abundances$siteID <- reorder(abundances$siteID, -abundances$total) #order by total community size
#abundances$siteID <- reorder(abundances$siteID, -abundances$cc) #order by community competence
#abundances$siteID <- reorder(abundances$siteID, -abundances$J) #order by Pielou's J

cc_plot <- abundances %>% dplyr::select(siteID, cc) %>% distinct() %>%
  ggplot(., aes(y=cc, x=siteID)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 0.8) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("CC") + xlab("site-months")
#labs(title = "Values of cc for sites ordered by community size (descending)") +

evenness_plot <- abundances %>% dplyr::select(siteID, J, richness) %>% distinct() %>%
  ggplot(., aes(y=J, x=siteID, fill = richness)) +
  scale_fill_viridis_c("Richness", direction = -1) +
  geom_bar(position = position_dodge(), stat = "identity", width = 0.8) +
  theme_minimal() +
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("Pielou's J")
#labs(title = "Values of evenness for sites ordered by community size (descending)") +

z <- abundances
z %<>% select(siteID,ABLow,AB26,AB42,AB21,AB9,J,total)
z %<>% rename(Low=ABLow)
z %<>% pivot_longer(cols=2:6,names_to="Species",values_to="abund") #every row is a site-month-species?
#z %<>% filter(J<0.6)
z$siteID <- reorder(z$siteID, -z$total) #order by total community size

AB_plot <- z %>% ggplot(.,aes(x=siteID,y=abund))+
  geom_bar(stat="identity", width = 0.8)+
  ylab("Abundance")+
  theme_minimal()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

RA_plot<- z %>% ggplot(.,aes(x=siteID,y=abund,fill=Species))+
  geom_bar(position="fill",stat="identity", width = 0.8)+
  scale_fill_manual(values = c("#238443", "#78C679", "#C2E699", "#FFFFB2", "black"))+
  theme_minimal()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("Rel. Abun.")

#still need to clean this up but final figure should have a similar format to this
final <- AB_plot/evenness_plot/RA_plot/cc_plot

final

########################################################################

# pd is dataframe made in phylo_diversity.R
# abundances %<>% left_join(.,pd)
# 
# 
# abundances$siteID <- reorder(abundances$siteID, -abundances$total) #order by total community size
# 
# pd_plot <- abundances %>% filter(is.na(cc)==F) %>% dplyr::select(siteID, PD) %>% distinct() %>%
#   ggplot(., aes(y=PD, x=siteID)) +
#   geom_bar(position = position_dodge(), stat = "identity") +
#   labs(title = "Values of phylogenetic diversity for sites ordered by community size (descending)") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=90))+
#   theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
# pd_plot
# 
# final <- evenness_plot/pd_plot/cc_plot
# 
# final
# 
# ggplot(abundances, aes(x=PD, y=J)) + geom_point() + geom_smooth()
# ggplot(abundances, aes(x=J, y=PD)) + geom_point() + geom_smooth()
