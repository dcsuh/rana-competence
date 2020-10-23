#Script Originator: Daniel Suh

#This script creates figure 6


library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(patchwork)
library(ggnewscale)
library(here)

data <- read_csv(here("data/weighted_prev_competence.csv"))

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

richness_cc <- full_join(tmp, richness_mat, by = "siteID") %>% remove_missing()



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


#barplots comparing community composition

competence_values <- data %>% select(Species, SQMean) %>% distinct() %>% group_by(Species) %>% summarize(comp = mean(SQMean))

log_abundances <- abundances %>% as_tibble() %>% column_to_rownames(., var = "siteID") %>% na_if(.,0) %>% log10(.) %>% 
  replace(is.na(.), 0) %>% rownames_to_column(., var = "siteID") %>% rename_with(., ~ gsub("AB","logAB",.x,fixed = "TRUE"))
log_abundances %<>% mutate(., logABLow = (logAB2 + logAB3 + logAB4 + logAB5 + logAB6 + logAB8 + 
                                            logAB20 + logAB24 + logAB27 + logAB28 + logAB29 + 
                                            logAB31 + logAB34 + logAB35 + logAB38 + logAB39 + logAB41))
log_abundances %<>% mutate(., logABhigh = (logAB9 + logAB21 + logAB26 + logAB42))
log_abundances %<>% mutate(., logtotal = rowSums(.[2:22]))

abundances <- full_join(abundances,lag_evenness, by = "siteID")
abundances <- mutate(abundances, ABHigh = (AB9 + AB21 + AB26 + AB42)) 
abundances <- mutate(abundances, ABLow = (AB2 + AB3 + AB4 + AB5 + AB6 + AB8 + AB20 + AB24 + AB27 + AB28 + AB29 + AB31 + AB34 + AB35 + AB38 + AB39 + AB41))
abundances <- mutate(abundances, total = ABHigh + ABLow)
log_abundances <- full_join(abundances,log_abundances, by = "siteID")


abundances$siteID <- reorder(abundances$siteID, -abundances$total) #order by total community size
log_abundances$siteID <- reorder(log_abundances$siteID, -log_abundances$total) #order by total community size
#abundances$siteID <- reorder(abundances$siteID, -abundances$cc) #order by community competence
#abundances$siteID <- reorder(abundances$siteID, -abundances$J) #order by Pielou's J
comp_plot_a <- abundances %>% remove_missing() %>% filter(J<0.6)  %>% dplyr::select(siteID, AB26, AB42, AB21, AB9, ABLow) %>% melt() %>%
  ggplot(., aes(y = log10(value), x = siteID, fill = variable)) +
  geom_bar(position ="stack", stat = "identity") +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#ffffb2", "black")) +
  labs(title = "Abundance of high and low competent species under Pielou's J of 0.6") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Abundance") +
  labs(fill = "log10(variable)")

comp_plot_b <- abundances %>% remove_missing() %>% filter(J<0.6)  %>% dplyr::select(siteID, ABHigh, ABLow) %>% melt() %>%
  ggplot(., aes(y = log10(value), x = siteID, fill = variable)) +
  geom_bar(position ="stack", stat = "identity") +
  scale_fill_manual(values = c("#238443","black")) +
  labs(title = "Abundance of high and low competent species under Pielou's J of 0.6") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Abundance") +
  labs(fill = "log10(variable)")

comp_plot_c <- log_abundances %>% remove_missing() %>% filter(J<0.6)  %>% dplyr::select(siteID, logAB26, logAB42, logAB21, logAB9, logABLow) %>% melt() %>%
  ggplot(., aes(y = (value), x = siteID, fill = variable)) +
  geom_bar(position ="stack", stat = "identity") +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#ffffb2", "black")) +
  labs(title = "Abundance of high and low competent species under Pielou's J of 0.6") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Abundance")

cc_plot <- abundances %>% remove_missing() %>% dplyr::select(siteID, cc) %>% distinct() %>%
  ggplot(., aes(y=cc, x=siteID)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  labs(title = "Values of cc for sites ordered by community size (descending)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))


evenness_plot <- abundances %>% dplyr::select(siteID, J) %>% distinct() %>%
  ggplot(., aes(y=J, x=siteID)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  labs(title = "Values of evenness for sites ordered by community size (descending)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

  
comp_plot_a / comp_plot_b / comp_plot_c

comp_plot_c / cc_plot


z <- abundances
z %<>% select(siteID,ABLow,AB26,AB42,AB21,AB9,J,total)
z %<>% rename(Low=ABLow)
z %<>% pivot_longer(cols=2:6,names_to="Species",values_to="abund")
#z %<>% filter(J<0.6)
z$siteID <- reorder(z$siteID, -z$total) #order by total community size
p1top <- z %>% ggplot(.,aes(x=siteID,y=abund))+
  geom_bar(stat="identity")+
  ylab("Abundance")+theme_minimal()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p1bot<- z %>% ggplot(.,aes(x=siteID,y=abund,fill=Species))+
  geom_bar(position="fill",stat="identity")+
  scale_fill_manual(values = c("#238443", "#78C679", "#C2E699", "#FFFFB2", "black"))+
  theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("Relative abundance")

#still need to clean this up but final figure should have a similar format to this
final <- p1top/evenness_plot/p1bot/cc_plot

final

########################################################################

#pd is dataframe made in phylo_diversity.R
abundances %<>% left_join(.,pd)


abundances$siteID <- reorder(abundances$siteID, -abundances$total) #order by total community size

pd_plot <- abundances %>% dplyr::select(siteID, PD) %>% distinct() %>%
  ggplot(., aes(y=PD, x=siteID)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  labs(title = "Values of phylogenetic diversity for sites ordered by community size (descending)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
pd_plot

final <- p1top/evenness_plot/p1bot/pd_plot/cc_plot

final
