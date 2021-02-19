#Script Originator: Daniel Suh

#This script creates figure 4


library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(patchwork)
library(ggnewscale)
library(here)

source(knitr::purl(here("/scripts/data_format.Rmd"), quiet=TRUE))


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
z %<>% select(siteID,ABLow,AB26,AB42,AB21,AB9,J,size)
z %<>% mutate(Low=ABLow)
z %<>% pivot_longer(cols=2:6,names_to="Species",values_to="abund") #every row is a site-month-species?
#z %<>% filter(J<0.6)
z$siteID <- reorder(z$siteID, -z$size) #order by total community size

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


#looking at correlations between total community abundance and community competence
evenness %>% ggplot(.,aes(x=cc,y=size)) +
  geom_point() +
  geom_smooth(method="lm")

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
