#Script Originator: Daniel Suh

#This script creates figure 4


library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(patchwork)
library(ggnewscale)
library(ggtree)
library(rotl)
library(here)

source(knitr::purl(here("/scripts/data_format.Rmd"), quiet=TRUE))

resolved <- tnrs_match_names(vl$tnrs_name)
tree <- tol_induced_subtree(resolved$ott_id)
tree$tip.label <- strip_ott_ids(tree$tip.label,remove_underscores=T)
vl %<>% select(tnrs_name,value) %>% mutate(.,label = tnrs_name) %>% select(-tnrs_name)
tree %<>% full_join(.,vl)

#horizontal tree
p1.h <- ggtree(tree,color="gray")+
  geom_tippoint(aes(color=value))+
  scale_colour_viridis_c("Viral load",direction = -1)+
  geom_tiplab(size=3,hjust=0,offset=0.2,fontface="italic")+
  xlim(0, 10)
#vertical tree
p1.v <- ggtree(tree,color="gray")+
  geom_tippoint(aes(color=value))+
  scale_colour_viridis_c("Viral load",direction = -1)+
  geom_tiplab(size=2.2,hjust=1,offset=-0.6,fontface="italic",angle=75,)+
  scale_x_reverse(limits=c(20,0))+
  coord_flip() 

abundances$siteID <- reorder(abundances$siteID, -abundances$cc)
cc_plot <- abundances %>% dplyr::select(siteID, cc) %>% distinct() %>%
  ggplot(., aes(y=cc, x=siteID)) +
  geom_bar(position = position_dodge(), stat = "identity", width = 0.8) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("CC") + xlab("site-months")
#labs(title = "Values of cc for sites ordered by community size (descending)") +

# evenness_plot <- abundances %>% dplyr::select(siteID, J, richness) %>% distinct() %>%
#   ggplot(., aes(y=J, x=siteID, fill = richness)) +
#   scale_fill_viridis_c("Richness", direction = -1) +
#   geom_bar(position = position_dodge(), stat = "identity", width = 0.8) +
#   theme_minimal() +
#   theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("Pielou's J")
#labs(title = "Values of evenness for sites ordered by community size (descending)") +

z <- abundances
z %<>% select(siteID,ABLow,AB26,AB42,AB21,AB9,J,size,cc)
z %<>% mutate(Low=ABLow)
z %<>% pivot_longer(cols=2:6,names_to="Species",values_to="abund") #every row is a site-month-species?
#z %<>% filter(J<0.6)
z$siteID <- reorder(z$siteID, -z$cc) #order by total community size

# AB_plot <- z %>% ggplot(.,aes(x=siteID,y=abund))+
#   geom_bar(stat="identity", width = 0.8)+
#   ylab("Abundance")+
#   theme_minimal()+
#   theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

RA_plot<- z %>% ggplot(.,aes(x=siteID,y=abund,fill=Species))+
  geom_bar(position="fill",stat="identity", width = 0.8)+
  scale_fill_manual(values = c("#238443", "#78C679", "#C2E699", "#FFFFB2", "black"))+
  theme_minimal()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  ylab("Rel. Abun.")


# final <- p1.h|(cc_plot/RA_plot)
final <- (cc_plot/RA_plot)/p1.v
final


#looking at correlations between total community abundance and community competence
# evenness %>% ggplot(.,aes(x=cc,y=size)) +
#   geom_point() +
#   geom_smooth(method="lm")

