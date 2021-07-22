#Script Originator: Daniel Suh

#This script creates figure 4


library(ggnewscale)
library(ggtree)
library(rotl)
library(here)

source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))
vl %<>% filter(tnrs_name!="Dryophytes avivoca") # as doesn't appear sampled in any communities? (AP)
#H. avivoca (ID=32) has no abundance measured but has one recording for viral load
#viral load for H. avivoca is actually the highest but should not be included because there is no abundance data


pruned_tree$tip.label <- gsub("_"," ",pruned_tree$tip.label)


vl %<>% select(species_name,value,ln_value,log10_value) %>% mutate(.,label = species_name) %>% select(-species_name)
vl$label <- gsub("_", " ", vl$label)
vl$label <- gsub("Anaxyrus terrestris","Bufo terrestris",vl$label)
vl$label <- gsub("Lithobates","Rana",vl$label)
vl$label <- gsub("catesbeianus","catesbeiana",vl$label)
vl$label <- gsub("sphenocephalus","sphenocephala",vl$label)
vl$label <- gsub("clamitans","clamitans",vl$label)


vl %<>% mutate(labelAlt=label)
top4sampled <- c("Pseudacris crucifer","Rana sphenocephala","Bufo terrestris","Notophthalmus viridescens")
vl %<>% mutate(labelAlt=as_factor(if_else(labelAlt %in% top4sampled,labelAlt,"Low competence")))
vl %<>% mutate(labelAlt=fct_relevel(labelAlt,c(top4sampled,"Low competence")))

vl_tree <- pruned_tree %>% full_join(.,vl)

abundances$siteID <- reorder(abundances$siteID, -abundances$cc)
z <- abundances
z %<>% select(siteID,ABLow,AB26,AB42,AB21,AB9,J,size,cc)
z %<>% mutate(Low=ABLow)
z %<>% pivot_longer(cols=2:6,names_to="Species",values_to="abund") #every row is a site-month-species?

z$siteID <- reorder(z$siteID, -z$cc) #order by community competence
z %<>% mutate(Species=if_else(Species=="AB21","Bufo terrestris",Species))
z %<>% mutate(Species=if_else(Species=="AB26","Pseudacris crucifer",Species))
z %<>% mutate(Species=if_else(Species=="AB42","Rana sphenocephala",Species))
z %<>% mutate(Species=if_else(Species=="AB9","Notophthalmus viridescens",Species))
z %<>% mutate(Species=if_else(Species=="ABLow","Other (low competence)",Species))
#z$Species <- factor(z$Species,levels=c(top4sampled,"Other (low competence)"))
z %<>% mutate(Species=fct_relevel(Species,c(top4sampled,"Other (low competence)")))
####### plot ########
cc_plot <- abundances %>% dplyr::select(siteID, cc) %>% distinct() %>%
  ggplot(., aes(y=cc, x=siteID)) +
  geom_bar(color="gray",fill="gray",position = position_dodge(), stat = "identity", width = 0.8) + 
  theme_minimal() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab("Community\ncompetence (CC)") + xlab("Communities")
RA_plot<- z %>% ggplot(.,aes(x=siteID,y=abund,fill=Species))+
  geom_bar(position="fill",stat="identity", width = 0.8)+
  scale_fill_manual(values = c("#238443", "#78C679", "#C2E699", "#FFFFB2", "gray"))+
  theme_minimal()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  ylab("Relative\nAbundance")
# newTree <- tree %>% ggtree(.,color="gray")+
#   geom_tippoint(aes(color=labelAlt,size=log10_value))+
#   scale_color_manual(name="Species",values = c("#238443", "#78C679", "#C2E699", "#FFFFB2", "gray"))+
#   geom_tiplab(size=2.2,hjust=,offset=7,angle=90)+#fontface="italic",
#   scale_size_continuous(name="log10(Viral load)")+
#   coord_flip()+guides(color=F)+xlim(0,300)
newTree <- vl_tree %>% ggtree(.,color="gray")+
  geom_tippoint(aes(color=labelAlt,size=log10_value))+
  scale_color_manual(name="Species",values = c("#238443", "#78C679", "#C2E699", "#FFFFB2", "gray"))+
  geom_tiplab(size=2.2,hjust=,offset=7,angle=0)+#fontface="italic",
  scale_size_continuous(name="log10(Viral load)")+
  guides(color=F)+xlim(0,350)


figure_4 <- cc_plot/RA_plot/newTree
figure_4

ggsave("fig_4.png",plot=figure_4,device="png",path=here("figures"))

