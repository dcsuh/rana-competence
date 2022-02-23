#Script Originator: Daniel Suh

#This script creates supplemental figure 6
#Objective: A plot that models the cophenetic distance of a community and relative abundances of species
#This provides a look at patterns related to species abundance and phylogenetic relatedness
#We hypothesize that the species in our data are most abundant when they are moderately related to their neighbor
#If their closest neighbor is very phylogenetically close then competitive exclusion is likely to take place
#If their closest neighbor is very phylogenetically distant then environmental filtering may reduce coexistence
#If their closest neighbor is moderately phylogenetically close then they may coexist at high abundance

library(here)

source(here("base","src.R"))


comm_data <- readRDS(here("processed_data","comm_data.rds"))
pruned_tree <- readRDS(here("processed_data","pruned_tree.rds"))


plot(pruned_tree)
pd <- cophenetic(pruned_tree)
rownames(pd) <- gsub("_"," ",rownames(pd))
#mpd and min.pd

y <- comm_data
y %<>% pivot_longer(starts_with("AB"),names_to="species",values_to="abund")
names2 <- names
names2 %<>% mutate(species=paste("AB",species_code,sep=""))
names2 %<>% mutate(Species=str_replace_all(species_name,"_"," "))
names2 %<>% select(species,Species)
names2 %<>% mutate(Species=if_else(Species=="Anaxyrus terrestris","Bufo terrestris",Species))
names2 %<>% mutate(Species=if_else(Species=="Lithobates catesbeianus","Rana catesbeiana",Species))
names2 %<>% mutate(Species=if_else(Species=="Lithobates sphenocephalus","Rana sphenocephala",Species))
names2 %<>% mutate(Species=if_else(Species=="Lithobates clamitans","Rana clamitans",Species))
y %<>% left_join(.,names2)
# relabundo
y %<>% rowwise() %>% mutate(relAbund=abund/size)
y %<>% drop_na(Species)
y %<>% mutate(mpd=-999,min.pd=-999,mpdAbund=-999,min.pdAbund=-999,neighbor=-999)
#optional focus on top4
#top4sampledSyn <- c(top4sampled,"Bufo terrestris")
#y %<>% filter(Species %in% top4sampledSyn)
#end optional focus

#y is a df where every row is one species in one site-month
#so row numbers should equal number of distinct species * number of site-months
#i.e. y = 16 species * 96 site-months = 1536 rows


y %<>% filter(abund>0) #shouldn't we filter out any species that were absent? (DS 5/18/21)
#the mean pd and min pd should just be calculated for species that are present but without 
#filtering out those absent then you calculate mean pd and min pd for all of the possible 16 species




for (k in 1:dim(y)[1]){ #loop through each row
  #get mean pd and min pd for each row
  tmp <- y %>% filter(siteID==y$siteID[k]) #filter for single site-month
  pd.row <- which(rownames(pd)==y$Species[k]) #select rows from pd matching site-month
  #select species that are not the species from that row
  otherSpp <- tmp %>% filter(Species!=y$Species[k]) %>% select(Species) %>% pull() 
  pd.cols <- which(rownames(pd) %in% otherSpp) #select pd for not focal species
  myDists <- pd[pd.row,pd.cols] #get distances between focal species and other species
  mpd <- mean(myDists) #calculate mean pd
  min.pd <- min(myDists) #calculate min pd
  neighbor <- names(which(myDists==min.pd))
  if(is.null(neighbor)){
    neighbor <- otherSpp
  }
  if(identical(neighbor, character(0))){
    neighbor <- "no neighbor"
  }
  y$mpd[k] <- mpd #add to df
  y$min.pd[k] <- min.pd #add to df
  y$neighbor[k] <- neighbor
  #do same thing for relative abundances of each species in each site-month
  tmp2 <- tmp %>% filter(Species!=y$Species[k]) #make tmp df for species that are not focal species
  tmp2 %<>% filter(relAbund == max(relAbund)) #filter for max relative abundance
  pd.colsAbund <- which(rownames(pd) %in% tmp2$Species) 
  myDistsAbund <- pd[pd.row,pd.colsAbund]
  mpdAbund <- mean(myDistsAbund)
  min.pdAbund <- min(myDistsAbund)
  y$mpdAbund[k] <- mpdAbund
  y$min.pdAbund[k] <- min.pdAbund
}

#each point is one species in one site-month
#min.pd is the 
y %>% filter(relAbund>0) %>% ggplot(.,aes(x=mpd,y=relAbund))+
  geom_point()+geom_smooth(method="loess",span=1.4)
#above line gives us the relative abundance of each species and the mean pd for their site-month
#hosts in the same site-month should be stacked vertically
supp6 <- y %>% filter(!is.na(mpd)) %>% ggplot(.,aes(x=min.pd,y=relAbund))+
  geom_point()+geom_smooth(method="loess",span=1.4)+
  theme_classic()+labs(x="Distance to Closest Neighbor", y = "Relative Abundance")
supp6
#supp6 gives us the relative abundance for each species and the pd for their closest phylogenetic neighbor in their site-month
#so each vertical stack should be an instance where two species co-occurred
#the species included are limited to those that were sampled for ranavirus
#need to figure out why there are 9 NA's

y %>% filter(relAbund>0) %>% ggplot(.,aes(x=mpdAbund,y=relAbund))+
  geom_point()+geom_smooth(method="loess",span=1.4)
y %>% filter(relAbund>0) %>% ggplot(.,aes(x=min.pdAbund,y=relAbund))+
  geom_point()+geom_smooth(method="loess",span=1.4)



#ggsave("supp6.png",plot=supp6,device="png",path=here("figures"))
