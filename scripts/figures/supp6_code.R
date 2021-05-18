#Script Originator: Daniel Suh

#This script creates supplemental figure 6


library(here)

source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))


plot(pruned_tree)
pd <- cophenetic(pruned_tree)
rownames(pd) <- gsub("_"," ",rownames(pd))
#mpd and min.pd
y <- abundances
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
y %<>% mutate(mpd=-999,min.pd=-999,mpdAbund=-999,min.pdAbund=-999)
#optional focus on top4
#top4sampledSyn <- c(top4sampled,"Bufo terrestris")
#y %<>% filter(Species %in% top4sampledSyn)
#end optional focus

#y is a df where every row is one species in one site-month
#so row numbers should equal number of distinct species * number of site-months
#i.e. y = 16 species * 96 site-months = 1536 rows


for (k in 1:dim(y)[1]){
  tmp <- y %>% filter(siteID==y$siteID[k])
  pd.row <- which(rownames(pd)==y$Species[k])
  otherSpp <- tmp %>% filter(Species!=y$Species[k]) %>% select(Species) %>% pull()
  pd.cols <- which(rownames(pd) %in% otherSpp)
  myDists <- pd[pd.row,pd.cols]
  mpd <- mean(myDists)
  min.pd <- min(myDists)
  y$mpd[k] <- mpd
  y$min.pd[k] <- min.pd
  tmp2 <- tmp %>% filter(Species!=y$Species[k])
  tmp2 %<>% filter(relAbund == max(relAbund))
  pd.colsAbund <- which(rownames(pd) %in% tmp2$Species)
  myDistsAbund <- pd[pd.row,pd.colsAbund]
  mpdAbund <- mean(myDistsAbund)
  min.pdAbund <- min(myDistsAbund)
  y$mpdAbund[k] <- mpdAbund
  y$min.pdAbund[k] <- min.pdAbund
}

y %>% filter(relAbund>0) %>% ggplot(.,aes(x=mpd,y=relAbund))+geom_point()+geom_smooth(method="loess",span=1.4)
supp6 <- y %>% filter(relAbund>0) %>% ggplot(.,aes(x=min.pd,y=relAbund))+geom_point()+geom_smooth(method="loess",span=1.4)
y %>% filter(relAbund>0) %>% ggplot(.,aes(x=mpdAbund,y=relAbund))+geom_point()+geom_smooth(method="loess",span=1.4)
y %>% filter(relAbund>0) %>% ggplot(.,aes(x=min.pdAbund,y=relAbund))+geom_point()+geom_smooth(method="loess",span=1.4)

ggsave("supp6.png",plot=supp6,device="png",path=here("figures"))
