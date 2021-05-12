library(ape)
tr2 <- ape::read.nexus("./data/asup_just_tree.txt")
mySpp <- c(union(names$tnrs_name,names$species_name),"Bufo terrestris")
mySpp <- gsub(" ","_",mySpp)
getRid <- setdiff(tr2$tip.label,mySpp)
tr2.pruned <- drop.tip(tr2,getRid)
plot(tr2.pruned)
pd <- cophenetic(tr2.pruned)
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
y %>% filter(relAbund>0) %>% ggplot(.,aes(x=min.pd,y=relAbund))+geom_point()+geom_smooth(method="loess",span=1.4)
y %>% filter(relAbund>0) %>% ggplot(.,aes(x=mpdAbund,y=relAbund))+geom_point()+geom_smooth(method="loess",span=1.4)
y %>% filter(relAbund>0) %>% ggplot(.,aes(x=min.pdAbund,y=relAbund))+geom_point()+geom_smooth(method="loess",span=1.4)