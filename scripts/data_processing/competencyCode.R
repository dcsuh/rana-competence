#Script originated by Andrew Park and edited by Daniel Suh

library(tidyverse)
library(magrittr)

#This script calculates the value for community competence at each sitexmonth (i.e. community)
#This assumes all individuals are exposed - so negative individuals count towards the average of a species

q <- read_csv("Ranavirus_Salamander_091819.csv")

v <- q %>% group_by(Species) %>% summarize(meanVL=mean(SQMean))

q %<>% mutate(vl.4=v$meanVL[v$Species==4])
q %<>% mutate(vl.9=v$meanVL[v$Species==9])
q %<>% mutate(vl.20=v$meanVL[v$Species==20])
q %<>% mutate(vl.21=v$meanVL[v$Species==21])
q %<>% mutate(vl.24=v$meanVL[v$Species==24])
q %<>% mutate(vl.26=v$meanVL[v$Species==26])
q %<>% mutate(vl.27=v$meanVL[v$Species==27])
q %<>% mutate(vl.28=v$meanVL[v$Species==28])
q %<>% mutate(vl.29=v$meanVL[v$Species==29])
q %<>% mutate(vl.31=v$meanVL[v$Species==31])
q %<>% mutate(vl.32=v$meanVL[v$Species==32])
q %<>% mutate(vl.34=v$meanVL[v$Species==34])
q %<>% mutate(vl.35=v$meanVL[v$Species==35])
q %<>% mutate(vl.38=v$meanVL[v$Species==38])
q %<>% mutate(vl.39=v$meanVL[v$Species==39])
q %<>% mutate(vl.41=v$meanVL[v$Species==41])
q %<>% mutate(vl.42=v$meanVL[v$Species==42])

q %<>% mutate(cc=RA4*vl.4+
                RA9*vl.9+
                RA20*vl.20+
                RA21*vl.21+
                RA24*vl.24+
                RA26*vl.26+
                RA27*vl.27+
                RA28*vl.28+
                RA29*vl.29+
                RA31*vl.31+
                RA34*vl.34+
                RA35*vl.35+
                RA38*vl.38+
                RA39*vl.39+
                RA41*vl.41+
                RA42*vl.42)

write_csv(q,path="competence.csv")
