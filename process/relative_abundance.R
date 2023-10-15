##Calculating Relative Abundance of amphibian species at each site
##Script Originator: Daniel Suh
##Date: 11/12/20


#Process this script first


#load packages
library(tidyverse)
library(magrittr)
library(here)

data <- read_csv(here("data/raw_data/Ranavirus_Salamander_091819.csv"))
data %<>% drop_na(ID)
data %<>% select(X1:AB42) #remove original RA values
data %<>% mutate(., total = AB9 + AB20 + AB21 + 
                   AB24 +  AB26 + AB28 + 
                   AB29 + AB34 + AB35 + 
                   AB38 + AB39 + AB41 + AB42) #sum is everything except for ABAmb and ABSal
#also remove 4, 27, 31 because low sampling for ranavirus
data %<>% mutate(.,
#                 RA4 = AB4/total,
                   RA9 = AB9/total,
                   RA20 = AB20/total,
                   RA21 = AB21/total,
                   RA24 = AB24/total,
                   RA26 = AB26/total,
#                   RA27 = AB27/total,
                   RA28 = AB28/total,
                   RA29 = AB29/total,
 #                  RA31 = AB31/total,
                   RA34 = AB34/total,
                   RA35 = AB35/total,
                   RA38 = AB38/total,
                   RA39 = AB39/total,
                   RA41 = AB41/total,
                   RA42 = AB42/total)

write_csv(data, "data/rv_data_111220.csv")
