#Script originator: Daniel Suh

#This script recalculates prevalence of ranavirus at each site by estimating a prevalence value that is weighted by overall abundances at that site
#i.e. sum of ((known prevalence in species n) X (relative abundance of species n)) / (total individuals in community) = weighted prevalence

#count number of individuals of species at a sitexmonth that were taken to lab for sampling
#count number of those that were found to be infected
#calculate %prev for that species at sitexmonth
#multiply that %prev by abundance of that species at sitexmonth
#repeat for all species at that sitexmonth
#add predicted number of infected individuals for sitexmonth
#divide that sum by total number of individuals for sitexmonth
#record prevalence for each site

#library
library(tidyverse)
library(magrittr)
library(here)

#load data
#data <- read_csv("competence.csv")
data <- read_csv("data/raw_data/competence_111220.csv")

tmp <- data %>% 
  group_by(WetAltID, Month.1, Species) %>% 
  summarize(infected = sum(RV.Status), total = n()) 
tmp %<>%
  add_column(., prev = tmp$infected/tmp$total) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) 

abundance <- data %>% 
  select(WetAltID, Month.1, Species, AB2:AB42) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>%
  distinct()

q <- full_join(tmp, abundance, by = "siteID") %>% 
  select(-c(WetAltID.y, Month.1.y, Species.y)) %>% 
  distinct() %>%
  filter(Species.x != 32)
q %<>%
  add_column(., ind_inf = q$prev) %>%
  mutate(., total = AB26 + AB42 + AB21 + AB9+ AB2 + AB3 + AB4 + AB5 + AB6 + AB8 + 
           AB20 + AB24 + AB27 + AB28 + AB29 + AB31 + AB34 + AB35 + AB38 + AB39 + AB41)

#this loops through each species at each sitexmonth
for (i in 1:nrow(q)) {
  col <- paste0("AB",q$Species.x[i]) #identifies relevant species and makes temporary string
  q$ind_inf[i] <- as.numeric(q[i,6]*q[i,col]) #calculates estimated prev of rv for relevant species
}

prevalence <- q %>%
  group_by(siteID) %>%
  summarize(total_infected = sum(ind_inf)) 

y <- full_join(prevalence, q, by = "siteID")

y %<>% select(siteID, total_infected, total) %>% 
  add_column(., prevalence = y$total_infected/y$total) %>% 
  distinct()

dat <- data %>% mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) 

new <- full_join(dat, y, by = "siteID")
new$diff <- new$prevalence - new$Prevalence

export <- full_join(dat, y, by = "siteID")
export %<>% select(-Prevalence)
export %<>% rename(Prevalence = prevalence)


#write_csv(export,path="weighted_prev_competence.csv")
write_csv(export,path="data/weighted_prev_competence_111220.csv")
