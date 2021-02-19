## -----------------------------------------------------------------------------------
library(tidyverse)
library(magrittr)
library(here)


## -----------------------------------------------------------------------------------
here()
data <- read_csv(here("/data/weighted_prev_competence_111220.csv"))
names <- read_csv(here("/data/species_names_ids.csv"))


## -----------------------------------------------------------------------------------
vl <- data %>% dplyr::select(vl.4:vl.42) %>% distinct() %>% pivot_longer(vl.4:vl.42)
vl$species_code <- as.double(gsub("vl.","",vl$name))
nrow(vl)
vl %<>% full_join(., names, by="species_code")
vl %<>% mutate(., ln_value = log(value)) %>% mutate(., log10_value = log10(value))
sum(vl$value>0)

