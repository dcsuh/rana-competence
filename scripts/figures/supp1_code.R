## Script originator: Daniel Suh
## Date created: Feb. 19, 2021

#This script creates supplementary figure 1: visual of the estimated viral load for each species

library(here)


source(knitr::purl(here("/scripts/data_format.Rmd"), quiet=TRUE))

vl %>% arrange(log10_value) %>% 
  mutate(abb_name = factor(abb_name, levels = abb_name)) %>% 
  ggplot(.,aes(x=abb_name,y=log10_value)) +
  geom_point() + 
  xlab("Species") + ylab("log10(Viral Load)") + 
  geom_hline(yintercept = log10(20)) + 
  labs(title = "Viral Loads") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
