## Script originator: Daniel Suh
## Date created: Mar. 29, 2021

#This script creates supplementary figure 4: prevalence ~ richness

library(here)


source(knitr::purl(here("/scripts/data_format.Rmd"), quiet=TRUE))


evenness %>% ggplot(.,aes(x=richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm")

cor.test(evenness$richness, evenness$Prevalence, method="spearm")
