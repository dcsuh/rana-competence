## Script originator: Daniel Suh
## Date created: Mar. 29, 2021

#This script creates supplementary figure 4: prevalence ~ richness

library(here)


source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))


evenness %>% ggplot(.,aes(x=richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm")

cor.test(evenness$richness, evenness$Prevalence, method="spearm")


supp4 <- lag_evenness %>% ggplot(.,aes(x=lag_richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm")

cor.test(lag_evenness$lag_richness, lag_evenness$Prevalence, method="spearm")

ggsave("supp4.png",plot=supp4,device="png",path=here("/figures"))
