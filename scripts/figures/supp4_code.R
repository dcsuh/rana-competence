## Script originator: Daniel Suh
## Date created: Mar. 29, 2021

#This script creates supplementary figure 4: prevalence ~ richness

library(here)


source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))


evenness %>% ggplot(.,aes(x=richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm")

cor.test(evenness$richness, evenness$Prevalence, method="spearm")

supp4a <- evenness %>% ggplot(.,aes(x=richness, y=cc)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Richness", y="Community Competence")

cor.test(evenness$richness, evenness$cc, method="spearm")

supp4b <- lag_evenness %>% ggplot(.,aes(x=lag_richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Richness", y="Prevalence (t+1)")

cor.test(lag_evenness$lag_richness, lag_evenness$Prevalence, method="spearm")

supp4c <- evenness %>% ggplot(.,aes(x=J, y=cc)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Evenness", y="Community Competence")

cor.test(evenness$J, evenness$cc, method="spearm")

supp4d <- lag_evenness %>% ggplot(.,aes(x=lag_J, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Evenness", y="Prevalence (t+1)")

cor.test(lag_evenness$lag_J, lag_evenness$Prevalence, method="spearm")

#ggsave("supp4a.png",plot=supp4a,device="png",path=here("figures"))
#ggsave("supp4b.png",plot=supp4b,device="png",path=here("figures"))

