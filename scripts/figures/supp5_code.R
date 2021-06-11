## Script originator: Daniel Suh
## Date created: Mar. 29, 2021

#This script creates supplementary figure 4: prevalence ~ richness

library(here)


source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))


supp5a <- evenness %>% ggplot(.,aes(x=richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Richness", y="Prevalence")

cor.test(evenness$richness, evenness$Prevalence, method="spearm")

supp5b <- evenness %>% ggplot(.,aes(x=richness, y=cc)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Richness", y="Community Competence")

cor.test(evenness$richness, evenness$cc, method="spearm")

supp5c <- lag_evenness %>% ggplot(.,aes(x=lag_richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Richness", y="Prevalence (t+1)")

cor.test(lag_evenness$lag_richness, lag_evenness$Prevalence, method="spearm")

evenness %>% ggplot(.,aes(x=J, y=cc)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Evenness", y="Community Competence")

evenness %>% ggplot(.,aes(x=J, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Evenness", y="Prevalence")

cor.test(evenness$J, evenness$cc, method="spearm")

supp5d <- lag_evenness %>% ggplot(.,aes(x=lag_J, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Evenness", y="Prevalence (t+1)")

cor.test(lag_evenness$lag_J, lag_evenness$Prevalence, method="spearm")

ggsave("supp5a.png",plot=supp5a,device="png",path=here("figures"))
ggsave("supp5b.png",plot=supp5b,device="png",path=here("figures"))
ggsave("supp5c.png",plot=supp5c,device="png",path=here("figures"))
ggsave("supp5d.png",plot=supp5d,device="png",path=here("figures"))

