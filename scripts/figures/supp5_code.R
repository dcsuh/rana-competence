## Script originator: Daniel Suh
## Date created: Mar. 29, 2021

#This script creates supplementary figure 4: prevalence ~ richness

library(here)

source(here("base","src.R"))

comm_data <- readRDS(here("processed_data","comm_data.rds"))


#create new df for ordered sites to test lagged richness and evenness on prevalence
lag_evenness <- comm_data[-c(4,14,47,90),] #brute force remove site_months out of sequence
lag_evenness$lag_richness <- c(0)
lag_evenness$lag_J <- c(0)
lag_evenness$lag_cc <- c(0)
lag_evenness$lag_size <- c(0)
for(n in 2:92){
  lag_evenness$lag_richness[n] <- lag_evenness$richness[n-1]
  lag_evenness$lag_J[n] <- lag_evenness$J[n-1]
  lag_evenness$lag_cc[n] <- lag_evenness$cc[n-1]
  lag_evenness$lag_size[n] <- lag_evenness$size[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
lag_evenness %<>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)

#fit two lines to cc_evenness_prev plot
lag_evenness$high <- c(0)
lag_evenness$low <- c(0)
for (i in 1:nrow(lag_evenness)){
  if (lag_evenness$lag_J[i] > 0.6){
    lag_evenness$high[i] <- 1
    lag_evenness$low[i] <- 1
  }
  else if (lag_evenness$lag_cc[i] < 45000){
    lag_evenness$low[i] <- 1
  }
  else if (lag_evenness$lag_cc[i] > 45000){
    lag_evenness$high[i] <- 1
  }
  else {
    lag_evenness$high[i] <- 0
    lag_evenness$low[i] <- 0
  }
}

supp5a <- comm_data %>% ggplot(.,aes(x=richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Richness", y="Prevalence")

cor.test(comm_data$richness, comm_data$Prevalence, method="spearm")

supp5b <- comm_data %>% ggplot(.,aes(x=richness, y=cc)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Richness", y="Community Competence")

cor.test(comm_data$richness, comm_data$cc, method="spearm")

supp5c <- lag_evenness %>% ggplot(.,aes(x=lag_richness, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Richness", y="Prevalence (t+1)")

cor.test(lag_evenness$lag_richness, lag_evenness$Prevalence, method="spearm")

comm_data %>% ggplot(.,aes(x=J, y=cc)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Evenness", y="Community Competence")

comm_data %>% ggplot(.,aes(x=J, y=Prevalence)) +
  geom_point() +
  geom_smooth(method="lm") +
  theme_classic() +
  labs(x="Evenness", y="Prevalence")

cor.test(comm_data$J, comm_data$cc, method="spearm")

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

