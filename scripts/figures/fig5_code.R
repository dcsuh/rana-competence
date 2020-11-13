#Daniel Suh and Andrew Park
#6/3/20

#Figure 5 in rv_cc manuscript

#install_github("thomasp85/patchwork")
library(tidyverse)
library(magrittr)
library(patchwork)
library(here)


data <- read_csv("data/weighted_prev_competence_111220.csv")

prev_cc <- data %>% dplyr::select(WetAltID, Month.1, cc, Prevalence, Month, AB2:AB42, MeanWaterTempPredC)
prev_cc %<>% dplyr::select(-ABAmb, -ABSal) %>%  mutate(., total = rowSums(.[6:26])) #get totals for community size
prev_cc %<>% dplyr::select(-AB2:-AB42)
prev_cc %<>% distinct()

order <- prev_cc[order(prev_cc$WetAltID, prev_cc$Month.1),]

#this removes months out of sequence
order <- order[-c(4,14,47,90),]

order$lag_cc <- order$cc
order$lag_size <- order$total
order$lag_temp <- order$MeanWaterTempPredC

#create new column that includes previous month's value for cc
for(n in 2:nrow(order)){
  order$lag_cc[n] <- order$cc[n-1]
  order$lag_size[n] <- order$total[n-1]
  order$lag_temp[n] <- order$MeanWaterTempPredC[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
clean <- order %>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)

#plot cleaned plot with lag
cc_corr <- clean %>% filter(Prevalence>0) %>% ggplot(.,aes(x=lag_cc, y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Community Competence", y = "Prevalence")

#correlation test
cor.test(clean$lag_cc,clean$Prevalence,method="spearman")


size_corr <- clean %>% filter(Prevalence>0) %>% ggplot(.,aes(x=log(lag_size), y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "ln(Community Size)", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())

cor.test(clean$lag_size,clean$Prevalence,method="spearman")


temp_corr <- clean %>% filter(Prevalence>0) %>% ggplot(.,aes(x=MeanWaterTempPredC, y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Mean Water Temp", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())


cor.test(clean$lag_temp,clean$Prevalence,method="spearman")


#plot everything together with patchwork
corr_plots <- cc_corr| size_corr| temp_corr
corr_plots
