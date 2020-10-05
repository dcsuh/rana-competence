#Daniel Suh and Andrew Park
#6/3/20

#Figure 5 in rv_cc manuscript

#install_github("thomasp85/patchwork")
library(tidyverse)
library(magrittr)
library(patchwork)


data <- read_csv("weighted_prev_competence.csv")

prev_cc <- data %>% dplyr::select(WetAltID, Month.1, cc, Prevalence, Month)
prev_cc %<>% distinct() %>% distinct() %>% remove_missing()

order <- prev_cc[order(prev_cc$WetAltID, prev_cc$Month.1),]

#this removes months out of sequence
order <- order[-c(4,37,38,45,85),]

order$lag <- order$cc

#create new column that includes previous month's value for cc
for(n in 2:nrow(order)){
  order$lag[n] <- order$cc[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
clean <- order %>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)

#plot cleaned plot with lag
cc_corr <- clean %>% filter(Prevalence>0) %>% ggplot(.,aes(x=lag, y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Community Competence", y = "Prevalence")

#correlation test
cor.test(clean$lag,clean$Prevalence,method="spearman")




#correlation plot for community size and prevalence
comm_size <- data %>% dplyr::select(WetAltID, Month.1, Prevalence, Month, AB2:AB42) %>% distinct() #subset data
comm_size %<>% dplyr::select(-ABAmb, -ABSal) %>%  mutate(., total = rowSums(.[5:25])) #get totals for community size
comm_size %<>% dplyr::select(-AB2:-AB42)

comm_size_order <- comm_size[order(comm_size$WetAltID, comm_size$Month.1),] #order by month
comm_size_order <- comm_size_order[-c(4,37,38,45,85),] #remove out of sequence month
comm_size_order$lag <- comm_size_order$total

#lag by one month
for(n in 2:nrow(comm_size_order)){
  comm_size_order$lag[n] <- comm_size_order$total[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
comm_size_clean <- comm_size_order %>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)

size_corr <- comm_size_clean %>% filter(Prevalence>0) %>% ggplot(.,aes(x=log(lag), y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "ln(Community Size)", y = "Prevalence")

cor.test(comm_size_clean$lag,comm_size_clean$Prevalence,method="spearman")



#correlation plot for mean water temp and prevalence
#mean water temp is a proxy for viral half-life
comm_half <- data %>% dplyr::select(WetAltID, Month.1, Prevalence, Month, MeanWaterTempPredC) %>% distinct() %>% remove_missing() #subset data

comm_half_order <- comm_half[order(comm_half$WetAltID, comm_half$Month.1),] #order by month
comm_half_order <- comm_half_order[-c(4,11,44),] #remove out of sequence months
comm_half_order$lag <- comm_half_order$MeanWaterTempPredC 

#lag by one month
for(n in 2:nrow(comm_half_order)){
  comm_half_order$lag[n] <- comm_half_order$MeanWaterTempPredC[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
comm_half_clean <- comm_half_order %>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)

temp_corr <- comm_half_clean %>% filter(Prevalence>0) %>% ggplot(.,aes(x=lag, y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Mean Water Temp", y = "Prevalence")

cor.test(comm_half_clean$lag,comm_half_clean$Prevalence,method="spearman")


#plot everything together with patchwork
corr_plots <- cc_corr/(size_corr | temp_corr)
corr_plots

cc_corr/ size_corr/ temp_corr
