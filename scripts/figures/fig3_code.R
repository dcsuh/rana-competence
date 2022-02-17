#Script Originator: Daniel Suh

#This script creates figure 3

library(here)

source(here("base","src.R"))

comm_data <- readRDS(here("processed_data","comm_data.rds"))

comm_data <- comm_data[order(comm_data$WetAltID),]


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

low_cc <- lag_evenness %>% filter(low==1)
high_cc <- lag_evenness %>% filter(high==1)

mySpan=1.5
myColor="gray80"
figure_3 <- ggplot(data=lag_evenness,aes(x=lag_J, y=lag_cc)) +
  geom_point(data = lag_evenness, aes(x=lag_J, y=lag_cc, color = Prevalence, size = lag_size)) +
  scale_color_viridis_c("Prevalence\n(t+1)", direction = -1) +
  new_scale_color() +
  geom_smooth(data = low_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  new_scale_color() +
  geom_smooth(data = high_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  theme_minimal()+ylim(0,3e+05)+
  labs(x="Evenness", y = "Community Competence (CC)", size = "Host\nAbundance") +
  annotate("text", label = "i", family="Times", fontface="italic", x = 0.134, y = 46000, size = 4, colour = "white")+
  annotate("text", label = "ii", family="Times", fontface="italic", x = 0.427, y = 61000, size = 3, colour = "white")+
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=12),
        legend.title = element_text(size=15),
        legend.text = element_text(size=10))

figure_3
ggsave("fig_3.png",plot=figure_3,device="png",path=here("figures"))



