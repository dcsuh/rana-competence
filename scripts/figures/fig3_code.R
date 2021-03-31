#Script Originator: Daniel Suh

#This script creates figure 3

library(here)

source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))


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
  labs(x="Evenness", y = "Community Competence (CC)", size = "Abundance") +
  annotate("text", label = "i", family="Times", fontface="italic", x = 0.134, y = 46000, size = 4, colour = "white")+
  annotate("text", label = "ii", family="Times", fontface="italic", x = 0.427, y = 61000, size = 3, colour = "white")


figure_3
ggsave("fig_3.png",plot=figure_3,device="png",path=here("/figures"))

