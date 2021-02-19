#Script Originator: Daniel Suh

#This script creates figure 3


library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(patchwork)
library(ggnewscale)
library(here)



p1 <- evenness %>% remove_missing() %>% ggplot(aes(x=richness, y=cc)) + 
  geom_smooth(method ="lm") +
  geom_point() +
  labs(x = "Richness", y = "CC", title = ) +
  theme_classic()

cor.test(richness_cc$cc, richness_cc$richness, method="spearman")

p2 <- lag_evenness %>% remove_missing() %>% ggplot(aes(x=lag_richness, y=Prevalence)) + 
  geom_smooth(method ="lm") +
  geom_point() +
  labs(x = "Richness", y = "Prevalence", title = ) +
  theme_classic()

cor.test(x = lag_evenness$lag_richness, y = lag_evenness$Prevalence, method = "spearman")

#make plot w/ linear model using regular (not lagged) J and cc
p3 <- lag_evenness %>% ggplot(aes(x = J, y = cc)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  labs(x="Evenness", y = "CC") +
  theme_minimal()+ylim(0,3e+05)
  
cor.test(y = lag_evenness$cc, x = lag_evenness$J, method = "spearman")

mySpan=1.5
myColor="gray80"
p4 <- ggplot(data=lag_evenness,aes(x=lag_J, y=lag_cc)) +
  geom_point(data = lag_evenness, aes(x=lag_J, y=lag_cc, color = Prevalence, size = lag_size)) +
  scale_color_viridis_c("Prevalence\n(t+1)", direction = -1) +
  new_scale_color() +
  geom_smooth(data = low_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  new_scale_color() +
  geom_smooth(data = high_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  theme_minimal()+ylim(0,3e+05)+
  labs(x="evenness", y = "community competence", size = "Abundance") +
  annotate("text", label = "i", family="Times", fontface="italic", x = 0.134, y = 46000, size = 4, colour = "white")+
  annotate("text", label = "ii", family="Times", fontface="italic", x = 0.427, y = 61000, size = 3, colour = "white")
  # geom_vline(xintercept = 0.6, linetype = 'dashed', color = "red") +
  # geom_hline(yintercept = 45000, linetype = 'dashed', color = "red")
p4

ggplot(data = lag_evenness, aes(x=richness, y=cc)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="richness", y="community competence")+
  theme_minimal()+ylim(0,3e+05)

ggplot(data=lag_evenness,aes(x=lag_J, y=lag_cc)) +
  geom_point() +
  theme_minimal()+ylim(0,3e+05)+
  labs(x="evenness", y = "community competence")+
  new_scale_color() +
  geom_smooth(data = low_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  new_scale_color() +
  geom_smooth(data = high_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F)+
  geom_vline(xintercept = 0.6, linetype = 'dashed', color = "red") +
  geom_hline(yintercept = 45000, linetype = 'dashed', color = "red")

#make plot w/ linear model using lagged J and cc
ggplot(data=lag_evenness,aes(x=lag_J, y=lag_cc)) +
  geom_point() +
  geom_smooth(method="loess")+
  theme_minimal()+ylim(0,3e+05)+
  labs(x="Evenness", y = "CC")
cor.test(y = lag_evenness$lag_cc, x = lag_evenness$lag_J, method = "spearman")


final <- p4
final
