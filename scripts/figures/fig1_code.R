#Daniel Suh and Andrew Park
#6/3/20

#Figure 1 in rv_cc manuscript

#install_github("thomasp85/patchwork")
library(tidyverse)
library(magrittr)
library(patchwork)
#install.packages("ggnewscale")
library(ggnewscale)


p1 <- reference %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(values = c(NA, "red")) + 
  labs(x = "Environmental", y = "Contact") +
  xlim(2,8) +
  ylim(2,10) + 
  labs(title = "Reference")


p2 <- composition %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(values = c(NA, "red")) + 
  new_scale_color()+
  geom_contour(data = reference, aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(values = c(NA, "gray")) +
  labs(x = "", y = "Contact") +
  xlim(2,8) +
  ylim(2,10) + 
  labs(title = "Composition")


p3 <- size %>% filter(tot != 150) %>% 
  ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(values = c(NA, "orange")) + 
  new_scale_color()+
  geom_contour(data = reference, aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(values = c(NA, "gray")) +
  labs(x = "", y = "") +
  xlim(2,8) +
  ylim(2,10) + 
  labs(title = "Size")


p4 <- halflife %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(values = c(NA, "green")) + 
  new_scale_color()+
  geom_contour(data = reference, aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(values = c(NA, "gray")) +
  labs(x = "Environmental", y = "Contact") +
  xlim(2,8) +
  ylim(2,10) + 
  labs(title = "Half-life")


p5 <- combined %>%
  ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(values = c(NA, "purple")) + 
  new_scale_color()+
  geom_contour(data = reference, aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(values = c(NA, "gray")) +
  labs(x = "Environmental", y = "") +
  xlim(2,8) +
  ylim(2,10) + 
  labs(title = "Combined")

figure_1 <- (p2 | p3)/(p4 | p5)
figure_1
