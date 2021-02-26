#Daniel Suh and Andrew Park
#6/3/20

#Figure 2 in rv_cc manuscript

library(here)

source(knitr::purl(here("/scripts/data_format.Rmd"), quiet=TRUE))


figure_2 <- site_scores %>% ggplot(.,aes(x=pc1Rank,y=cc))+
  geom_point(aes(color = Month, shape = factor(WetAltID), size = size))+
  labs(x="Principal Component 1 Rank", y = "Community Competence (CC)", size = "Size") + 
  scale_shape_manual(values = rep(1:20, len = 20)) +
  scale_size_continuous(range = c(2,10)) +
  guides(shape=F) +
  theme_classic() +
  theme(aspect.ratio=1/1.618, legend.position = c(0.85,0.4), legend.box = "horizontal") +
  scale_color_viridis_d(direction = -1)

figure_2
  
ggsave("fig_2.png",plot=figure_2,device="png",path=here("/figures"))
