#Daniel Suh and Andrew Park
#6/3/20

#Figure 1 in rv_cc manuscript

library(here)

source(here("base","src.R"))

comm_data <- readRDS(here("processed_data","comm_data.rds"))


figure_1 <- comm_data %>% ggplot(.,aes(x=pc1Rank,y=cc))+
  geom_point(aes(color = Month, shape = factor(WetAltID), size = size))+
  labs(x="Principal Component 1 Rank", y = "Community Competence", size = "Host\nAbundance") + 
  scale_shape_manual(values = rep(1:20, len = 20)) +
  scale_size_continuous(range = c(2,10)) +
  guides(shape=F) +
  theme_classic() +
  theme(aspect.ratio=1/1.618, legend.position = c(0.85,0.85), legend.box = "horizontal") +
  scale_color_viridis_d(direction = -1) +
  theme(text = element_text(size=18))

comm_data %>% ggplot(., aes(x=pc1Rank, y=cc)) +
  geom_point(aes(color = factor(WetAltID), shape = factor(WetAltID), size = size)) + 
  facet_wrap(vars(Month.1), nrow = 6)+
  scale_color_manual(values = rep(1:20, len = 20)) +
  scale_shape_manual(values = rep(1:20, len = 20)) +
  scale_size_continuous(range = c(2,10)) +
  guides(shape=F) +
  labs(x="Principal Component 1 Rank (Community Similarity)", y = "Community Competence", size = "Host\nAbundance") + 
  theme_classic() +
  theme(text = element_text(size=18))
  
figure_1
  
ggsave("fig1.png",plot=figure_1,device="png",height = 20, width = 1.618*20, units = "cm",path=here("figures"))
