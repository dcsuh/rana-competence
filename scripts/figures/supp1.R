## Daniel Suh

#This script creates supplementary figure 1: visual of the estimated viral load for each species

library(here)

source(here("base","src.R"))


vl <- readRDS(here("processed_data","vl.rds"))

vl %<>% filter(m>3)

supp1 <- vl %>% filter(!is.na(mean)) %>% arrange(mean) %>% 
  mutate(tnrs_name = factor(tnrs_name, levels = tnrs_name)) %>% 
  ggplot(.,aes(x=tnrs_name,y=log10(mean))) +
  geom_point(aes(size=1)) + 
  xlab("Species") + ylab("Viral Load") + 
  geom_hline(yintercept = 1) + 
  labs(title = "Viral Loads") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,face = "italic", size = 12),
        title = element_text(size = 18),
        legend.position = "none") +
  scale_y_continuous(breaks=seq(0,6,1))

supp1_se <- vl %>% filter(!is.na(mean)) %>% arrange(mean) %>% 
  mutate(abb_name = factor(abb_name, levels = abb_name)) %>% 
  ggplot(.,aes(x=abb_name,y=log10(mean))) +
  geom_point(size = 2) + 
  xlab("Species") + ylab("log(Viral Load)") + 
  geom_linerange(aes(ymin=log10(mean)-log10(se),ymax=log10(mean)+log10(se))) +
  geom_hline(yintercept = 1) + 
  labs(title = "Viral Load") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,face = "italic", size = 12),
        title = element_text(size = 18),
        legend.position = "none") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

supp1_raw <- vl %>% arrange(mean) %>% 
  mutate(abb_name = factor(abb_name, levels = abb_name)) %>% 
  ggplot(.,aes(x=abb_name,y=mean)) +
  geom_point() + 
  xlab("Species") + ylab("Viral Load") + 
  geom_linerange(aes(ymin=mean-se,ymax=mean+se)) +
  labs(title = "Viral Loads") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

supp1_log10 <- vl %>% arrange(mean) %>% 
  mutate(abb_name = factor(abb_name, levels = abb_name)) %>% 
  ggplot(.,aes(x=abb_name,y=log1p_mean)) +
  geom_point() + 
  xlab("Species") + ylab("log1p(Viral Load)") + 
  geom_linerange(aes(ymin=log1p_mean-log1p_se,ymax=log1p_mean+log1p_se)) +
  geom_hline(yintercept = 1) + 
  labs(title = "Viral Loads") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_y_continuous(breaks=seq(0,6,1))

ggsave("supp1.png",plot=supp1_se,device="png",path=here("figures"))
