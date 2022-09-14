#Daniel Suh and Andrew Park
#6/3/20

#Figure 5 in rv_cc manuscript

library(here)

source(here("base","src.R"))


comm_data <- readRDS(here("processed_data","comm_data.rds"))

axis_title_size = 18
axis_text_size = 15

prev_cc <- comm_data

order <- prev_cc[order(prev_cc$WetAltID, prev_cc$Month.1),]

#this removes months out of sequence
order <- order[-c(4,14,47,90),]

order %<>% add_column(lag_cc = NA, lag_size = NA, lag_temp = NA)


#create new column that includes previous month's value for cc
for(n in 2:nrow(order)){
  order$lag_cc[n] <- order$cc[n-1]
  order$lag_size[n] <- order$size[n-1]
  order$lag_temp[n] <- order$MeanWaterTempPredC[n-1]
}

#remove the first entry for each wetland to remove the carryover from the last wetland
clean <- order %>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)

#plot cleaned plot with lag
cc_corr <- clean %>% ggplot(.,aes(x=lag_cc, y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Community Competence (CC)", y = "Prevalence (t+1)") +
  ylim(0,0.65) +
  theme(axis.title = element_text(size=axis_title_size),
        axis.text = element_text(size=axis_text_size))



m1 <- cor.test(clean$lag_cc,clean$Prevalence,method="spearman")


size_corr <- clean %>% ggplot(.,aes(x=log10(lag_size), y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "log10(Host Abundance)", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), 
        axis.title = element_text(size=axis_title_size),
        axis.text = element_text(size=axis_text_size)) +
  ylim(0,0.65)

m2 <- cor.test(clean$lag_size,clean$Prevalence,method="spearman")


temp_corr <- clean %>% ggplot(.,aes(x=lag_temp, y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Mean Water Temp", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), 
        axis.title = element_text(size=axis_title_size),
        axis.text = element_text(size=axis_text_size)) +
  ylim(0,0.65)

m3 <- cor.test(clean$lag_temp,clean$Prevalence,method="spearman")

#multiple comparisons test
pvals <- c(m1$p.value,m2$p.value,m3$p.value)
pvals
p.adjust(pvals,method="holm")

#plot everything together with patchwork
figure_5 <- cc_corr| size_corr| temp_corr
figure_5
ggsave("fig_5.png",plot=figure_5+ proj_theme,width = outwidth[1], height = outwidth[1]/golden,device="png",path=here("figures"))



#test with 0's filtered out

# tmp <- clean %>% filter(Prevalence>0)
# cor.test(tmp$lag_cc,tmp$Prevalence,method="spearm")
# cor.test(tmp$lag_size,tmp$Prevalence,method="spearm")
# cor.test(tmp$lag_temp,tmp$Prevalence,method="spearm")
# 
# tmp %>% ggplot(.,aes(x=lag_temp, y=Prevalence)) +
#   geom_point() +
#   theme_classic() + geom_smooth(method = "lm") +
#   labs(title = "", x = "Mean Water Temp", y = "") +
#   theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
#   ylim(0,0.65)

