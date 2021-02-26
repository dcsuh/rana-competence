#Daniel Suh and Andrew Park
#6/3/20

#Figure 5 in rv_cc manuscript

library(here)


source(knitr::purl(here("/scripts/data_format.Rmd"), quiet=TRUE))


#plot cleaned plot with lag
cc_corr <- clean %>% ggplot(.,aes(x=lag_cc, y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Community Competence", y = "Prevalence") +
  ylim(0,0.65)

#cor.test(clean$lag_cc,clean$Prevalence,method="spearman")


size_corr <- clean %>% ggplot(.,aes(x=log10(lag_size), y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "log(Community Size)", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  ylim(0,0.65)

#cor.test(clean$lag_size,clean$Prevalence,method="spearman")


temp_corr <- clean %>% ggplot(.,aes(x=lag_temp, y=Prevalence)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Mean Water Temp", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) +
  ylim(0,0.65)

#cor.test(clean$lag_temp,clean$Prevalence,method="spearman")


#plot everything together with patchwork
figure_5 <- cc_corr| size_corr| temp_corr
ggsave("fig_5.png",plot=figure_5,device="png",path=here("/figures"))

