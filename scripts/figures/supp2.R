#Daniel Suh and Andrew Park

#Supplementary figure 2

library(here)

source(here("base","src.R"))



comm_data <- readRDS(here("processed_data","comm_data.rds"))

axis_title_size = 15
axis_text_size = 12

prev_cc <- comm_data

order <- prev_cc[order(prev_cc$WetAltID, prev_cc$Month.1),]

#this removes months out of sequence
order <- order[-c(4,14,47,90),]

order %<>% add_column(prev_ratio = NA)

order %<>% group_by(WetAltID) %>% mutate(max_prev = max(Prevalence)) %>% ungroup()

#order %<>% filter(!WetAltID %in% c(8,11,20))


for(n in 2:nrow(order)){
  if (order$Prevalence[n] > order$Prevalence[n-1]) {
    order$prev_ratio[n-1] <- (order$Prevalence[n] - order$Prevalence[n-1]) / ((order$max_prev[n-1]) - order$Prevalence[n-1])
  } else if (order$Prevalence[n] < order$Prevalence[n-1]) {
    order$prev_ratio[n-1] <- 1 - ((order$Prevalence[n-1] - order$Prevalence[n]) / (order$Prevalence[n-1]))
  } else {
    order$prev_ratio[n-1] <- 0
  }
}



#remove the first entry for each wetland to remove the carryover from the last wetland
clean <- order %>% 
  mutate(month_n = gsub("Month", "", Month.1)) %>% 
  group_by(WetAltID) %>% 
  filter(month_n != max(month_n))

order %>% ggplot(.,aes(x=Month.1, y=Prevalence, size = prev_ratio)) + geom_point() + facet_wrap(vars(WetAltID))
clean %>% ggplot(.,aes(x=Month.1, y=Prevalence, size = prev_ratio)) + geom_point() + facet_wrap(vars(WetAltID))

clean %>% ggplot(.,aes(x=Prevalence, y=prev_ratio, color = Month.1)) + geom_point() + facet_wrap(vars(WetAltID))


#plot cleaned plot with lag
cc_corr <- clean %>% ggplot(.,aes(x=cc, y=prev_ratio)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Community Competence", y = expression("Prevalence Ratio" ~(Theta))) +
  ylim(0,0.65) +
  theme(axis.title = element_text(size=axis_title_size),
        axis.text = element_text(size=axis_text_size))



m1 <- cor.test(clean$cc,clean$prev_ratio,method="spearman")


size_corr <- clean %>% ggplot(.,aes(x=log10(size), y=prev_ratio)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "log10(Host Abundance)", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), 
        axis.title = element_text(size=axis_title_size),
        axis.text = element_text(size=axis_text_size)) +
  ylim(0,0.65)

m2 <- cor.test(log10(clean$size),clean$prev_ratio,method="spearman")


temp_corr <- clean %>% ggplot(.,aes(x=MeanWaterTempPredC, y=prev_ratio)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Mean Water Temp", y = "") +
  theme(axis.text.y = element_blank(), axis.ticks.y=element_blank(), 
        axis.title = element_text(size=axis_title_size),
        axis.text = element_text(size=axis_text_size)) +
  ylim(0,0.65)

m3 <- cor.test(clean$MeanWaterTempPredC,clean$prev_ratio,method="spearman")

#multiple comparisons test
pvals <- c(m1$p.value,m2$p.value,m3$p.value)
rho <- c(m1$estimate, m2$estimate, m3$estimate)
rho
pvals
p.adjust(pvals,method="holm")

#plot everything together with patchwork
supp_2 <- cc_corr| size_corr| temp_corr
supp_2
ggsave("supp2a.png",plot=supp_2,width = outwidth[1], height = outwidth[1]/golden,device="png",path=here("figures"))

m4 <- lm(prev_ratio ~ MeanWaterTempPredC + size + cc, data = clean)






clean_alt <- clean %>% filter(WetAltID %in% c(3, 8, 10, 11, 14, 15, 16, 17, 20))

m1 <- cor.test(clean_alt$cc,clean_alt$prev_ratio,method="spearman")

m2 <- cor.test(log10(clean_alt$size),clean_alt$prev_ratio,method="spearman")

m3 <- cor.test(clean_alt$MeanWaterTempPredC,clean_alt$prev_ratio,method="spearman")

#multiple comparisons test
pvals <- c(m1$p.value,m2$p.value,m3$p.value)
rho <- c(m1$estimate, m2$estimate, m3$estimate)
rho
pvals
p.adjust(pvals,method="holm")



clean_alt <- order %>% 
  mutate(month_n = gsub("Month", "", Month.1)) %>% 
  group_by(WetAltID) %>% 
  mutate(prev_ratio_alt = case_when(month_n < max(month_n) ~ prev_ratio,
                                    month_n == max(month_n) ~ 0)) #set last month's prev ratio to zero always


supp_prev <- clean_alt %>% 
  ggplot(.,aes(x=Month, y=Prevalence, size = prev_ratio_alt)) + 
  geom_point() + 
  facet_wrap(vars(WetAltID)) +
  labs(size = "Prevalence\nRatio") + 
  theme(text = element_text(size=10))

supp_prev

ggsave("supp2b.png",
       plot=supp_prev,
       width = outwidth[1], 
       height = outwidth[1]/golden,
       device="png",
       path=here("figures"))









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
