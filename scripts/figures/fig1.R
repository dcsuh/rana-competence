#Daniel Suh and Andrew Park

#Figure 1 in rv_cc manuscript

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

#create new column that is a ratio between this and last month's prevalence and scales from 0 to 1
#if next month prev is higher, then range is (0,1)
#if next month prev is lower, then range is still (0,1)
# for(n in 2:nrow(order)){
#   if (order$Prevalence[n] > order$Prevalence[n-1] && order$Prevalence[n-1] > 0) {
#     order$prev_ratio[n-1] <- (order$Prevalence[n] - order$Prevalence[n-1]) / order$Prevalence[n-1]
#   } else if (order$Prevalence[n-1] == 0 && order$Prevalence[n] > 0) {
#     order$prev_ratio[n-1] <- order$Prevalence[n]
#   } else if (order$Prevalence[n] < order$Prevalence[n-1] && order$Prevalence[n] > 0) {
#     order$prev_ratio[n-1] <- (order$Prevalence[n-1] - order$Prevalence[n]) / order$Prevalence[n]
#   } else if (order$Prevalence[n] == 0 && order$Prevalence[n-1] > 0) {
#     order$prev_ratio[n-1] <- 1-order$Prevalence[n-1]
#   } else if (order$Prevalence[n-1] > 0 && order$Prevalence[n] == order$Prevalence[n-1]) {
#     order$prev_ratio[n-1] <- 1
#   } else {
#     order$prev_ratio[n-1] <- 0
#   }
# }

# for(n in 2:nrow(order)){
#   if (order$Prevalence[n] > order$Prevalence[n-1]) {
#     order$prev_ratio[n-1] <- (order$Prevalence[n] - order$Prevalence[n-1]) / (order$max_prev[n-1])
#   } else if (order$Prevalence[n] < order$Prevalence[n-1]) {
#     order$prev_ratio[n-1] <- 1 - ((order$Prevalence[n-1] - order$Prevalence[n]) / order$max_prev[n-1])
#   } else {
#     order$prev_ratio[n-1] <- 0
#   }
# }

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
clean <- order %>% mutate(month_n = gsub("Month", "", Month.1)) %>% group_by(WetAltID) %>% filter(month_n != max(month_n))

order %>% ggplot(.,aes(x=Month.1, y=Prevalence, size = prev_ratio)) + geom_point() + facet_wrap(vars(WetAltID))
clean %>% ggplot(.,aes(x=Month.1, y=Prevalence, size = prev_ratio)) + geom_point() + facet_wrap(vars(WetAltID))


#plot cleaned plot with lag
cc_corr <- clean %>% ggplot(.,aes(x=cc, y=prev_ratio)) +
  geom_point() +
  theme_classic() + geom_smooth(method = "lm") +
  labs(title = "", x = "Community Competence", y = "Prevalence Ratio") +
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

m2 <- cor.test(clean$size,clean$prev_ratio,method="spearman")


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
pvals
p.adjust(pvals,method="holm")

#plot everything together with patchwork
figure_2 <- cc_corr| size_corr| temp_corr
figure_2
ggsave("fig2.png",plot=figure_2,width = outwidth[1], height = outwidth[1]/golden,device="png",path=here("figures"))

m4 <- lm(prev_ratio ~ MeanWaterTempPredC + size + cc, data = clean)


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
