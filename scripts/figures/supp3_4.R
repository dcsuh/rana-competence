## Daniel Suh


#supplementary figures for cc values by month and size


library(here)

source(here("base","src.R"))


comm_data <- readRDS(here("processed_data","comm_data.rds"))

comm_data$Month <- factor(comm_data$Month,levels=c("Feb","Mar","Apr","May","Jun","Jul"))

supp3a <- comm_data %>% ggplot(.,aes(x=Month,y=cc)) + 
  geom_boxplot() + 
  stat_summary(fun = mean, geom = "point", size = 4, color = "red") + 
  ylab("Community Competence")

comm_data %>% ggplot(.,aes(x=Month,y=cc)) + geom_boxplot() + geom_point() + geom_line(aes(group=WetAltID))

supp3b <- comm_data %>% ggplot(.,aes(x=log10(size),y=cc)) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  xlab("log10(Host Abundance)") + 
  ylab("Community Competence")
m1 <- cor.test(log10(comm_data$size), comm_data$cc, method = "spearman")

supp3c <- comm_data %>% ggplot(.,aes(x=MeanWaterTempPredC,y=cc)) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  xlab("Mean Water Temperature") + 
  ylab("Community\nCompetence")
m2 <- cor.test(comm_data$MeanWaterTempPredC, comm_data$cc, method = "spearman")

supp3d <- comm_data %>% ggplot(.,aes(x=MeanWaterTempPredC,y=log10(size))) + 
  geom_point() + 
  geom_smooth(method="lm") + 
  xlab("Mean Water Temperature") + 
  ylab("log10(Host Abundance)")
m3 <- cor.test(comm_data$MeanWaterTempPredC, log10(comm_data$size), method = "spearman")


pvals <- c(m1$p.value,m2$p.value, m3$p.value)
pvals
p.adjust(pvals,method="holm")
  
ggsave("supp3.png",plot=supp3a,width=outwidth[1],scale=golden,units=unit,device="png",path=here("figures"))
# ggsave("supp3b.png",plot=supp3b,width=outwidth[1],scale=golden,units=unit,device="png",path=here("figures"))
# ggsave("supp3c.png",plot=supp3c,width=outwidth[1],scale=golden,units=unit,device="png",path=here("figures"))
# ggsave("supp3d.png",plot=supp3d,width=outwidth[1],scale=golden,units=unit,device="png",path=here("figures"))

supp4 <- supp3b + supp3c / supp3d
ggsave("supp4.png",plot=supp4,width=outwidth[1],scale=golden,units=unit,device="png",path=here("figures"))
