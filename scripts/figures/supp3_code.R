## Daniel Suh


#supplementary figures for cc values by month and size


library(here)

source(here("base","src.R"))


comm_data <- readRDS(here("processed_data","comm_data.rds"))

comm_data$Month <- factor(comm_data$Month,levels=c("Feb","Mar","Apr","May","Jun","Jul"))

supp3a <- comm_data %>% ggplot(.,aes(x=Month,y=cc)) + geom_boxplot() + geom_point() + ylab("Community Competence")

comm_data %>% ggplot(.,aes(x=Month,y=cc)) + geom_boxplot() + geom_point() + geom_line(aes(group=WetAltID))

supp3b <- comm_data %>% ggplot(.,aes(x=log10(size),y=cc)) + geom_point() + geom_smooth(method="lm")
cor.test(log10(comm_data$size), comm_data$cc, method = "spearman")

supp3c <- comm_data %>% ggplot(.,aes(x=MeanWaterTempPredC,y=cc)) + geom_point() + geom_smooth(method="lm")
cor.test(comm_data$MeanWaterTempPredC, comm_data$cc, method = "spearman")

  
ggsave("supp3a.png",plot=supp3a,device="png",path=here("figures"))
ggsave("supp3b.png",plot=supp3b,device="png",path=here("figures"))
ggsave("supp3c.png",plot=supp3c,device="png",path=here("figures"))
