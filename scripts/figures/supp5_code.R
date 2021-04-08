## Daniel Suh


#supplementary figures for cc values by month and size


library(here)


source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))

evenness$Month <- factor(evenness$Month,levels=c("Feb","Mar","Apr","May","Jun","Jul"))

supp5 <- evenness %>% ggplot(.,aes(x=Month,y=cc)) + geom_boxplot() + geom_point() + ylab("Community Competence")

evenness %>% ggplot(.,aes(x=Month,y=cc)) + geom_boxplot() + geom_point() + geom_line(aes(group=WetAltID))

supp6 <- evenness %>% ggplot(.,aes(x=log10(size),y=cc)) + geom_point() + geom_smooth(method="lm")
cor.test(log10(evenness$size), evenness$cc, method = "spearman")
  
ggsave("supp5.png",plot=supp5,device="png",path=here("/figures"))
ggsave("supp6.png",plot=supp6,device="png",path=here("/figures"))
