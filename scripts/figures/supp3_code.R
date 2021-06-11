## Daniel Suh


#supplementary figures for cc values by month and size


library(here)


source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))

evenness$Month <- factor(evenness$Month,levels=c("Feb","Mar","Apr","May","Jun","Jul"))

supp5a <- evenness %>% ggplot(.,aes(x=Month,y=cc)) + geom_boxplot() + geom_point() + ylab("Community Competence")

evenness %>% ggplot(.,aes(x=Month,y=cc)) + geom_boxplot() + geom_point() + geom_line(aes(group=WetAltID))

supp5b <- evenness %>% ggplot(.,aes(x=log10(size),y=cc)) + geom_point() + geom_smooth(method="lm")
cor.test(log10(evenness$size), evenness$cc, method = "spearman")
  
ggsave("supp3a.png",plot=supp3a,device="png",path=here("figures"))
ggsave("supp3b.png",plot=supp3b,device="png",path=here("figures"))
