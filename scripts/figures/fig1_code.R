#Daniel Suh and Andrew Park
#6/3/20

#Figure 1 in rv_cc manuscript


library(here)

source(here("base", "src.R"))

#define function contour(): inputs are parameters and output is values for grid search for R0
contour <- function(trans_1_min, trans_1_max, trans_3_min, trans_3_max, mort1, mort2, degr) {
  #expand data
  data <- expand_grid("trans1" = seq(from = trans_1_min, to = trans_1_max, by = 0.00001),
                      "trans2" = c(0.0001),
                      "trans3" = seq(from = trans_3_min, to = trans_3_max, by = 0.00001),
                      "shed1" = c(rep(1/2, length("trans1"))),
                      "shed2" = c(rep(1/2, length("trans1"))),
                      "reco1" = c(rep(1/10, length("trans1"))),
                      "reco2" = c(rep(1/10, length("trans1"))),
                      "mort1" = c(mort1), 
                      "mort2" = c(mort2), 
                      "degr" = c(degr), 
                      "birth1" = c(rep(5/3, length("trans1"))),
                      "birth2" = c(rep(5/3, length("trans1"))))
  
  #calculate disease-free equilibria for use in next-gen matrix
  data %<>% add_column(dfe1 = (data$birth1/data$mort1),dfe2 = (data$birth2/data$mort2), eigen = NA) 
  #make matrix
  mat <- matrix(1, nrow = 3, ncol = 3)
  
  m <- nrow(data)
  
  #calculate eigenvalue using next-gen matrix. max eigenvalue is R0
  for(n in 1:m){
    mat[1] = (data$trans1[n]*data$dfe1[n])/(data$reco1[n]+data$mort1[n])
    mat[2] = (data$trans1[n]*data$dfe2[n])/(data$reco2[n]+data$mort2[n])
    mat[3] = (data$shed1[n]/data$degr[n])
    mat[4] = (data$trans2[n]*data$dfe1[n])/(data$reco1[n]+data$mort1[n])
    mat[5] = (data$trans2[n]*data$dfe2[n])/(data$reco2[n]+data$mort2[n])
    mat[6] = (data$shed2[n]/data$degr[n])
    mat[7] = (data$trans3[n]*data$dfe1[n])/(data$reco1[n]+data$mort1[n])
    mat[8] = (data$trans3[n]*data$dfe2[n])/(data$reco2[n]+data$mort2[n])
    mat[9] = 0
    eigen <- max(eigen(mat)$values)
    data$eigen[n] <- eigen
  }
  
  data$span1 <- 1/data$mort1 #lifespan for species 1
  data$span2 <- 1/data$mort2 #lifespan for species 2
  data$half <- log(1/2)/-data$degr #viral half-life
  data$epid <- data$eigen > 1 #invasion potential (i.e. R0 > 1)
  data$prop1 <- data$trans1/data$trans2 #proportion of transmission of species 2 over species 1
  data$prop3 <- data$trans3/data$trans2 #proportion of environmental transmission over species 1
  data$prop_mort <- data$mort1/data$mort2 #relative mortality for species 1 and 2
  data$cap1 <- data$birth1/data$mort1 #carrying capacity for species 1
  data$cap2 <- data$birth2/data$mort2 #carrying capacity for species 2
  data$cap_prop <- data$cap1/data$cap2 #relative carrying capacities for species 1 and 2
  data$tot <- data$cap1 + data$cap2 #total abundance for species 1 and 2
  data$cc <- data$trans1*(data$cap1/data$tot) + data$trans2*(data$cap2/data$tot) #cc defined as transmission of a species multiplied by its relative abundance
  output <- data %>% dplyr::select(trans1, trans2, prop1, prop3, trans3, eigen, prop_mort, mort1, mort2, degr, half, cap1, cap2, cap_prop, tot, cc) %>% distinct()
  return(output)
}

reference <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0.0001,trans_3_max = 0.001,mort1 = 1/45,mort2 = 1/45,degr = 1/1.947799)

composition <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0.0001,trans_3_max = 0.001,mort1 = 1/60,mort2 = 1/30,degr = 1/1.947799)

size <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0.0001,trans_3_max = 0.001,mort1 = c(1/52.5),mort2 = c(1/52.5),degr = 1/1.947799)

halflife <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0.0001,trans_3_max = 0.001,mort1 = c(1/45),mort2 = 1/45,degr = 1/3.895598)

combined <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0.0001,trans_3_max = 0.001,mort1 = c(1/70),mort2 = c(1/35),degr = 1/3.895598)


axis_text_size = 8
plot_label_size = 20

p1 <- reference %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(na.value = NA, values = c(NA, "red")) + 
  labs(x = "Environmental", y = "Contact") +
  xlim(2,7) +
  ylim(2,10) + 
  labs(title = "Reference")


p2 <- composition %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(na.value = NA, values = c(NA, "red")) + 
  new_scale_color()+
  geom_contour(data = reference, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(na.value = NA, values = c(NA, "gray")) +
  labs(x = "", y = "") +
  xlim(2,7) +
  ylim(2,10) + 
  labs(title = "Composition") + 
  theme(axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=plot_label_size))
#p2

p3 <- size %>% filter(tot != 150) %>% 
  ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(na.value = NA, values = c(NA, "orange")) + 
  new_scale_color()+
  geom_contour(data = reference, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(na.value = NA, values = c(NA, "gray")) +
  labs(x = "", y = "") +
  xlim(2,7) +
  ylim(2,10) + 
  labs(title = "Abundance") + 
  theme(plot.title = element_text(size=plot_label_size))
#p3

p4 <- halflife %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  scale_colour_manual(na.value = NA, values = c(NA, "forestgreen")) + 
  new_scale_color()+
  geom_contour(data = reference, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(na.value = NA, values = c(NA, "gray")) +
  labs(x="Relative Environmental  ", y = "                                                                             Relative Contact Transmission Rate", title = "Half-life") +
  xlim(2,7) +
  ylim(2,10) +
  theme_classic() +
  theme(axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=plot_label_size),
        axis.title.x=element_text(hjust=1.05))
p4

p5 <- combined %>%
  ggplot(.,aes(x=prop3,y=prop1)) +
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE) +
  scale_colour_manual(na.value = NA, values = c(NA, "purple")) + 
  new_scale_color()+
  geom_contour(data = reference, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F) +
  scale_colour_manual(na.value = NA, values = c(NA, "gray")) +
  labs(x = "                           Transmission Rate", y = "", title = "All") +
  xlim(2,7) +
  ylim(2,10) +
  theme_classic() +
  theme(axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=plot_label_size),
        axis.title.x=element_text(hjust=-0.625))
#p5

figure_1 <- (p2 | p3)/(p4 | p5)
figure_1
ggsave("fig_1.png",plot=figure_1+ proj_theme,width = outwidth[1], height = outwidth[1]/golden,device="png",path=here("figures"))









