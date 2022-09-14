#Daniel Suh and Andrew Park
#1/24/22

#another version of figure 1 from manuscript

library(tidyverse)
library(magrittr)
#install_github("thomasp85/patchwork")
library(patchwork)
library(ggnewscale)
library(here)

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

reference <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0,trans_3_max = 0.001,
                     mort1 = 1/45,mort2 = 1/45,degr = 1/1.947799)

composition <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0,trans_3_max = 0.001,
                       mort1 = 1/60,mort2 = 1/30,degr = 1/1.947799)

size <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0,trans_3_max = 0.001,
                mort1 = c(1/52.5),mort2 = c(1/52.5),degr = 1/1.947799)

halflife <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0,trans_3_max = 0.001,
                    mort1 = c(1/40),mort2 = 1/40,degr = 1/3.895598)

combined <- contour(trans_1_min = 0.0001,trans_1_max = 0.001,trans_3_min = 0,trans_3_max = 0.001,
                    mort1 = c(1/70),mort2 = c(1/35),degr = 1/3.895598)


axis_text_size = 8
plot_label_size = 20

p1 <- reference %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(na.value = NA, values = c(NA, "red")) + 
  labs(x = "Environmental", y = "Contact") +
  xlim(1,8) +
  ylim(1,14) + 
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
  xlim(1,8) +
  ylim(1,14) + 
  labs(title = "A") + 
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
  xlim(1,8) +
  ylim(1,14) + 
  labs(title = "B") + 
  theme(plot.title = element_text(size=plot_label_size))
#p3

p4 <- halflife %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  scale_colour_manual(na.value = NA, values = c(NA, "forestgreen")) + 
  new_scale_color()+
  geom_contour(data = reference, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(na.value = NA, values = c(NA, "gray")) +
  labs(x="Relative Environmental", y = "                                                               Relative Contact Transmission Rate", title = "C") +
  xlim(1,8) +
  ylim(1,14) + 
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
  labs(x = "Transmission Rate", y = "", title = "D") +
  xlim(1,8) +
  ylim(1,14) + 
  theme_classic() +
  theme(axis.title = element_text(size=axis_text_size),
        plot.title = element_text(size=plot_label_size),
        axis.title.x=element_text(hjust=-0.625))
#p5

figure_1 <- (p2 | p3)/(p4 | p5)
figure_1

fig_2 <- reference %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  theme_classic()+
  scale_colour_manual(na.value = NA, values = c(NA, "gray")) + 
  new_scale_color()+
  geom_contour(data = composition, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(na.value = NA, values = c(NA, "red")) +
  new_scale_color()+
  geom_contour(data = size, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(na.value = NA, values = c(NA, "orange")) +
  new_scale_color()+
  geom_contour(data = halflife, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(na.value = NA, values = c(NA, "forestgreen")) +
  new_scale_color()+
  geom_contour(data = combined, 
               aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = F)+
  scale_colour_manual(na.value = NA, values = c(NA, "purple")) +
  geom_vline(xintercept = 1)+
  geom_hline(yintercept = 1)+
  labs(x = "env. trans : low beta", y = "high beta : low beta") +
  xlim(1,8) +
  ylim(1,14) + 
  labs(title = "gray=reference; red=composition; orange=abundance; green=halflife; purple=combined")

reference %>% ggplot(.,aes(x=prop3,y=prop1,fill=eigen)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red")
composition %>% ggplot(.,aes(x=prop3,y=prop1,fill=eigen)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red")
size %>% ggplot(.,aes(x=prop3,y=prop1,fill=eigen)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red")
halflife %>% ggplot(.,aes(x=prop3,y=prop1,fill=eigen)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red")
combined %>% ggplot(.,aes(x=prop3,y=prop1,fill=eigen)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red")


#heatmaps and contour plots of difference between control and treatment

ref_sel <- reference %>% select(prop3,prop1,eigen) %>% rename(.,ref_eigen = eigen)
comp_sel <- composition %>% select(prop3,prop1,eigen) %>% rename(.,comp_eigen = eigen)
size_sel <- size %>% select(prop3,prop1,eigen) %>% rename(.,size_eigen = eigen)
halflife_sel <- halflife %>% select(prop3,prop1,eigen) %>% rename(.,half_eigen = eigen)
combined_sel <- combined %>% select(prop3,prop1,eigen) %>% rename(.,all_eigen = eigen)


ref_comp <- inner_join(ref_sel,comp_sel, by = c("prop3","prop1"))
ref_comp %<>% inner_join(.,size_sel, by =c("prop3","prop1"))
ref_comp %<>% inner_join(.,halflife_sel, by =c("prop3","prop1"))
ref_comp %<>% inner_join(.,combined_sel, by =c("prop3","prop1"))


ref_comp %<>% mutate(.,comp_diff = comp_eigen-ref_eigen,
                     size_diff = size_eigen-ref_eigen,
                     half_diff = half_eigen-ref_eigen,
                     all_diff = all_eigen-ref_eigen)

#heatmaps with contours

comp_fig <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1,fill=comp_diff)) + 
  geom_tile() +
  scale_fill_gradient2() + 
  geom_contour(aes(z=comp_diff))

size_fig <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1,fill=size_diff)) + 
  geom_tile() +
  scale_fill_gradient2() + 
  geom_contour(aes(z=size_diff))

half_fig <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1,fill=half_diff)) + 
  geom_tile() +
  scale_fill_gradient2() + 
  geom_contour(aes(z=half_diff))

combined_fig <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1,fill=all_diff)) + 
  geom_tile() +
  scale_fill_gradient2() + 
  geom_contour(aes(z=all_diff))

(comp_fig | size_fig)/(half_fig | combined_fig)

#filled contour with standardized breaks

breaks <- c(-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,2)

comp_fig <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=comp_diff),breaks = breaks) +
  scale_fill_viridis_d(drop=F) +
  theme(legend.position = "none") + 
  labs(x="environmental transmission", y = "contact transmission", title = "composition")

size_fig <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=size_diff),breaks = breaks) +
  scale_fill_viridis_d(drop=F) +
  theme(legend.position = "none") + 
  labs(x="environmental transmission", y = "contact transmission", title = "abundance")

half_fig <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=half_diff),breaks = breaks) +
  scale_fill_viridis_d(drop=F) +
  theme(legend.position = "none") + 
  labs(x="environmental transmission", y = "contact transmission", title = "halflife")

combined_fig <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=all_diff),breaks = breaks) +
  scale_fill_viridis_d(drop=F) +   
  labs(x="environmental transmission", y = "contact transmission", title = "combined")

(comp_fig | size_fig)/(half_fig | combined_fig)

#comparison between combined treatment and the additive effect

ref_comp %<>% rowwise %>% mutate(diff_sum = comp_diff + size_diff + half_diff)
ref_comp %<>% rowwise %>% mutate(effect = all_diff-diff_sum)

ref_comp %>% ggplot(.,aes(x=prop3,y=prop1,fill=effect)) + 
  geom_tile() +
  scale_fill_gradient2() + geom_contour(aes(z=effect))

ref_comp %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=effect))

breaks <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

comb <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=all_diff),breaks = breaks) +
  scale_fill_viridis_d(drop=F) +   
  theme(legend.position = "none") + 
  labs(x="", y = "", title = "combined")

sum <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=diff_sum),breaks = breaks) +
  scale_fill_viridis_d(drop=F) + 
  labs(x = "", y = "contact transmission", title = "sum of its parts")

diff <- ref_comp %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=effect),breaks = breaks) +
  scale_fill_viridis_d(drop=F) + 
  theme(legend.position = "none") +
  labs(x = "environmental", y = "", title = "difference (combined - sum of its parts)")

comb / sum / diff

#comparison between combined treament and additive effect by percent change

ref_percent <- ref_comp %>% mutate(., comp_per = comp_diff/ref_eigen,
                                   size_per = size_diff/ref_eigen,
                                   half_per = half_diff/ref_eigen,
                                   all_per = all_eigen/ref_eigen,
                                   sum_per = diff_sum/ref_eigen)

breaks <- c(0:20)*0.1

comb_per <- ref_percent %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=all_per),breaks = breaks) +
  scale_fill_viridis_d(drop=F) +   
  theme(legend.position = "none") + 
  labs(x="", y = "", title = "combined")

sum_per <- ref_percent %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=sum_per),breaks = breaks) +
  scale_fill_viridis_d(drop=F) + 
  labs(x = "", y = "contact transmission", title = "sum of its parts")

diff_per <- ref_percent %>% ggplot(.,aes(x=prop3,y=prop1)) + 
  geom_contour_filled(aes(z=all_per-sum_per),breaks = breaks) +
  scale_fill_viridis_d(drop=F) + 
  theme(legend.position = "none") +
  labs(x = "environmental", y = "", title = "difference (combined - sum of its parts)")

comb_per / sum_per / diff_per


