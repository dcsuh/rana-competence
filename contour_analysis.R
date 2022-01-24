##figure out whether changes in "combined" treatment are additive of all other treaments or greater than that

##sum values where x and y intercept = 1 and compare to combined
##compare proportional change due to each individual treatment and compare to combined


#define function contour(): inputs are parameters and output is values for grid search for R0
contour <- function(trans_1_min, trans_1_max, trans_3_min, trans_3_max, mort1, mort2, degr) {
  #expand data
  data <- expand_grid("trans1" = seq(from = trans_1_min, to = trans_1_max, by = 0.0000001),
                      "trans2" = c(0.0001),
                      "trans3" = seq(from = trans_3_min, to = trans_3_max, by = 0.0000001),
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

# data$prop1 <- data$trans1/data$trans2 #proportion of transmission of species 2 over species 1
# data$prop3 <- data$trans3/data$trans2 #proportion of environmental transmission over species 1
#"trans2" = c(0.0001),
reference <- contour(trans_1_min = 0.001,trans_1_max = 0.002,trans_3_min = 0.0001,trans_3_max = 0.0001,mort1 = 1/45,mort2 = 1/45,degr = 1/1.947799)
reference_1 <- contour(trans_1_min = 0.0001,trans_1_max = 0.0001,trans_3_min = 0.0001,trans_3_max = 0.0008,mort1 = 1/45,mort2 = 1/45,degr = 1/1.947799)

composition <- contour(trans_1_min = 0.0001,trans_1_max = 0.002,trans_3_min = 0.0001,trans_3_max = 0.0001,mort1 = 1/60,mort2 = 1/30,degr = 1/1.947799)
composition_1 <- contour(trans_1_min = 0.0001,trans_1_max = 0.0001,trans_3_min = 0.0001,trans_3_max = 0.0008,mort1 = 1/60,mort2 = 1/30,degr = 1/1.947799)

size <- contour(trans_1_min = 0.001,trans_1_max = 0.002,trans_3_min = 0.0001,trans_3_max = 0.0001,mort1 = c(1/52.5),mort2 = c(1/52.5),degr = 1/1.947799)
size_1 <- contour(trans_1_min = 0.0001,trans_1_max = 0.0001,trans_3_min = 0.0001,trans_3_max = 0.0008,mort1 = c(1/52.5),mort2 = c(1/52.5),degr = 1/1.947799)

halflife <- contour(trans_1_min = 0.001,trans_1_max = 0.002,trans_3_min = 0.0001,trans_3_max = 0.0001,mort1 = c(1/40),mort2 = 1/40,degr = 1/3.895598)
halflife_1 <- contour(trans_1_min = 0.0001,trans_1_max = 0.0001,trans_3_min = 0.0001,trans_3_max = 0.0008,mort1 = c(1/40),mort2 = 1/40,degr = 1/3.895598)

combined <- contour(trans_1_min = 0.0001,trans_1_max = 0.002,trans_3_min = 0.0001,trans_3_max = 0.0001,mort1 = c(1/70),mort2 = c(1/35),degr = 1/3.895598)
combined_1 <- contour(trans_1_min = 0.0001,trans_1_max = 0.0001,trans_3_min = 0.0001,trans_3_max = 0.0008,mort1 = c(1/70),mort2 = c(1/35),degr = 1/3.895598)


ref_prop1 <- reference_1 %>% filter(prop1==1 & eigen > 1)
a1 <- ref_prop1[1,]
a1 %<>% select(trans1:eigen)
a1$ID <- "A1"
ref_prop3 <- reference %>% filter(prop3==1 & eigen > 1)
a2 <- ref_prop3[1,]
a2 %<>% select(trans1:eigen)
a2$ID <- "A2"

comp_prop1 <- composition_1 %>% filter(prop1==1 & eigen > 1)
b1 <- comp_prop1[1,]
b1 %<>% select(trans1:eigen)
b1$ID <- "B1"
comp_prop3 <- composition %>% filter(prop3==1 & eigen > 1)
b2 <- comp_prop3[1,]
b2 %<>% select(trans1:eigen)
b2$ID <- "B2"

size_prop1 <- size_1 %>% filter(prop1==1 & eigen > 1)
c1 <- size_prop1[1,]
c1 %<>% select(trans1:eigen)
c1$ID <- "C1"
size_prop3 <- size %>% filter(prop3==1 & eigen > 1)
c2 <- size_prop3[1,]
c2 %<>% select(trans1:eigen)
c2$ID <- "C2"

half_prop1 <- halflife_1 %>% filter(prop1==1 & eigen > 1)
d1 <- half_prop1[1,]
d1 %<>% select(trans1:eigen)
d1$ID <- "D1"
half_prop3 <- halflife %>% filter(prop3==1 & eigen > 1)
d2 <- half_prop3[1,]
d2 %<>% select(trans1:eigen)
d2$ID <- "D2"

all_prop1 <- combined_1 %>% filter(prop1==1 & eigen > 1)
e1 <- all_prop1[1,]
e1 %<>% select(trans1:eigen)
e1$ID <- "E1"
all_prop3 <- combined %>% filter(prop3==1 & eigen > 1)
e2 <- all_prop3[1,]
e2 %<>% select(trans1:eigen)
e2$ID <- "E2"

intercepts <- rbind(a1,a2,b1,b2,c1,c2,d1,d2,e1,e2)

#if true, combined approach has a greater change in parameter space compared to additive effect of rest for contact
#when env. trans is very low, how well is contact trans doing?
(a2$prop1-e2$prop1) > (a2$prop1-b2$prop1) + (a2$prop1-c2$prop1) + (a2$prop1-d2$prop1)
(a2$prop1-e2$prop1) - ((a2$prop1-b2$prop1) + (a2$prop1-c2$prop1) + (a2$prop1-d2$prop1))

#if true, combined approach has a greater change in parameter space compared to additive effect of rest for environmental
#if false, the combined effect is actually less than the additive effect
#when contact trans is low., how well is env. trans doing?
(a1$prop3-e1$prop3) > (a1$prop3-b1$prop3) + (a1$prop3-c1$prop3) + (a1$prop3-d1$prop3)
(a1$prop3-e1$prop3) - ((a1$prop3-b1$prop3) + (a1$prop3-c1$prop3) + (a1$prop3-d1$prop3))


a2$prop1-d2$prop1
#interesting thing is that the intercept for halflife is actually higher than the intercept for the reference
#this means that when halflife is extended, favoring the parasite, that contact transmission rate has to be higher than the contact transmission rate for the referecne


tmp <- reference %>% ggplot(.,aes(x=prop3,y=prop1))+
  geom_contour(aes(x=prop3,y=prop1,z=eigen,colour = factor(..level.. == 1,levels = c(F,T)),group=mort1), show.legend = FALSE)+
  scale_colour_manual(na.value = NA, values = c(NA, "red")) + 
  labs(x = "Environmental", y = "Contact") +
  labs(title = "Reference")

tmp

intercept <- reference %>% filter(prop1==1 & eigen > 1)
intercept_1 <- reference %>% filter(prop3==1 & eigen > 1)

intercept_2 <- reference %>% filter(prop3==1 & prop1==1)
