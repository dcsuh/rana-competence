
library(here)

source(here("base", "src.R"))
library(deSolve)

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

dat <- reference %>% dplyr::select(prop3, prop1, eigen) %>% mutate(reference = eigen) %>% dplyr::select(-eigen)

get_invade <- function(df1, df2){
  df2 %<>% dplyr::select(prop3, prop1, eigen) %>% mutate(manip = eigen) %>% dplyr::select(-eigen)
  df1 %<>% left_join(., df2)
  
  output <- df1 %>% mutate(invade = case_when(reference > 1 & manip > 1 ~ "yes",
                                              reference < 1 & manip > 1 ~ "maybe",
                                              reference < 1 & manip < 1 ~ "no"))
  
  return(output)
}

get_plot <- function(df){
  df %>% ggplot(., aes(x=prop3, y=prop1, fill = invade)) + 
    geom_tile(show.legend = F) + 
    scale_fill_manual(values = c("orange", "red", "forestgreen")) + 
    xlim(2, 7) + 
    ylim(2, 7) + 
    theme_classic()
}

dat_comp <- get_invade(dat, composition)
dat_size <- get_invade(dat, size)
dat_half <- get_invade(dat, halflife)
dat_all <- get_invade(dat, combined)

p1 <- get_plot(dat_comp) + labs(title = "composition", x = "", y = "")
p2 <- get_plot(dat_size) + labs(title = "abundance", x = "", y = "")
p3 <- get_plot(dat_half) + labs(title = "half-life", x = "", y = "")
p4 <- get_plot(dat_all) + labs(title = "combined", x = "", y = "")

p1 / p2 / p3 / p4
