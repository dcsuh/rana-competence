#Daniel Suh and Andrew Park
#6/22/23

#Figure 1a in rv_cc manuscript


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
  data$half <- log(2)/data$degr #viral half-life
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



rm0<-function(t,x,params){
  Sa <- x[1]  
  Ia <- x[2]
  Ra <- x[3]
  Sb <- x[4]  
  Ib <- x[5]
  Rb <- x[6]
  V <- x[7]
  with(as.list(params),{
    dSa <- -beta*Sa*Ia - betab*Sa*Ib - phi*Sa*V + lambda - mu_a*Sa
    dIa <- beta*Sa*Ia + betab*Sa*Ib + phi*Sa*V - alpha*Ia - mu_a*Ia
    dRa <- alpha*Ia - mu_a*Ra
    
    dSb <- -betab*Sb*Ib - beta*Sb*Ia - phi*Sb*V + lambda - mu_b*Sb
    dIb <- betab*Sb*Ib + beta*Sb*Ia + phi*Sb*V - alpha*Ib - mu_b*Ib
    dRb <- alpha*Ib - mu_b*Rb
    
    dV <- sigma*Ia + sigma*Ib - epsilon*V
    res<-c(dSa,dIa,dRa,dSb,dIb,dRb,dV)
    list(res)
  })}
#no demography
rm1<-function(t,x,params){
  Sa <- x[1]  
  Ia <- x[2]
  Ra <- x[3]
  Sb <- x[4]  
  Ib <- x[5]
  Rb <- x[6]
  V <- x[7]
  with(as.list(params),{
    dSa <- -beta*Sa*Ia - betab*Sa*Ib - phi*Sa*V
    dIa <- beta*Sa*Ia + betab*Sa*Ib + phi*Sa*V - alpha*Ia
    dRa <- alpha*Ia
    
    dSb <- -betab*Sb*Ib - beta*Sb*Ia - phi*Sb*V
    dIb <- betab*Sb*Ib + beta*Sb*Ia + phi*Sb*V - alpha*Ib
    dRb <- alpha*Ib
    
    dV <- sigma*Ia + sigma*Ib - epsilon*V
    res<-c(dSa,dIa,dRa,dSb,dIb,dRb,dV)
    list(res)
  })}


maxTime <- 300.0 # time is in years - run model for this time
times<-seq(0,maxTime,by=1) # how long we run the model for
# notes on params
# mu =     <- density-independent mortality rate
# beta =       <- contact transmission rate by infected A
# betab =       <- contact transmission rate by infected B
# phi =      <- environmental transmission rate
# lambda =       <- density-dependent growth rate
# alpha =     <- recovery rate
# sigma =    <- shedding rate
# epsilon =       <- viral degradation rate
params<-c(mu_a=1/50,
          mu_b=1/50,
          beta=0.001,
          betab=0.0001,
          phi=0.001,
          lambda=5/3,
          alpha=1/10,
          sigma=1/2,
          epsilon=1/2)  # model parameters

xstart<-c(Sa=116,
          Ia=1,
          Ra=0,
          Sb=83,
          Ib=1,
          Rb=0,
          V=0)  # initial conditions

output<-as.data.frame(lsoda(xstart,times,rm0,params)) # tells computer to solve (integrate) equations


spec_2 <- 0.0001
spec_1 <- 0.00055
env_rate <- 0.00065
param_num <- which(reference$eigen>1.103065 & reference$eigen<1.103192)


params_n <- c(mu_a=1/45,
              mu_b=1/45,
              beta=spec_1,
              betab=spec_2,
              phi=env_rate,
              lambda=5/3,
              alpha=1/10,
              sigma=1/2,
              epsilon=1/1.947799)  # model parameters

xstart_n <- c(Sa=params_n[[6]]/params_n[[1]],
              Ia=1,
              Ra=0,
              Sb=params_n[[6]]/params_n[[2]],
              Ib=1,
              Rb=0,
              V=0)  # initial conditions

output_n <- as.data.frame(lsoda(xstart_n,times,rm0,params_n)) # tells computer to solve (integrate) equations
output_n1<- as.data.frame(lsoda(xstart_n,times,rm1,params_n)) # tells computer to solve (integrate) equations

reference_eigen <- reference$eigen[param_num]


params_c <- c(mu_a=1/60,
              mu_b=1/30,
              beta=spec_1,
              betab=spec_2,
              phi=env_rate,
              lambda=5/3,
              alpha=1/10,
              sigma=1/2,
              epsilon=1/1.947799)  # model parameters

xstart_c <- c(Sa=params_c[[6]]/params_c[[1]],
              Ia=1,
              Ra=0,
              Sb=params_c[[6]]/params_c[[2]],
              Ib=1,
              Rb=0,
              V=0)  # initial conditions

output_c <- as.data.frame(lsoda(xstart_c,times,rm0,params_c)) # tells computer to solve (integrate) equations
output_c1 <- as.data.frame(lsoda(xstart_c,times,rm1,params_c)) # tells computer to solve (integrate) equations

community_eigen <- composition$eigen[param_num]



params_a <- c(mu_a=1/52.5,
              mu_b=1/52.5,
              beta=spec_1,
              betab=spec_2,
              phi=env_rate,
              lambda=5/3,
              alpha=1/10,
              sigma=1/2,
              epsilon=1/1.947799)  # model parameters

xstart_a <- c(Sa=params_a[[6]]/params_a[[1]],
              Ia=1,
              Ra=0,
              Sb=params_a[[6]]/params_a[[2]],
              Ib=1,
              Rb=0,
              V=0)  # initial conditions

output_a <- as.data.frame(lsoda(xstart_a,times,rm0,params_a)) # tells computer to solve (integrate) equations
output_a1 <- as.data.frame(lsoda(xstart_a,times,rm1,params_a)) # tells computer to solve (integrate) equations

abundance_eigen <- size$eigen[param_num]



params_h <- c(mu_a=1/45,
              mu_b=1/45,
              beta=spec_1,
              betab=spec_2,
              phi=env_rate,
              lambda=5/3,
              alpha=1/10,
              sigma=1/2,
              epsilon=1/3.895598)  # model parameters

xstart_h <- c(Sa=params_h[[6]]/params_h[[1]],
              Ia=1,
              Ra=0,
              Sb=params_h[[6]]/params_h[[2]],
              Ib=1,
              Rb=0,
              V=0)  # initial conditions

output_h <- as.data.frame(lsoda(xstart_h,times,rm0,params_h)) # tells computer to solve (integrate) equations
output_h1 <- as.data.frame(lsoda(xstart_h,times,rm1,params_h)) # tells computer to solve (integrate) equations


halflife_eigen <- halflife$eigen[param_num]



params_x <- c(mu_a=1/70,
              mu_b=1/35,
              beta=spec_1,
              betab=spec_2,
              phi=env_rate,
              lambda=5/3,
              alpha=1/10,
              sigma=1/2,
              epsilon=1/3.895598)  # model parameters



xstart_x <- c(Sa=params_x[[6]]/params_x[[1]],
              Ia=1,
              Ra=0,
              Sb=params_x[[6]]/params_x[[2]],
              Ib=1,
              Rb=0,
              V=0)  # initial conditions



output_x1 <- as.data.frame(lsoda(xstart_x,times,rm1,params_x)) # tells computer to solve (integrate) equations



combined_eigen <- combined$eigen[param_num]

output_a1 %<>% mutate(I = Ia+Ib)
output_c1 %<>% mutate(I = Ia+Ib)
output_h1 %<>% mutate(I = Ia+Ib)
output_n1 %<>% mutate(I = Ia+Ib)
output_x1 %<>% mutate(I = Ia+Ib)

output_a1 %<>% mutate(X = I + V)
output_c1 %<>% mutate(X = I + V)
output_h1 %<>% mutate(X = I + V)
output_n1 %<>% mutate(X = I + V)
output_x1 %<>% mutate(X = I + V)


color_x <- paste("Combined\nR0=", abbreviate(combined_eigen))

color_c <- paste("Composition\nR0=", abbreviate(community_eigen))
color_a <- paste("Abundance\nR0=", abbreviate(abundance_eigen))
color_h <- paste("Halflife\nR0=", abbreviate(halflife_eigen))
color_n <- paste("Reference\nR0=", abbreviate(reference_eigen))

both <- ggplot()+
  geom_line(data=output_x1,mapping = aes(y=I,x=time, color=color_x))+
  scale_color_manual(name="",values = "purple")+
  new_scale_color()+
  geom_line(data=output_c1, mapping=aes(y=I,x=time, color=color_c))+
  scale_color_manual(name="",values = "red")+
  new_scale_color()+
  geom_line(output_a1, mapping=aes(y=I,x=time, color=color_a))+
  scale_color_manual(name="",values = "orange")+
  new_scale_color()+
  geom_line(output_h1, mapping=aes(y=I,x=time, color=color_h))+
  scale_color_manual(name="",values = "forestgreen")+
  new_scale_color()+
  geom_line(output_n1, mapping=aes(y=I,x=time, color=color_n))+
  scale_color_manual(name="",values = "gray")+
  labs(y="Number of Infected Individuals",x="")+
  ylim(0,45)+
  ggtitle(paste("Both"))+
  theme_bw() + 
  theme(legend.position = "none")

spec_A <- ggplot()+
  geom_line(data=output_x1,mapping = aes(y=Ia,x=time, color=color_x))+
  scale_color_manual(name="",values = "purple")+
  new_scale_color()+
  geom_line(data=output_c1, mapping=aes(y=Ia,x=time, color=color_c))+
  scale_color_manual(name="",values = "red")+
  new_scale_color()+
  geom_line(output_a1, mapping=aes(y=Ia,x=time, color=color_a))+
  scale_color_manual(name="",values = "orange")+
  new_scale_color()+
  geom_line(output_h1, mapping=aes(y=Ia,x=time, color=color_h))+
  scale_color_manual(name="",values = "forestgreen")+
  new_scale_color()+
  geom_line(output_n1, mapping=aes(y=Ia,x=time, color=color_n))+
  scale_color_manual(name="",values = "gray")+
  labs(y="",x="Time")+
  ylim(0,45)+
  ggtitle(paste("High competence"))+
  theme_bw()+
  theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())

spec_B <- ggplot()+
  geom_line(data=output_x1,mapping = aes(y=Ib,x=time, color=color_x))+
  scale_color_manual(name="",values = "purple")+
  new_scale_color()+
  geom_line(data=output_c1, mapping=aes(y=Ib,x=time, color=color_c))+
  scale_color_manual(name="",values = "red")+
  new_scale_color()+
  geom_line(output_a1, mapping=aes(y=Ib,x=time, color=color_a))+
  scale_color_manual(name="",values = "orange")+
  new_scale_color()+
  geom_line(output_h1, mapping=aes(y=Ib,x=time, color=color_h))+
  scale_color_manual(name="",values = "forestgreen")+
  new_scale_color()+
  geom_line(output_n1, mapping=aes(y=Ib,x=time, color=color_n))+
  scale_color_manual(name="",values = "gray")+
  labs(y="",x="")+
  ylim(0,45)+
  ggtitle(paste("Low competence"))+
  theme_bw()+
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

dynamics <- both + spec_A + spec_B



ggsave("fig_1a.png",plot=dynamics,width = outwidth[1], height = outwidth[1]/golden,device="png",path=here("figures"))

