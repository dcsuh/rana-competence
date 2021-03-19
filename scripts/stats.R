### Daniel Suh


### stats tests

library(here)

source(knitr::purl(here("/scripts/data_format.Rmd"), quiet=TRUE))


tmp <- data %>% 
  group_by(WetAltID, Month.1) %>% 
  summarize(infected = sum(RV.Status), total = n()) 
tmp$Month <-   gsub("Month", "", tmp$Month.1)
tmp %<>%
  add_column(., prev = tmp$infected/tmp$total) %>%
  mutate(., siteID = paste(WetAltID, Month, sep = "_")) 

M <- tmp[which(tmp$siteID %in% clean$siteID),]
M %<>% ungroup() %>% select(prev, total, siteID)

lagged <- clean %>% ungroup() %>% select(siteID, lag_cc, lag_size, lag_temp) %>% drop_na()

M <- M[which(M$siteID %in% lagged$siteID),]

M <- M[order(M$siteID),]
M %<>% select(-siteID)
lagged <- lagged[order(lagged$siteID),]

Model <- glm(as.matrix(M)~lagged$lag_cc+lagged$lag_size+lagged$lag_temp, family="quasibinomial")

summary(Model)


library(mgcv)
library(biogeo)

coords <- read_csv(here("/data/srs_wetland_coordinates.csv"))
coords$dd_x <- dms2dd(coords$xdd,coords$xmm,coords$xss,coords$xdir)
coords$dd_y <- dms2dd(coords$ydd,coords$ymm,coords$yss,coords$ydir)



gam1<-gam(cc~te(x,y)+s(pop2016)+s(mean)+s(hi)+s(statefips,bs="re"),family="binomial",data=comm_summ)
summary(gam1)