### Daniel Suh


### stats tests

library(here)

source(knitr::purl(here("scripts/data_format.Rmd"), quiet=TRUE))


#### multivariate test for effect of cc, abundance, and temp on future month's infection prevalence

tmp <- data %>% 
  group_by(WetAltID, Month.1) %>% 
  summarize(infected = sum(RV.Status), total = n()) 
tmp$Month <-gsub("Month", "", tmp$Month.1)
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


### univariate version of multivariate glm

Model_cc <- glm(as.matrix(M)~lagged$lag_cc, family="quasibinomial")
summary(Model_cc)
Model_size <- glm(as.matrix(M)~lagged$lag_size, family="quasibinomial")
summary(Model_size)
Model_temp <- glm(as.matrix(M)~lagged$lag_temp, family="quasibinomial")
summary(Model_temp)






#### regression analysis for cc~size for each month

corr_plots <- comm_summ %>% ggplot(., aes(x = log10(size), y = cc)) + geom_point() + geom_smooth(method="lm") + facet_wrap(vars(Month.1))

by_month <- comm_summ %>% group_by(Month.1) %>% nest()

p_fun <- function(df){
  return((cor.test(log10(df$size), df$cc, method = "spearman")$p.value))
}
est_fun <- function(df){
  return((cor.test(log10(df$size), df$cc, method = "spearman")$estimate))
}

by_month <- mutate(by_month, p.val = purrr::map(data, cor_fun), estimate = purrr::map(data, est_fun))


#### geo-additive gam for cc~space,month,size
#### space found to not be significant so instead just plotted regressions for each month (seen above)
#### this can go into supp methods

library(mgcv)
library(biogeo)

coords <- read_csv(here("data/srs_wetland_coordinates.csv"))
coords$dd_x <- dms2dd(coords$xdd,coords$xmm,coords$xss,coords$xdir)
coords$dd_y <- dms2dd(coords$ydd,coords$ymm,coords$yss,coords$ydir)
coords %<>% dplyr::select(WetAltID, dd_x, dd_y)

geo_gam <- full_join(comm_summ, coords, by = "WetAltID") %>% dplyr::select(siteID, WetAltID, Month.1, cc, size, dd_x, dd_y)
geo_gam$Month <- gsub("Month", "", geo_gam$Month.1)
geo_gam$Month <- as.factor(geo_gam$Month)

gam1<-gam(cc~te(dd_x,dd_y)+s(Month,bs="re")+s(log(size)),family=Gamma(link=log),data=geo_gam)
gam2<-gam(cc~s(Month,bs="re")+s(log(size)),family=Gamma(link=log),data=geo_gam) #no geo-additive term

summary(gam1)
plot.gam(gam1)



