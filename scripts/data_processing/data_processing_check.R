library(here)


old <- read_csv("data/weighted_prev_competence.csv")
new <- read_csv("data/weighted_prev_competence_111220.csv")

vl <- new %>% select(vl.4:vl.42) %>% distinct() %>% pivot_longer(vl.4:vl.42)

old %<>% select(WetAltID, Month.1, RA2:RA42, cc) %>% distinct() 
new %<>% select(WetAltID, Month.1, RA2:RA42, cc) %>% distinct()

old_12_4 <- old %>% filter(WetAltID == 12 & Month.1 == "Month4") %>% pivot_longer(RA2:RA42) %>% select(name, value)
new_12_4 <- new %>% filter(WetAltID == 12 & Month.1 == "Month4") %>% pivot_longer(RA2:RA42) %>% select(name, value)

old %<>% select(WetAltID, Month.1, cc) %>% mutate(old_cc = cc, siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() %>% select(-c(WetAltID, Month.1, cc))
new %<>% select(WetAltID, Month.1, cc) %>% mutate(new_cc = cc, siteID = paste(WetAltID, Month.1, sep = "_")) %>% select(-c(WetAltID, Month.1, cc))

join <- full_join(old, new, by = "siteID")
join %<>% mutate(old_minus_new = old_cc-new_cc)
join %>% ggplot(., aes(x=old_cc, y=new_cc)) + geom_point()
join %>% ggplot(.,aes(x=siteID,y=old_minus_new)) + geom_point() + theme(axis.text.x = element_text(angle=90))



#largest discrepancy is site 12 month 4 and then the 5 site-months that had na for cc because of na's in old RA values
join_12_4 <- full_join(old_12_4, new_12_4, by = "name")
join_12_4 %<>% mutate(.,old_minus_new = value.x-value.y)
join_12_4 %>% ggplot(.,aes(x=name,y=old_minus_new)) + geom_point() + theme(axis.text.x = element_text(angle=90))
#old 12_4 is much higher because it has a relative abundance of species 21 that is 93%
#new 12_4 is diluted by species 5, comprising 25% of the community. Species 5 has no vl so it counts as 0
