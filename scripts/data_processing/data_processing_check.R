library(here)


old <- read_csv("data/weighted_prev_competence.csv")
new <- read_csv("data/weighted_prev_competence_111220.csv")

old %<>% select(WetAltID, Month.1, RA2:RA42, cc) %>% distinct() 
new %<>% select(WetAltID, Month.1, RA2:RA42, cc) %>% distinct()


old %<>% select(WetAltID, Month.1, cc) %>% mutate(old_cc = cc, siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() %>% select(-c(WetAltID, Month.1, cc))
new %<>% select(WetAltID, Month.1, cc) %>% mutate(new_cc = cc, siteID = paste(WetAltID, Month.1, sep = "_")) %>% select(-c(WetAltID, Month.1, cc))

join <- full_join(old, new, by = "siteID")
join %<>% mutate(diff = old_cc-new_cc)
join %>% ggplot(., aes(x=old_cc, y=new_cc)) + geom_point()

