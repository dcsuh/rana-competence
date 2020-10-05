#Community competence and Diversity

#Species richness
#Pielou's J for species evenness


library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
library(patchwork)
library(ggnewscale)


data <- read_csv("prev_test.csv")

#Get species richness from matrix

community_mat <- data %>% dplyr::select(WetAltID, Month.1, AB2:AB8, AB9, AB20:AB42) %>%
  mutate(., siteID = paste(WetAltID, Month.1, sep = "_")) %>% distinct() 

community_mat <- community_mat[order(community_mat$WetAltID, community_mat$Month.1),]

community_mat %<>% arrange(.,WetAltID) %>% dplyr::select(-WetAltID, -Month.1) %>% dplyr::select(siteID, AB2:AB42)

abundances <- community_mat

community_mat %<>% column_to_rownames(., var = "siteID")

presence_mat <- ifelse(community_mat>0, 1, 0)

siteID <- rownames(presence_mat)

richness_mat <- tibble()
for (i in 1:length(siteID)){
  richness_mat[i,1] <- siteID[i]
  richness_mat[i,2] <- sum(presence_mat[i,])
}
richness_mat %<>% rename(siteID = 1, richness = 2)

order <- richness_mat

#Get species richness for each site and community competence values

tmp <- data %>% dplyr::select(WetAltID, Month.1, cc, Prevalence,AB2:AB8, AB9, AB20:AB42) %>% 
    distinct() %>% 
    mutate(., siteID = paste(WetAltID, Month.1, sep = "_"))
tmp %<>% mutate(.,size = rowSums(.[5:25])) %>% select(-c(AB2:AB42))
tmp <- tmp[order(tmp$WetAltID, tmp$Month.1),]
tmp %<>% dplyr::select(-Month.1)

richness_cc <- full_join(tmp, richness_mat, by = "siteID")



#J = H'/ln(S)
#H' = Shannon diversity index
#S = Species Richness

h <- vegan::diversity(community_mat, index = "shannon")
h <- melt(as.matrix(h))
h %<>% rename(siteID = Var1, h = value) %>% dplyr::select(-Var2)

evenness <- richness_cc %>% mutate(lnS = log(richness)) %>% rowwise()

evenness %<>% full_join(h, evenness, by = "siteID")

evenness %<>% mutate(J = h/lnS) %>% rowwise()

#create new df for ordered sites to test lagged richness and evenness on prevalence
lag_evenness <- evenness[-c(4,14,47,90),]
lag_evenness$lag_richness <- c(0)
lag_evenness$lag_J <- c(0)
lag_evenness$lag_cc <- c(0)
lag_evenness$lag_size <- c(0)
for(n in 2:91){
  lag_evenness$lag_richness[n] <- lag_evenness$richness[n-1]
  lag_evenness$lag_J[n] <- lag_evenness$J[n-1]
  lag_evenness$lag_cc[n] <- lag_evenness$cc[n-1]
  lag_evenness$lag_size[n] <- lag_evenness$size[n-1]
}
#remove the first entry for each wetland to remove the carryover from the last wetland
lag_evenness %<>% group_by(WetAltID) %>% filter(duplicated(WetAltID) | n()==1)

p1 <- evenness %>% remove_missing() %>% ggplot(aes(x=richness, y=cc)) + 
  geom_smooth(method ="lm") +
  geom_point() +
  labs(x = "Richness", y = "CC", title = ) +
  theme_classic()

cor.test(richness_cc$cc, richness_cc$richness, method="spearman")

p2 <- lag_evenness %>% remove_missing() %>% ggplot(aes(x=lag_richness, y=Prevalence)) + 
  geom_smooth(method ="lm") +
  geom_point() +
  labs(x = "Richness", y = "Prevalence", title = ) +
  theme_classic()

cor.test(x = lag_evenness$lag_richness, y = lag_evenness$Prevalence, method = "spearman")

p3 <- lag_evenness %>% ggplot(aes(x = lag_J, y = lag_cc, color = Prevalence, size = lag_size)) +
  geom_point() +
  labs(x="Evenness", y = "CC", size = "Size") +
  geom_smooth(method = "lm", show.legend = F) +
  scale_color_gradient(low="blue", high="red") + 
  theme_minimal()

cor.test(x=lag_evenness$lag_J, y = lag_evenness$lag_cc, method = "spearman")
  
final <- (p1|p3)/p2
p1
p2
p3
final

#barplots comparing community composition

competence_values <- data %>% select(Species, SQMean) %>% distinct() %>% group_by(Species) %>% summarize(comp = mean(SQMean))

log_abundances <- abundances %>% as_tibble() %>% column_to_rownames(., var = "siteID") %>% na_if(.,0) %>% log10(.) %>% 
  replace(is.na(.), 0) %>% rownames_to_column(., var = "siteID") %>% rename_with(., ~ gsub("AB","logAB",.x,fixed = "TRUE"))
log_abundances %<>% mutate(., logABLow = (logAB2 + logAB3 + logAB4 + logAB5 + logAB6 + logAB8 + 
                                            logAB20 + logAB24 + logAB27 + logAB28 + logAB29 + 
                                            logAB31 + logAB34 + logAB35 + logAB38 + logAB39 + logAB41))
log_abundances %<>% mutate(., logABhigh = (logAB9 + logAB21 + logAB26 + logAB42))
log_abundances %<>% mutate(., logtotal = rowSums(.[2:22]))

abundances <- full_join(abundances,lag_evenness, by = "siteID")
abundances <- mutate(abundances, ABHigh = (AB9 + AB21 + AB26 + AB42)) 
abundances <- mutate(abundances, ABLow = (AB2 + AB3 + AB4 + AB5 + AB6 + AB8 + AB20 + AB24 + AB27 + AB28 + AB29 + AB31 + AB34 + AB35 + AB38 + AB39 + AB41))
abundances <- mutate(abundances, total = ABHigh + ABLow)
log_abundances <- full_join(abundances,log_abundances, by = "siteID")


abundances %>% filter(J < 0.5) %>% dplyr::select(siteID, ABHigh, ABLow) %>% melt() %>%
  ggplot(., aes(y = value, x = siteID, fill = variable)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  labs(title = "Abundance of high and low competence for Pielou's J < 0.5") +
  theme(axis.text.x = element_text(angle=90))

abundances %>% filter(J < 0.5) %>% filter(cc > 45000) %>% dplyr::select(siteID, ABHigh, ABLow) %>% melt() %>%
  ggplot(., aes(y = value, x = siteID, fill = variable)) +
  geom_bar(position = position_dodge(), stat = "identity") + 
  labs(title = "Abundance of high and low competence for Pielou's J < 0.5 and cc > 45000") +
  theme(axis.text.x = element_text(angle=90))

abundances %>% filter(J < 0.5) %>% filter(cc < 45000) %>% dplyr::select(siteID, ABHigh, ABLow) %>% melt() %>%
  ggplot(., aes(y = value, x = siteID, fill = variable)) +
  geom_bar(position = position_dodge(), stat = "identity") + 
  labs(title = "Abundance of high and low competence for Pielou's J < 0.5 and cc < 45000") +
  theme(axis.text.x = element_text(angle=90))

abundances %>% filter (J<0.5) %>% mutate(.,ratio = ABHigh/ABLow)%>%
  ggplot(., aes(y = ratio, x = cc)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  labs(title = "Abundance of high and low competence for Pielou's J < 0.5", y = "ratio(high:low)") +
  theme(axis.text.x = element_text(angle=90))




abundances$siteID <- reorder(abundances$siteID, -abundances$total) #order by total community size
log_abundances$siteID <- reorder(log_abundances$siteID, -log_abundances$total) #order by total community size
#abundances$siteID <- reorder(abundances$siteID, -abundances$cc) #order by community competence
#abundances$siteID <- reorder(abundances$siteID, -abundances$J) #order by Pielou's J
comp_plot_a <- abundances %>% filter(J<0.6)  %>% dplyr::select(siteID, AB26, AB42, AB21, AB9, ABLow) %>% melt() %>%
  ggplot(., aes(y = log10(value), x = siteID, fill = variable)) +
  geom_bar(position ="stack", stat = "identity") +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#ffffb2", "black")) +
  labs(title = "Abundance of high and low competent species under Pielou's J of 0.6") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Abundance") +
  labs(fill = "log10(variable)")

comp_plot_b <- abundances %>% filter(J<0.6)  %>% dplyr::select(siteID, ABHigh, ABLow) %>% melt() %>%
  ggplot(., aes(y = log10(value), x = siteID, fill = variable)) +
  geom_bar(position ="stack", stat = "identity") +
  scale_fill_manual(values = c("#238443","black")) +
  labs(title = "Abundance of high and low competent species under Pielou's J of 0.6") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Abundance") +
  labs(fill = "log10(variable)")

comp_plot_c <- log_abundances %>% filter(J<0.6)  %>% dplyr::select(siteID, logAB26, logAB42, logAB21, logAB9, logABLow) %>% melt() %>%
  ggplot(., aes(y = (value), x = siteID, fill = variable)) +
  geom_bar(position ="stack", stat = "identity") +
  scale_fill_manual(values = c("#238443", "#78c679", "#c2e699", "#ffffb2", "black")) +
  labs(title = "Abundance of high and low competent species under Pielou's J of 0.6") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Abundance")

cc_plot <- abundances %>% dplyr::select(siteID, cc) %>% distinct() %>%
  ggplot(., aes(y=cc, x=siteID)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  labs(title = "Values of cc for sites ordered by community size (descending)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))


evenness_plot <- abundances %>% dplyr::select(siteID, J) %>% distinct() %>%
  ggplot(., aes(y=J, x=siteID)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  labs(title = "Values of evenness for sites ordered by community size (descending)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

  
comp_plot_a / comp_plot_b / comp_plot_c

comp_plot_c / cc_plot


z <- abundances
z %<>% select(siteID,ABLow,AB26,AB42,AB21,AB9,J,total)
z %<>% rename(Low=ABLow)
z %<>% pivot_longer(cols=2:6,names_to="Species",values_to="abund")
#z %<>% filter(J<0.6)
z$siteID <- reorder(z$siteID, -z$total) #order by total community size
p1top <- z %>% ggplot(.,aes(x=siteID,y=abund))+
  geom_bar(stat="identity")+
  ylab("Abundance")+theme_minimal()+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
p1bot<- z %>% ggplot(.,aes(x=siteID,y=abund,fill=Species))+
  geom_bar(position="fill",stat="identity")+
  scale_fill_manual(values = c("#238443", "#78C679", "#C2E699", "#FFFFB2", "black"))+
  theme_minimal()+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+ylab("Relative abundance")
p1top/evenness_plot/p1bot/cc_plot


abundances  %>% dplyr::select(siteID, ABLow, AB9, AB21, AB26, AB42) %>% melt() %>%
  ggplot(., aes(y = value, x = siteID, fill = variable)) +
  geom_bar(position ="stack", stat = "identity") +
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_manual(values = c("black", "#ffffb2","#fecc5c", "#fd8d3c", "#e31a1c")) +
  labs(title = "Abundance of high and low competent species for all sites")

abundances %>% dplyr::select(siteID, cc) %>% distinct() %>%
  ggplot(., aes(y=cc, x=siteID)) +
  geom_bar(position = position_dodge(), stat = "identity") +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Values of cc for all sites ordered by community size (descending)")
  


 #fit two lines to cc_evenness_prev plot
lag_evenness$high <- c(0)
lag_evenness$low <- c(0)
for (i in 1:nrow(lag_evenness)){
  if (lag_evenness$lag_J[i] > 0.6){
    lag_evenness$high[i] <- 1
    lag_evenness$low[i] <- 1
  }
  else if (lag_evenness$lag_cc[i] < 45000){
    lag_evenness$low[i] <- 1
  }
  else if (lag_evenness$lag_cc[i] > 45000){
    lag_evenness$high[i] <- 1
  }
  else {
    lag_evenness$high[i] <- 0
    lag_evenness$low[i] <- 0
  }
}
low_cc <- lag_evenness %>% filter(low==1)
cor.test(x=low_cc$lag_J, y=low_cc$lag_cc,method = "spearman")
high_cc <- lag_evenness %>% filter(high==1)
cor.test(x=high_cc$lag_J, y=high_cc$lag_cc,method = "spearman")

p4 <- ggplot() + 
  geom_point(data = lag_evenness, aes(x=lag_J, y=lag_cc, color = Prevalence, size = lag_size)) +
  scale_color_gradient(low="blue", high="red") + 
  new_scale_color() +
  geom_smooth(data = low_cc, aes(x=lag_J, y=lag_cc), method = "lm", show.legend = F) +
  new_scale_color() +
  geom_smooth(data = high_cc, aes(x=lag_J, y=lag_cc), method = "lm", show.legend = F) +
  theme_minimal()
p4  


mySpan=1.5
myColor="gray80"
p4 <- ggplot() +
  geom_point(data = lag_evenness, aes(x=lag_J, y=lag_cc, color = Prevalence, size = lag_size)) +
  scale_color_gradient(low="blue", high="red") +
  new_scale_color() +
  geom_smooth(data = low_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  new_scale_color() +
  geom_smooth(data = high_cc, aes(x=lag_J, y=lag_cc), method = "loess",span=mySpan,se=F,col=myColor, show.legend = F) +
  theme_minimal()+ylim(0,3e+05)+
  annotate("text", label = "i", family="Times", fontface="italic", x = 0.134, y = 46000, size = 4, colour = "white")+
  annotate("text", label = "ii", family="Times", fontface="italic", x = 0.427, y = 61000, size = 3, colour = "white")
p4
