# script to look at number of sightings alone per individual vs age
# load packages
library(tidyverse)

# load data
nodes <- read_csv('data_processed/motnp_elenodes_22.01.13.csv') %>% 
  filter(dem_class == 'AM' | dem_class == 'PM')

sightings <- read_csv('data_processed/motnp_eles_long_22.01.13.csv') %>% 
  select (-number) %>% 
  filter(elephant %in% nodes$id)

solo <- sightings[sightings$total_elephants == 1,]
solo <- solo[solo$total_id_hkm == 1,]
solo <- separate(solo, elephant, into = c('sex','num'), sep = 1, remove = F) %>% 
  filter(sex == 'M') %>% 
  select(-sex, -num)
length(unique(solo$elephant))

true_ages <- readRDS('data_processed/motnp_ageestimates_mcmcoutput_22.07.13.rds')
true_ages <- true_ages[, which(colnames(true_ages) %in% nodes$id)]

age_long <- pivot_longer(data = true_ages, cols = everything(),
                         names_to = 'elephant', values_to = 'age_draw')

true_ages <- as.matrix(true_ages)
age_mean <- data.frame(elephant = colnames(true_ages),
                       mean_age = NA,
                       stdv_age = NA)
for(i in 1:nrow(age_mean)){
  age_mean$mean_age[i] <- mean(true_ages[,i])
  age_mean$stdv_age[i] <- sd(true_ages[,i])
}

#ggplot(age_long, aes(x = elephant, y = age_draw))+
#  geom_point(shape = 1, size = 2, colour = rgb(0,0,1,0.2))

age_long_solo <- age_long[age_long$elephant %in% solo$elephant,]
ggplot(age_long_solo, aes(x = elephant, y = age_draw))+
  geom_point()

# ages vs group size
s <- sightings[,c('encounter','elephant','total_elephants_numeric')]
age_gs_full <- full_join(age_long, s, by = 'elephant')
age_gs_mean <- full_join(age_mean, s, by = 'elephant')
age_gs_mean$lwr_age <- age_gs_mean$mean_age - age_gs_mean$stdv_age
age_gs_mean$upr_age <- age_gs_mean$mean_age + age_gs_mean$stdv_age
age_gs_mean$max_gs <- NA
age_gs_mean$mean_gs <- NA
age_gs_mean$stdv_gs <- NA
for(i in 1:nrow(age_gs_mean)){
  ele <- age_gs_mean[age_gs_mean$elephant == age_gs_mean$elephant[i],]
  age_gs_mean$max_gs[i] <- ifelse(age_gs_mean$total_elephants_numeric[i] == max(ele$total_elephants_numeric), 
                                  'max','not_max')
  age_gs_mean$mean_gs[i] <- mean(ele$total_elephants_numeric, na.rm = T)
  age_gs_mean$stdv_gs[i] <- sd(ele$total_elephants_numeric, na.rm = T)
  rm(ele)
}

plot(age_gs_mean$total_elephants_numeric ~ age_gs_mean$mean_age,
     pch = 19, col = rgb(0,0,1,0.3), las = 1,
     xlab = 'mean age estimate',
     ylab = 'total elephants in group')
for(i in 1:nrow(age_gs_mean)){
  lines(x = c(age_gs_mean$lwr_age[i], age_gs_mean$upr_age[i]),
        y = c(age_gs_mean$total_elephants_numeric[i], age_gs_mean$total_elephants_numeric[i]),
        col = rgb(0,0,1,0.3))
}

ggplot(age_gs_mean[age_gs_mean$max_gs == 'max',], aes(x = mean_age, y = total_elephants_numeric))+
  geom_point(colour = rgb(0,0,1,0.3))

plot(age_gs_mean$mean_gs ~ jitter(x = age_gs_mean$mean_age, amount = 2),
     pch = 19, col = rgb(0,0,1,0.3), las = 1,
     xlab = 'mean age estimate',
     ylab = 'total elephants in group')
for(i in 1:nrow(age_gs_mean)){
  lines(x = c(age_gs_mean$lwr_age[i], age_gs_mean$upr_age[i]),
        y = c(age_gs_mean$total_elephants_numeric[i], age_gs_mean$total_elephants_numeric[i]),
        col = rgb(0,0,1,0.3))
}

plot_means <- age_gs_mean[,c('elephant','mean_age','stdv_age','mean_gs','stdv_gs')] %>% distinct()
ggplot(plot_means)+
  geom_errorbar(aes(x = jitter(mean_age, amount = 2), ymin = mean_gs - stdv_gs, ymax = mean_gs + stdv_gs),
                colour = rgb(0.5,0,1,0.5))+
  labs(x = 'mean age estimate (jittered by 2)', y = 'mean group size')+
  theme_classic()

plot_means$age_cat <- ifelse(plot_means$mean_age < 10, 1,
                             ifelse(plot_means$mean_age < 15, 2,
                                    ifelse(plot_means$mean_age < 20, 3,
                                           ifelse(plot_means$mean_age < 25, 4,
                                                  ifelse(plot_means$mean_age < 40, 5, 6)))))
plot_means$index <- NA
for(i in 1:nrow(plot_means)){
  ele <- plot_means[plot_means$age_cat == plot_means$age_cat[i],]
  plot_means$index[i] <- ifelse(is.na(ele$index[1]) == TRUE, 1, max(ele$index, na.rm = T) + 1)
  rm(ele)
}
ggplot(plot_means)+
  geom_point(aes(x = index, y = mean_gs), colour = rgb(0.5,0,1,0.5))+
  geom_errorbar(aes(x = index, ymin = mean_gs - stdv_gs, ymax = mean_gs + stdv_gs),
                colour = rgb(0.5,0,1,0.5))+
  labs(x = 'individual', y = 'mean group size')+
  theme_classic()+
  theme(axis.ticks.x = element_blank(),axis.text.x = element_blank())+
  facet_wrap(. ~ age_cat, scales = 'free_x')











