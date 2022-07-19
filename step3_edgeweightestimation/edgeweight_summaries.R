library(tidyverse)
##### MOTNP #####
motnp_sightings <- read_csv('data_processed/motnp_eles_long_22.01.13.csv') %>% 
  separate(elephant, into = c('sex','casename'), 1, remove = F) %>% 
  filter(sex == 'M') %>% 
  select(-casename, -sex, -number) %>% 
  distinct()
colnames(motnp_sightings)[15] <- 'id'
min(motnp_sightings$date) ; max(motnp_sightings$date)
max(motnp_sightings$date) - min(motnp_sightings$date)

motnp_eles <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
male_sightings <- left_join(motnp_sightings, motnp_eles, by = 'id')%>% 
  filter(dem_class == 'AM' | dem_class == 'PM')

length(unique(male_sightings$id))
length(unique(male_sightings$encounter))
counts_per_male <- as.data.frame(table(male_sightings$id_no))
summary(counts_per_male$Freq) ; sd(counts_per_male$Freq)

male_obs <- male_sightings[,c(1:5,7,14)] %>% distinct()
table(male_obs$herd_type, male_obs$total_elephants_numeric)
summary(male_obs$total_elephants_numeric) ; sd(male_obs$total_elephants_numeric)
quantile(male_obs$total_elephants_numeric, seq(0,1,0.01))

male_obs %>% filter(herd_type == 'MO') %>% summary(total_elephants_numeric)

rm(motnp_eles, motnp_sightings, male_obs, male_sightings, counts_per_male)

# loads counts_df
counts_df <- read_delim('data_processed/motnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',')
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
counts_df$age_category_1 <- ifelse(counts_df$age_category_1 == '1-2','0-3',counts_df$age_category_1)
counts_df$age_cat_id_1 <- ifelse(counts_df$age_category_1 == '0-3', 1,
                                 ifelse(counts_df$age_category_1 == '3-4', 1,
                                        ifelse(counts_df$age_category_1 == '4-5', 1,
                                               ifelse(counts_df$age_category_1 == '5-6', 2,
                                                      ifelse(counts_df$age_category_1 == '6-7', 2,
                                                             ifelse(counts_df$age_category_1 == '7-8', 2,
                                                                    ifelse(counts_df$age_category_1 == '8-9', 2,
                                                                           ifelse(counts_df$age_category_1 == '9-10', 2,
                                                                                  ifelse(counts_df$age_category_1 == '10-15', 3,
                                                                                         ifelse(counts_df$age_category_1 == '15-19', 4,
                                                                                                ifelse(counts_df$age_category_1 == '20-25', 5,
                                                                                                       ifelse(counts_df$age_category_1 == '20-35', 5,
                                                                                                              ifelse(counts_df$age_category_1 == '25-40', 6,
                                                                                                                     ifelse(counts_df$age_category_1 == '35-50', 6,
                                                                                                                            ifelse(counts_df$age_category_1 == '40+', 7,
                                                                                                                                   ifelse(counts_df$age_category_1 == '50+', 7, counts_df$age_category_1))))))))))))))))
counts_df$age_cat_id_1[which(is.na(counts_df$age_cat_id_1))] <- 1 # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_class_1 <- ifelse(counts_df$age_cat_id_1 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_1 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_1 > 4, 'Adult','Pubescent')))
unique(counts_df$age_category_2) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
counts_df$age_category_2 <- ifelse(counts_df$age_category_2 == '1-2','0-3',counts_df$age_category_2)
counts_df$age_cat_id_2 <- ifelse(counts_df$age_category_2 == '0-3', 1,
                                 ifelse(counts_df$age_category_2 == '3-4', 1,
                                        ifelse(counts_df$age_category_2 == '4-5', 1,
                                               ifelse(counts_df$age_category_2 == '5-6', 2,
                                                      ifelse(counts_df$age_category_2 == '6-7', 2,
                                                             ifelse(counts_df$age_category_2 == '7-8', 2,
                                                                    ifelse(counts_df$age_category_2 == '8-9', 2,
                                                                           ifelse(counts_df$age_category_2 == '9-10', 2,
                                                                                  ifelse(counts_df$age_category_2 == '10-15', 3,
                                                                                         ifelse(counts_df$age_category_2 == '15-19', 4,
                                                                                                ifelse(counts_df$age_category_2 == '20-25', 5,
                                                                                                       ifelse(counts_df$age_category_2 == '20-35', 5,
                                                                                                              ifelse(counts_df$age_category_2 == '25-40', 6,
                                                                                                                     ifelse(counts_df$age_category_2 == '35-50', 6,
                                                                                                                            ifelse(counts_df$age_category_2 == '40+', 7,
                                                                                                                                   ifelse(counts_df$age_category_2 == '50+', 7, counts_df$age_category_2))))))))))))))))
counts_df$age_cat_id_2[which(is.na(counts_df$age_cat_id_2))] <- 1   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_class_2 <- ifelse(counts_df$age_cat_id_2 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_2 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_2 > 4, 'Adult','Pubescent')))
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))
counts_df$age_class_id_1 <- ifelse(counts_df$age_class_1 == 'Adult',4,
                                   ifelse(counts_df$age_class_1 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_1 == 'Juvenile',2,1)))
counts_df$age_class_id_2 <- ifelse(counts_df$age_class_2 == 'Adult',4,
                                   ifelse(counts_df$age_class_2 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_2 == 'Juvenile',2,1)))
counts_df$dem_type <- ifelse(counts_df$age_class_id_1 > counts_df$age_class_id_2,
                             paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'),
                             ifelse(counts_df$age_class_id_1 < counts_df$age_class_id_2,
                                    paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_'),
                                    paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_')))
sort(unique(counts_df$dem_type))
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.
counts_df$apart <- counts_df$count_dyad - counts_df$all_events
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1

# load model output
draws_motnp2.2 <- readRDS('data_processed/motnp_bayesian_edgedistributions_a2.b2_22.02.07.rds') %>% 
  data.matrix()

# summarise -- look for any anomalies in draw values or chain variation
summaries <- data.frame(dyad = colnames(draws_motnp2.2[,2:106954]),
                        min = rep(NA, ncol(draws_motnp2.2)-1),
                        max = rep(NA, ncol(draws_motnp2.2)-1),
                        mean = rep(NA, ncol(draws_motnp2.2)-1),
                        median = rep(NA, ncol(draws_motnp2.2)-1),
                        sd = rep(NA, ncol(draws_motnp2.2)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_motnp2.2[,i+1])
  summaries$max[i]    <- max(draws_motnp2.2[,i+1])
  summaries$mean[i]   <- mean(draws_motnp2.2[,i+1])
  summaries$median[i] <- median(draws_motnp2.2[,i+1])
  summaries$sd[i]     <- sd(draws_motnp2.2[,i+1])
}

summary(summaries$median)
summary(summaries$mean) ; sd(summaries$mean)

all_data <- left_join(counts_df, summaries, by = 'dyad')
str(all_data)

sort(unique(all_data$dem_type))
male_data <- all_data %>% filter(dem_type == 'AM_AM' | dem_type == 'AM_PM' | dem_type == 'PM_PM')
#non_males_males <- anti_join(all_data, male_data) %>% filter(dem_class_1 == 'AM')
#sort(unique(non_males_males$dem_class_2))
#non_males_males <- anti_join(all_data, male_data) %>% filter(dem_class_2 == 'AM')
#sort(unique(non_males_males$dem_class_1))

summary(male_data$median) ; sd(male_data$median)
summary(male_data$mean) ; sd(male_data$mean)
summary(male_data$sd)

summaries$period <- 1
medians_motnp <- summaries[,c('period','median')]
ggplot(medians_motnp)+
  geom_density(aes(x = median, colour = factor(period)), lwd = 0.4)+
  xlab('median dyad draw value from MCMC chains')+
  theme_light()+
  theme(legend.position = 'none')

rm(counts_df, summaries, all_data, male_data, non_males, non_males_males, male_sightings, draws_motnp2.2, counts_df, motnp_eles, motnp_sightings, summaries)

##### MPNP #####
mpnp <- read_csv('data_processed/mpnp_eles_long_22.03.08.csv') %>% 
  separate(elephant, into = c('BT','num'), 1, remove = F) %>% 
  filter(BT == 'B') %>% 
  select(-BT, -num)
length(unique(mpnp$elephant))
counts_per_male <- as.data.frame(table(mpnp$elephant))
summary(counts_per_male$Freq) ; sd(counts_per_male$Freq)
summary(mpnp$total_elephants) ; sd(mpnp$total_elephants)
length(which(mpnp$total_elephants == 0))
check <- mpnp[mpnp$total_elephants == 0,] # these make no sense at all...

mpnp <- read_csv('data_processed/mpnp_eles_long_22.03.08.csv') %>% 
  separate(elephant, into = c('BT','num'), 1, remove = F) %>% 
  filter(BT == 'B') %>% 
  select(-BT, -num) %>% 
  filter(total_elephants > 0)

periods <- read_csv('data_processed/mpnp_dyad_weightdistributions_2.2_allperiods_22.06.01.csv')
unique(periods$period_start)
periods$period_start[which(periods$period_start == "2013-11-26")[1]] - periods$period_start[1] # 571 days
periods$period_start[which(periods$period_start == "2018-08-07")[1]] + 571 # "2020-02-29"
mpnp <- mpnp %>% filter(date < periods$period_start[which(periods$period_start == "2018-08-07")[1]] + 571)

length(unique(mpnp$elephant))
counts_per_male <- as.data.frame(table(mpnp$elephant))
summary(counts_per_male$Freq) ; sd(counts_per_male$Freq)
summary(mpnp$total_elephants) ; sd(mpnp$total_elephants)

min(mpnp$date) ; max(mpnp$date)

mpnp_sightings <- mpnp[,c(4:8)] %>% distinct()

length(which(mpnp_sightings$total_elephants == 1))
length(which(mpnp_sightings$total_elephants > 1))

summaries <- read_csv('data_processed/mpnp_dyad_weightdistributions_2.2_allperiods_22.06.01.csv')
summary(summaries[which(summaries$period == 1),]$median)
summary(summaries[which(summaries$period == 2),]$median)
summary(summaries[which(summaries$period == 3),]$median)
summary(summaries[which(summaries$period == 4),]$median)
summary(summaries[which(summaries$period == 5),]$median)
sd(summaries[which(summaries$period == 1),]$median)
sd(summaries[which(summaries$period == 2),]$median)
sd(summaries[which(summaries$period == 3),]$median)
sd(summaries[which(summaries$period == 4),]$median)
sd(summaries[which(summaries$period == 5),]$median)

medians_mpnp <- summaries[,c('period','median')]
ggplot(medians_mpnp)+
  geom_density(aes(x = median, colour = factor(period)), lwd = 0.4)+
  xlab('median dyad draw value from MCMC chains')+
  theme_light()+
  theme(legend.position = 'none')

rm(mpnp, mpnp_sightings, counts_per_male, check, summaries)
##### ANP #####
anp <- read_csv('data_processed/anp_sightings_updated_22.06.22.csv')
min(anp$obs_date) ; max(anp$obs_date)
length(unique(anp$casename))
counts_per_male <- as.data.frame(table(anp$casename))
mean(counts_per_male$Freq) ; sd(counts_per_male$Freq)
mean(anp$grp_size) ; sd(anp$grp_size)

anp$obs_type # M = Mixed, B = Bull
anp <- anp[,c(4,5,7,13)] %>% distinct()
table(anp$obs_type)

ate <- readxl::read_excel(path = 'data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1) %>% janitor::clean_names()
ate$id <- paste('M',ate$casename, sep = '')
ate$node_id <- as.integer(as.factor(ate$casename))
ate$obs_date <- lubridate::as_date(ate$obs_date)
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)
lu <- function(x) { length(unique(x)) }
ate_nums <- tapply(X = ate$obs_num, INDEX = ate$obs_date, FUN = lu )
ate$obs_num <- ifelse(ate$obs_num == '0','00', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '0a','0A', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '0b','0B', ate$obs_num)
ate$obs_num <- ifelse(ate$obs_num == '1','01', ate$obs_num)
ate$obs_num_std <- NA
for(i in 1:length(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_std <- as.integer(as.factor(sort(date_row$obs_num)))
  ate$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}
ate <- ate[,c(1:3,25,26,4,27,28,8,29,9:24)]
sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_std','grid_code')] %>% distinct()
sightings <- sightings[c(1:12,15:24176),]
periods <- data.frame(period = 1:31,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = 32)[1:31],
                      period_end = seq(from = min(sightings$obs_date),
                                       to = max(sightings$obs_date),
                                       length.out = 32)[2:32])
periods$period_end[1] - periods$period_start[1] # 580 days

summaries <- read_csv('data_processed/anp_dyad_weightdistributions_2.2_period1to28_22.05.06.csv')
summary(summaries)

for(i in 1:28){
  print(paste0(i, ' min: ', round(min(summaries[which(summaries$period == i),]$median),3)))
  print(paste0(i, ' median: ', round(median(summaries[which(summaries$period == i),]$median),3)))
  print(paste0(i, ' mean: ', round(mean(summaries[which(summaries$period == i),]$median),3)))
  print(paste0(i, ' max: ', round(max(summaries[which(summaries$period == i),]$median),3)))
  print(paste0(i, ' sd: ', round(sd(summaries[which(summaries$period == i),]$median),3)))
}

distributions <- data.frame(period = 1:28,
                            min = rep(NA, 28),
                            median = rep(NA, 28),
                            mean = rep(NA, 28),
                            max = rep(NA, 28))
for(i in 1:nrow(distributions)){
  data <- summaries[summaries$period == i,]
  distributions$min[i]  <- min(data$median)
  distributions$median[i]  <- median(data$median)
  distributions$mean[i] <- mean(data$median)
  distributions$max[i]  <- max(data$median)
}

plot_dists <- pivot_longer(distributions, cols = c('min','median','mean','max'),
                           names_to = 'measure', values_to = 'draw_value')

ggplot(plot_dists, aes(x = draw_value))+
  geom_density()+
  facet_wrap(~ measure, nrow = 2, scales = 'free_y')+
  theme(strip.background = element_rect(fill = rgb(0,0,0.8,0.5)))

plot_dists_all <- pivot_longer(summaries[,c(7,31:35)], cols = c('min','median','mean','max'),
                               names_to = 'measure', values_to = 'draw_value')
ggplot(plot_dists_all, aes(x = draw_value, colour = factor(period)))+
  geom_density()+
  facet_wrap(~ measure, nrow = 2, scales = 'free_y')+
  theme_light()+
  theme(legend.position = 'none', strip.background = element_rect(fill = rgb(0,0,1,0.5)))
ggplot(plot_dists_all)+
  geom_density(aes(x = draw_value, colour = factor(period)), lwd = 0.4)+
  facet_wrap(~ measure, nrow = 2, scales = 'free_y')+
  xlab('median dyad draw value from MCMC chains')+
  theme_light()+
  theme(legend.position = 'none', strip.background = element_rect(fill = rgb(0,0,1,0.5)))+
  geom_density(aes(x = draw_value), lwd = 1.2, lty = 3)

medians_anp <- summaries[,c('period','median')]
ggplot(medians_anp)+
  geom_density(aes(x = median, colour = factor(period)), lwd = 0.4)+
  xlab('median dyad draw value from MCMC chains')+
  theme_light()+
  theme(legend.position = 'none')

rm(anp, counts_per_male, ate, sightings, periods)

#### plots ####
medians_motnp$park <- 'MOTNP'
medians_mpnp$park <- 'MPNP'
medians_anp$park <- 'ANP'

par(mai = c(2,2,2,2),
    mar = c(10,10,10,5))
medians <- rbind(medians_motnp, medians_mpnp, medians_anp)
medians$park <- factor(medians$park, levels = c('MOTNP','MPNP','ANP'))
ggplot(medians)+
  geom_density(aes(x = median, colour = factor(period)))+
  geom_density(aes(x = median))+
  facet_wrap(. ~ park)+
  scale_x_continuous(limits = c(0,0.5))+
  theme_light()+
  xlab('median dyad draw value from MCMC chains')+
  theme(legend.position = 'none',
        strip.background = element_rect(fill = rgb(0,0,1,0.5)),
        panel.spacing = unit(0.5, 'cm'),
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10))

mid_averages <- data.frame(period = 1:28,
                           period_median = rep(NA, 28),
                           period_mean = rep(NA, 28))
for(i in 1:nrow(mid_averages)){
  p <- medians_anp[medians_anp$period == mid_averages$period[i],]
  mid_averages$period_median[i] <- median(p$median)
  mid_averages$period_mean[i] <- mean(p$median)
}
anp <- left_join(medians_anp, mid_averages, by = 'period')

ggplot(anp)+
  geom_point(aes(x = period, y = median), colour = rgb(0,0,1,0.3))+
  geom_line(aes(x = period, y = period_median),  colour = 'red', lwd = 2)+
  geom_point(aes(x = period, y = period_median), colour = 'red')+
  theme_light()+
  xlab('time window')+
  ylab('median of median dyad association strength distributions')

mid_averages <- data.frame(period = 1:5,
                           period_median = rep(NA, 5),
                           period_mean = rep(NA, 5))
for(i in 1:nrow(mid_averages)){
  p <- medians_mpnp[medians_anp$period == mid_medians$period[i],]
  mid_averages$period_median[i] <- median(p$median)
  mid_averages$period_mean[i] <- median(p$mean)
}
mpnp <- left_join(medians_mpnp, mid_medians, by = 'period')

ggplot(mpnp)+
  geom_point(aes(x = period, y = median), colour = rgb(0,0,1,0.3))+
  geom_line(aes(x = period, y = period_median), colour = 'red', lwd = 2)+
  theme_light()+
  xlab('time window')+
  ylab('median of median dyad association strength distributions')











