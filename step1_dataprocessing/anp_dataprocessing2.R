# ATE data model prep
### Information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022
### Set up ####
library(tidyverse)
library(rstan)
library(rethinking)
library(cmdstanr)
library(lubridate)

### read in sightings data ####
# new data frame, sent August
ate <- readxl::read_xlsx('../data_raw/Raw_ATE_MaleSightingsCleaned_Fishlock220808.xlsx') %>% 
  janitor::clean_names()

# old data frame, sent February -- add location data
old <- readxl::read_xlsx('../data_raw/old_anp/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% janitor::clean_names()
locations <- old[,c('obs_id','casename','obs_num','utm_lat','utm_long','bull_q_r','obs_type')]
ate <- left_join(ate, locations, by = c('obs_id','casename'))
ate <- ate[,c('obs_id','casename','musth','obs_date','obs_time','obs_num','grid_code','utm_lat','utm_long','hab_code_1','hab_code_2','bull_q_r','obs_type','grp_q_c','grp_size','bull_q_c','num_bulls','bulls_1_2','bulls_3_5','musth_male','oestrus_fem','act_code')]
colnames(ate)[c(6,8,9,12,13)] <- c("obs_num_old","utm_lat_old","utm_long_old","bull_q_r_old","obs_type_old")
rm(locations,old)

# add information
ate$id <- paste('M',ate$casename, sep = '')
ate$node_id <- as.integer(as.factor(ate$casename))
ate$bull_q_c <- as.character(ate$bull_q_c)
ate$grp_q_c <- as.character(ate$grp_q_c)
ate$obs_date <- lubridate::as_date(ate$obs_date)
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)

# match up group numbers between old and new data frames
lu <- function(x) { length(unique(x)) }
ate_nums <- tapply(X = ate$obs_num_old, INDEX = ate$obs_date, FUN = lu )
ate$obs_num_old <- ifelse(ate$obs_num_old == '0','00', ate$obs_num_old)
ate$obs_num_old <- ifelse(ate$obs_num_old == '0a','0A', ate$obs_num_old)
ate$obs_num_old <- ifelse(ate$obs_num_old == '0b','0B', ate$obs_num_old)
ate$obs_num_old <- ifelse(ate$obs_num_old == '1','01', ate$obs_num_old)
ate$obs_num_old_std <- NA
for(i in 1:nrow(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_old_std <- as.integer(as.factor(date_row$obs_num_old))
  ate$obs_num_old_std[i] <- date_row$obs_num_old_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}
ate <- ate[,c(1:3,25:26,4,27:28,8,29,9:24)]
head(ate)

eles <- ate[,c('id','node_id')] %>% distinct()

## clean environment
rm(ate_nums, date_row, i)

## progress report
print(paste0('sightings data imported at ', Sys.time()))

#### Import nodes data ####
nodes <- readxl::read_excel('../data_raw/Raw_ATE_LifeHistoryData_Fishlock220808.xlsx') %>%
  janitor::clean_names() %>%
  distinct()
colnames(nodes)

## make character string of ID number
nodes$id <- paste('M', nodes$casename, sep = '')

## check ID numbers match in both data frames
#nodes_id <- sort(unique(nodes$id))
#ate_id <- sort(unique(ate$id))
#length(nodes_id) ; length(ate_id) ## DO NOT MATCH, OR EVEN CLOSE!

## progress report
print(paste0('nodes data imported at ', Sys.time()))

## read in processed data files
all <- readRDS(file = '../data_processed/anp_bayesian_allpairwiseevents.RDS')

## merge sightings information into dyad data
sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_old_std','grid_code')] %>% distinct()
sightings$obs_id_intfact <- as.integer(as.factor(sightings$obs_id))
nrow(sightings)              # 24386
length(unique(ate$obs_id))   # 24386

colnames(all)[4] <- 'obs_id_intfact'
all <- left_join(x = all, y = sightings, by = 'obs_id_intfact') %>% distinct()
rm(ate)

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
all$dyad <- paste(all$node_1, all$node_2, sep = '_')
all$dyad_id <- as.integer(as.factor(all$dyad))
all$location_id <- as.integer(as.factor(all$grid_code))
head(all)

max(sightings$obs_date) - min(sightings$obs_date) # 17999 days (call that 18000 including both ends!). Time period for ALERT data = 581 days --> split into 31 time windows (so need 32 boundaries)
periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)
events <- all[all$social_event == 1,] # cuts down to 114209 dyad pairs from 38098504 -- assigning to time window doesn't take forever
events$period <- NA
for(i in 1:nrow(events)){
  events$period[i] <- which(periods <= events$obs_date[i])[length(which(periods <= events$obs_date[i]))] # take last value in vector
  if(i %% 1000 == 0) {print(i)}
}
periods ; View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

# check elephants all match up
length(unique(all$node_1))
length(unique(all$node_2))
length(unique(events$node_1))
length(unique(events$node_2))

### convert to Binomial model data format -- aggregate all sightings of each dyad together into a count ####
## all time:
df_agg <- all %>%
  group_by(node_1, node_2) %>%
  summarise(event_count=sum(social_event), dyad_id=cur_group_id())
length(unique(df_agg$node_1)) ; length(unique(df_agg$node_2))
length(df_agg$node_1) == cumsum(1:length(unique(all$node_1)))[length(unique(all$node_1))] # number should be the nth value of the triangular number sequence in which n = total number of elephants in analysis (693). If TRUE, correct number of pairs.
rm(all)

## by time window:
df_split <- events %>%
  group_by(node_1, node_2, period) %>%
  summarise(event_count=sum(social_event),
            dyad_id=cur_group_id())
head(df_split)

colnames(eles) <- c('id_1','node_1')
df <- left_join(df_split, eles, by = 'node_1') %>% distinct()
colnames(eles) <- c('id_2','node_2')
df <- left_join(df, eles, by = 'node_2') %>% distinct()
head(df) ; tail(df)

df$dyad <- paste(df$id_1, df$id_2, sep = '_')
df$dyad_id_period <- df$dyad_id              # every dyad has it's own ID number, including if same dyad in a different time window
df$dyad_id <- as.integer(as.factor(df$dyad)) # every dyad has it's own ID number, but same dyad in different windows share ID number

### create dyad row for all pairs per period ####
dyads <- data.frame(id_1 = rep(sort(eles$id_2), each = nrow(eles)),
                    id_2 = rep(sort(eles$id_2), nrow(eles)))

colnames(eles) <- c('id_1','node_1')
dyads <- left_join(dyads, eles, by = 'id_1')
colnames(eles) <- c('id_2','node_2')
dyads <- left_join(dyads, eles, by = 'id_2')
dyads <- dyads[dyads$node_1 < dyads$node_2,]

dyads <- data.frame(id_1 = rep(dyads$id_1, length(unique(df$period))),
                    id_2 = rep(dyads$id_2, length(unique(df$period))),
                    node_1 = rep(dyads$node_1, length(unique(df$period))),
                    node_2 = rep(dyads$node_2, length(unique(df$period))),
                    period = rep(sort(unique(df$period)), each = nrow(dyads)))
head(df) ; head(dyads)

#rm(date_row, events, df_split) # may need to clear environment a little before running next step
data <- left_join(x = dyads, y = df, by = c('id_1','id_2','period'))
data <- data[,c(1:5,8:10)]
colnames(data)[3:4] <- c('node_1','node_2')
data$event_count <- ifelse(is.na(data$event_count) == TRUE, 0, data$event_count)
table(data$event_count)
data$dyad <- paste(data$id_1, data$id_2, sep = '_')
data$dyad_id <- as.integer(as.factor(data$dyad))
head(data, 20)

periods <- data.frame(period_start = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)[1:31],
                      period = 1:31)
data <- left_join(x = data, y = periods, by = 'period')
head(data)

### add data about nodes ####
males <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste('M',males$casename, sep = '')
males <- males[,c(21,2,5,6,8,9,14,18)]
colnames(males) ; colnames(data)

colnames(males)[1] <- 'id_1'
data <- left_join(x = data, y = males, by = 'id_1')
colnames(data)[c(10:16)] <- c('id_no_1','bmo_1','byr_1','dmo_1','dyr_1','indyr_1','musthyr_1')
colnames(males)[1] <- 'id_2'
data <- left_join(x = data, y = males, by = 'id_2')
colnames(data)[c(17:23)] <- c('id_no_2','bmo_2','byr_2','dmo_2','dyr_2','indyr_2','musthyr_2')
data <- data[,c(8,7,1:5,9,6,10,17,11,18,12,19,13,20,14,21,15,22,16,23)]
head(data)

### write to file ####
readr::write_delim(data, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_pairwiseevents.csv', delim = ',')

## clean environment
rm(list = ls())

### adding in count data ####
# long version (put into data processing script): -- do not rerun unnecessarily, have commented out data read-in lines ####
### add column for total number of sightings per pair
ate <- readxl::read_excel(path = '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1) %>% janitor::clean_names()
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
head(ate)
rm(date_row, ate_nums, lu, i)

sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_std','grid_code')] %>% distinct()
sightings <- sightings[c(1:12,15:24176),]

### import data for aggregated model (binomial)
#counts_df <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_pairwiseevents_22.03.28.csv', delim = ',') %>% 
select(-c('indyr_1','indyr_2','musthyr_1','musthyr_2'))
str(counts_df)

# check all dyads present
length(unique(counts_df$id_1)) ; length(unique(counts_df$id_2))
ids <- unique(c(counts_df$id_1, counts_df$id_2))
nrow(counts_df) == (cumsum(1:length(ids))[length(ids)-1])*length(unique(counts_df$period)) # number of dyads per time period should be the cumulative sum of n-1 where n is the number of elephants in the dataset. counts_df contains this for 31 time periods.

# remove impossible dyads -- dyads where one individual or the other is dead/unborn before end/at start of period
head(counts_df)
counts_df$birth_1 <- ifelse(counts_df$bmo_1 > 0 & counts_df$bmo_1 < 10,
                            paste0(counts_df$byr_1,'-0',counts_df$bmo_1,'-01'),
                            paste0(counts_df$byr_1,'-',counts_df$bmo_1,'-01'))
counts_df$birth_1 <- ifelse(is.na(counts_df$bmo_1) == TRUE,
                            paste0(counts_df$byr_1,'-01-01'),
                            counts_df$birth_1)
counts_df$birth_1 <- lubridate::as_date(counts_df$birth_1)
summary(counts_df$birth_1)

counts_df$birth_2 <- ifelse(counts_df$bmo_2 > 0 & counts_df$bmo_2 < 10,
                            paste0(counts_df$byr_2,'-0',counts_df$bmo_2,'-01'),
                            paste0(counts_df$byr_2,'-',counts_df$bmo_2,'-01'))
counts_df$birth_2 <- ifelse(is.na(counts_df$bmo_2) == TRUE,
                            paste0(counts_df$byr_2,'-01-01'),
                            counts_df$birth_2)
counts_df$birth_2 <- lubridate::as_date(counts_df$birth_2)
summary(counts_df$birth_2)

counts_df$death_1 <- ifelse(counts_df$dmo_1 > 0 & counts_df$dmo_1 < 10 & counts_df$dyr_1 > 100,
                            paste0(counts_df$dyr_1,'-0',counts_df$dmo_1,'-01'),
                            ifelse(counts_df$dmo_1 > 0 & counts_df$dmo_1 >= 10 & counts_df$dyr_1 > 100,
                                   paste0(counts_df$dyr_1,'-',counts_df$dmo_1,'-01'),
                                   'living'))
counts_df$death_1 <- ifelse(is.na(counts_df$dmo_1) == TRUE & counts_df$dyr_1 > 100,
                            paste0(counts_df$dyr_1,'-01-01'),
                            counts_df$death_1)
counts_df$death_1[which(counts_df$death_1 == 'living')] <- '2022-03-28' # set date to when data frame was created, so definitely after last sighting. 
counts_df$death_1 <- lubridate::as_date(counts_df$death_1)
summary(counts_df$death_1)

counts_df$death_2 <- ifelse(counts_df$dmo_2 > 0 & counts_df$dmo_2 < 10 & counts_df$dyr_2 > 100,
                            paste0(counts_df$dyr_2,'-0',counts_df$dmo_2,'-01'),
                            ifelse(counts_df$dmo_2 > 0 & counts_df$dmo_2 >= 10 & counts_df$dyr_2 > 100,
                                   paste0(counts_df$dyr_2,'-',counts_df$dmo_2,'-01'),
                                   'living'))
counts_df$death_2 <- ifelse(is.na(counts_df$dmo_2) == TRUE & counts_df$dyr_2 > 100,
                            paste0(counts_df$dyr_2,'-01-01'),
                            counts_df$death_2)
counts_df$death_2[which(counts_df$death_2 == 'living')] <- '2022-03-28' # set date to when data frame was created, so definitely after last sighting. 
counts_df$death_2 <- lubridate::as_date(counts_df$death_2)
summary(counts_df$death_2)

periods <- data.frame(period = 1:31,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = 32)[1:31],
                      period_end = seq(from = min(sightings$obs_date),
                                       to = max(sightings$obs_date),
                                       length.out = 32)[2:32])

counts_df <- left_join(counts_df, periods, by = 'period')
counts_df <- counts_df[,c(1:7,24,25,9:23)]
colnames(counts_df)[8] <- 'period_start'

counts_df$alive_1 <- ifelse(counts_df$birth_1 > counts_df$period_start, 'No',
                            ifelse(counts_df$death_1 < counts_df$period_end, 'No', 'Yes'))
length(which(is.na(counts_df$alive_1) == TRUE))
table(counts_df$alive_1)
counts_df$alive_2 <- ifelse(counts_df$birth_2 > counts_df$period_start, 'No',
                            ifelse(counts_df$death_2 < counts_df$period_end, 'No', 'Yes'))
length(which(is.na(counts_df$alive_2) == TRUE))
table(counts_df$alive_2)

counts_df <- counts_df[counts_df$alive_1 == 'Yes' & counts_df$alive_2 == 'Yes',1:23] ; counts_df <- counts_df[!is.na(counts_df$dyad),]

# add count data for both elephants to calculate times observed apart
#males <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA ; males$p7 <- NA
males$p8 <- NA ; males$p9 <- NA ; males$p10 <- NA ; males$p11 <- NA ; males$p12 <- NA ; males$p13 <- NA
males$p14 <- NA ; males$p15 <- NA ; males$p16 <- NA ; males$p17 <- NA ; males$p18 <- NA ; males$p19 <- NA
males$p20 <- NA ; males$p21 <- NA ; males$p22 <- NA ; males$p23 <- NA ; males$p24 <- NA ; males$p25 <- NA
males$p26 <- NA ; males$p27 <- NA ; males$p28 <- NA ; males$p29 <- NA ; males$p30 <- NA ; males$p31 <- NA
periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)
for(i in 1:nrow(males)){
  counts <- ate[ate$id == males$id[i],]
  males$count_all[i] <- nrow(counts)
  males$p1[i] <- length(which(counts$obs_date >= periods[1] & counts$obs_date < periods[2]))
  males$p2[i] <- length(which(counts$obs_date >= periods[2] & counts$obs_date < periods[3]))
  males$p3[i] <- length(which(counts$obs_date >= periods[3] & counts$obs_date < periods[4]))
  males$p4[i] <- length(which(counts$obs_date >= periods[4] & counts$obs_date < periods[5]))
  males$p5[i] <- length(which(counts$obs_date >= periods[5] & counts$obs_date < periods[6]))
  males$p6[i] <- length(which(counts$obs_date >= periods[6] & counts$obs_date < periods[7]))
  males$p7[i] <- length(which(counts$obs_date >= periods[7] & counts$obs_date < periods[8]))
  males$p8[i] <- length(which(counts$obs_date >= periods[8] & counts$obs_date < periods[9]))
  males$p9[i] <- length(which(counts$obs_date >= periods[9] & counts$obs_date < periods[10]))
  males$p10[i] <- length(which(counts$obs_date >= periods[10] & counts$obs_date < periods[11]))
  males$p11[i] <- length(which(counts$obs_date >= periods[11] & counts$obs_date < periods[12]))
  males$p12[i] <- length(which(counts$obs_date >= periods[12] & counts$obs_date < periods[13]))
  males$p13[i] <- length(which(counts$obs_date >= periods[13] & counts$obs_date < periods[14]))
  males$p14[i] <- length(which(counts$obs_date >= periods[14] & counts$obs_date < periods[15]))
  males$p15[i] <- length(which(counts$obs_date >= periods[15] & counts$obs_date < periods[16]))
  males$p16[i] <- length(which(counts$obs_date >= periods[16] & counts$obs_date < periods[17]))
  males$p17[i] <- length(which(counts$obs_date >= periods[17] & counts$obs_date < periods[18]))
  males$p18[i] <- length(which(counts$obs_date >= periods[18] & counts$obs_date < periods[19]))
  males$p19[i] <- length(which(counts$obs_date >= periods[19] & counts$obs_date < periods[20]))
  males$p20[i] <- length(which(counts$obs_date >= periods[20] & counts$obs_date < periods[21]))
  males$p21[i] <- length(which(counts$obs_date >= periods[21] & counts$obs_date < periods[22]))
  males$p22[i] <- length(which(counts$obs_date >= periods[22] & counts$obs_date < periods[23]))
  males$p23[i] <- length(which(counts$obs_date >= periods[23] & counts$obs_date < periods[24]))
  males$p24[i] <- length(which(counts$obs_date >= periods[24] & counts$obs_date < periods[25]))
  males$p25[i] <- length(which(counts$obs_date >= periods[25] & counts$obs_date < periods[26]))
  males$p26[i] <- length(which(counts$obs_date >= periods[26] & counts$obs_date < periods[27]))
  males$p27[i] <- length(which(counts$obs_date >= periods[27] & counts$obs_date < periods[28]))
  males$p28[i] <- length(which(counts$obs_date >= periods[28] & counts$obs_date < periods[29]))
  males$p29[i] <- length(which(counts$obs_date >= periods[29] & counts$obs_date < periods[30]))
  males$p30[i] <- length(which(counts$obs_date >= periods[30] & counts$obs_date < periods[31]))
  males$p31[i] <- length(which(counts$obs_date >= periods[31] & counts$obs_date < periods[32]))
}

# counts per dyad per time period
periods <- data.frame(period = 1:31,
                      period_start = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)[1:31],
                      period_end   = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)[2:32])
counts_df$period_count_1 <- NA ; counts_df$period_count_2 <- NA
for(i in 1:nrow(counts_df)){
  counts_1 <- ate[ate$id == counts_df$id_1[i] &
                    ate$obs_date >= periods$period_start[which(counts_df$period[i] == periods$period)] &
                    ate$obs_date < periods$period_end[which(counts_df$period[i] == periods$period)],]
  counts_2 <- ate[ate$id == counts_df$id_2[i] &
                    ate$obs_date >= periods$period_start[which(counts_df$period[i] == periods$period)] &
                    ate$obs_date < periods$period_end[which(counts_df$period[i] == periods$period)],]
  counts_df$period_count_1[i] <- nrow(counts_1)
  counts_df$period_count_2[i] <- nrow(counts_2)
}
str(counts_df)
summary(counts_df$period_count_1)
summary(counts_df$period_count_2)

counts_df$period_count_dyad <- (counts_df$period_count_1 + counts_df$period_count_2) - counts_df$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

# delete dyads where one or other is not seen at any point during a time window -- if never seen, you don't know that they were actually within the population at all for that time, and therefore they should be treated as if you don't know they exist
table(counts_df$period_count_1)
#    0     1     2     3     4     5     6     7    14 
#17169  5740  2293  1935   413   237   233   210   211 
table(counts_df$period_count_2)
#    0     1     2     3     4     5     6     7    14 
#26623  1162   325   207    63     1     5    28    27 
counts_df_non0 <- counts_df[counts_df$period_count_1 > 0,]           # 1581181 to 1051755
table(counts_df_non0$period_count_2) # 649056 where one never seen
counts_df_non0 <- counts_df_non0[counts_df_non0$period_count_2 > 0,] # 1051755 down to 402699

# add column for total number of sightings per pair where they were NOT together
counts_df_non0$apart <- counts_df_non0$period_count_dyad - counts_df_non0$event_count

# create variable for age difference (no point doing the same for sex difference as all are males)
counts_df_non0$age_start_1 <- round(as.numeric(as.character((counts_df_non0$period_start - counts_df_non0$birth_1)/365.25)),0) # approx. age id_1 at start of period
counts_df_non0$age_start_2 <- round(as.numeric(as.character((counts_df_non0$period_start - counts_df_non0$birth_2)/365.25)),0) # approx. age id_2 at start of period
counts_df_non0$age_diff <- abs(counts_df_non0$age_start_1 - counts_df_non0$age_start_2)

# write new data frame to file including counts
write_csv(counts_df_non0, path = '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_pairwiseevents.csv')

## clean environment
rm(list = ls())
