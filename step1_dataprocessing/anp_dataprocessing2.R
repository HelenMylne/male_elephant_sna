# ATE data model prep
### information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022
### set up ####
library(tidyverse)
library(rstan)
library(rethinking)
library(cmdstanr)
library(lubridate)

#load('ate_dataprocessing1_preparallelisation.RData')

############### produce data for 500 day time windows (match MOTNP) ####
### read in data ####
# raw sightings data
#old <- readxl::read_excel(path = '../data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1) %>%
#  janitor::clean_names()
#old$id <- paste('M',old$casename, sep = '')
#old$node_id <- as.integer(as.factor(old$casename))
#old$obs_date <- lubridate::as_date(old$obs_date)
#old <- separate(old, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
#old$correct_time_hms <- hms::as_hms(old$correct_time)
#old$corrected_time <- lubridate::hour(old$correct_time_hms)*60*60 + lubridate::minute(old$correct_time_hms) + lubridate::second(old$correct_time_hms)
#lu <- function(x) { length(unique(x)) }
#old_nums <- tapply(X = old$obs_num, INDEX = old$obs_date, FUN = lu )
#old$obs_num <- ifelse(old$obs_num == '0','00', old$obs_num)
#old$obs_num <- ifelse(old$obs_num == '0a','0A', old$obs_num)
#old$obs_num <- ifelse(old$obs_num == '0b','0B', old$obs_num)
#old$obs_num <- ifelse(old$obs_num == '1','01', old$obs_num)
#old$obs_num_std <- NA
#for(i in 1:length(old)){
#  date_row <- old[old$obs_date == old$obs_date[i],]
#  date_row$obs_num_std <- as.integer(as.factor(sort(date_row$obs_num)))
#  old$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == old$obs_id[i])[1]]
#}
#old <- old[,c(1:3,25,26,4,27,28,8,29,9:24)]
#head(old)
#rm(date_row, old_nums, lu, i)

ate <- readxl::read_excel(path = '../data_raw/Raw_ATE_MaleSightingsCleaned_Fishlock220808.xlsx', sheet = 1) %>%
  janitor::clean_names()
ate$id <- paste('M',ate$casename, sep = '')
ate$node_id <- as.integer(as.factor(ate$casename))
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'),
                remove = F, sep = ' ')
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)
ate$obs_num_std <- NA
for(i in 1:nrow(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_std <- as.integer(as.factor(sort(date_row$correct_time)))
  ate$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}
ate <- ate[,c(8,1:3,21,22,4,23,24,25,13,9)]
head(ate)
rm(date_row, i) ; gc()

## read in processed data files
anp <- readRDS('../data_processed/anp_allpairwiseevents/anp_bayesian_allpairwiseevents_block1.RDS') %>% as.data.frame()

# check that data processing has worked properly
block_test <- data.frame(block = 1:488, check = NA)
for(i in 2:488){
  # read in block data RDS and append to rest of data
  new_block <- readRDS(paste0('../data_processed/anp_allpairwiseevents/anp_bayesian_allpairwiseevents_block',i,'.RDS')) %>% as.data.frame()
  anp <- rbind(anp, new_block)
  # compare number of rows in new_block to total number that there should be according to gbi_matrix
  gbi_check <- as.data.frame(gbi_matrix[((i*50)-49):(i*50),])
  gbi_check$group_size <- rowSums(gbi_check)
  gbi_check$gs_x_totaleles <- gbi_check$group_size * 689
  x <- as.data.frame(table(new_block$obs_id))
  x$check <- gbi_check$gs_x_totaleles
  block_test$check <- length(which(x$Freq != x$check)) # all have the right number of elephants
}
table(block_test$check) # should all be 0 if it has worked correctly

# save as a single data file
saveRDS(anp, '../data_processed/step1_dataprocessing/anp_allpairwiseevents/anp_bayesian_allpairwiseevents.RDS')
#anp <- readRDS(file = '../data_processed/step1_dataprocessing/anp_allpairwiseevents/anp_bayesian_allpairwiseevents.RDS')

## merge sightings information into dyad data
sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_std', 'grid_code')] %>% distinct()
sightings$obs_id_intfact <- as.integer(as.factor(sightings$obs_id))
length(unique(sightings$obs_id_intfact))  # 24386
nrow(sightings)                           # 24386
length(unique(ate$obs_id))                # 24386

colnames(anp)[4] <- 'obs_id_intfact'
anp <- left_join(x = anp, y = sightings, by = 'obs_id_intfact') %>% distinct()

# clean environment
rm(ate_asnipe, block_table, block_test, factor_blocks, gbi_check, new_block, x, block_size, factor_levels, i, length_na) ; gc()

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
anp$dyad <- paste(anp$node_1, anp$node_2, sep = '_')
anp$dyad_id <- as.integer(as.factor(anp$dyad))
anp$location_id <- as.integer(as.factor(anp$grid_code))
head(anp)

max(sightings$obs_date) - min(sightings$obs_date) # 18197 days. Time period for ALERT data = 504 days --> split into 36 time windows (so need 37 boundaries)
num_periods <- 36
periods <- data.frame(period = 1:num_periods,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = (num_periods+1))[1:num_periods],
                      period_end = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = (num_periods+1))[2:(num_periods+1)])
periods$period_start <- as.Date(periods$period_start)
periods$period_end <- as.Date(periods$period_end)

events <- anp[anp$social_event == 1,] ; events <- events[!is.na(events$node_1),] # cuts down to 114188 dyad pairs from 38098504 -- assigning to time window doesn't take forever
events$period <- NA
for(i in 1:nrow(events)){
  x <- which(periods$period_start <= events$obs_date[i]) # vector of which periods started before observation
  events$period[i] <- x[length(x)] # take last value in vector (most recent = current)
  if(i %% 1000 == 0) {print(i)}
}
periods$period_start ; View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

# check elephants all match up
length(unique(anp$node_1))    ; length(unique(anp$node_2))    # 690 elephants for both variables
length(unique(events$node_1)) ; length(unique(events$node_2)) # don't match -- this is possible without an error if by chance of numbering, some individuals only interact with ones with lower numbers (e.g. elephant 3 is only seen once and that is with node 2 so he is in events as node_2 but not as node_1).

sort(unique(events$node_1)) ; sort(unique(events$node_2)) # node 3 only ever in node 2, node 5 never in there, 18 only in node 1
test3 <- anp[anp$node_1 == 3,] # all are no interaction
test3 <- anp[anp$node_2 == 3,] # only one interaction

test5 <- anp[anp$node_1 == 5,] # all are no interaction
test5 <- anp[anp$node_2 == 5,] # all are no interaction

test18 <- anp[anp$node_1 == 18,] # only with 19 and 49
test18 <- anp[anp$node_2 == 18,] # all are no interaction

rm(test3, test5, test18) ; gc()

### convert to Binomial model data format -- aggregate all sightings of each dyad together into a count ####
df_split <- events %>%
  group_by(node_1, node_2, period) %>%
  summarise(event_count = sum(social_event),
            dyad_id = cur_group_id())
head(df_split)

eles <- ate[,c('id','node_id')] %>% distinct()
colnames(eles) <- c('id_1','node_1')
df <- left_join(df_split, eles, by = 'node_1') %>% distinct()
colnames(eles) <- c('id_2','node_2')
df <- left_join(df, eles, by = 'node_2') %>% distinct()
head(df) ; tail(df)

df$dyad <- paste(df$id_1, df$id_2, sep = '_')
#df$dyad_id_period <- df$dyad_id              # every dyad has it's own ID number, including if same dyad in a different time window
#df$dyad_id <- as.integer(as.factor(df$dyad)) # every dyad has it's own ID number, but same dyad in different windows share ID number

### create dyad row for all pairs per period ####
# create data frame containing all dyad pairs
colnames(eles) <- c('id','node') # unnecessary, but makes it clearer that these are not part of a dyad but rather all of the elephants in the population
dyads <- data.frame(id_1 = rep(sort(eles$id), each = nrow(eles)),
                    id_2 = rep(sort(eles$id), nrow(eles)),
                    node_1 = rep(sort(eles$node), each = nrow(eles)),
                    node_2 = rep(sort(eles$node), nrow(eles))
                    )

# duplicate data frame for every period
dyads <- data.frame(id_1 = rep(dyads$id_1, num_periods),
                    id_2 = rep(dyads$id_2, num_periods),
                    node_1 = rep(dyads$node_1, num_periods),
                    node_2 = rep(dyads$node_2, num_periods),
                    period = rep(sort(unique(df$period)), each = nrow(dyads)))
head(df) ; head(dyads)

#rm(events, df_split) # may need to clear environment a little before running next step
data <- left_join(x = dyads, y = df, by = c('id_1','id_2','period'))
data <- data[,c(1:5,8:10)]
colnames(data)[3:4] <- c('node_1','node_2')
data$event_count <- ifelse(is.na(data$event_count) == TRUE, 0, data$event_count)
table(data$event_count)
data$dyad <- paste(data$id_1, data$id_2, sep = '_')
data$dyad_id <- as.integer(as.factor(data$dyad))
head(data, 20)

data <- left_join(x = data, y = periods, by = 'period')
head(data)

### add data about nodes ####
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste('M',males$casename, sep = '')
males <- males[,c(21,2,5,6,8,9,14,18)]
colnames(males)
colnames(data)

colnames(males)[1] <- 'id_1'
data <- left_join(x = data, y = males, by = 'id_1')
colnames(data)[c(11:17)] <- c('id_no_1','bmo_1','byr_1','dmo_1','dyr_1','indyr_1','musthyr_1')
colnames(males)[1] <- 'id_2'
data <- left_join(x = data, y = males, by = 'id_2')
colnames(data)[c(18:24)] <- c('id_no_2','bmo_2','byr_2','dmo_2','dyr_2','indyr_2','musthyr_2')
data <- data[,c(8,7,1:5,9,10,6,11,18,12,19,13,20,14,21,15,22,16,23,17,24)]
head(data)

### write to file ####
data <- data[data$id_1 != data$id_2,]
data <- data[data$node_1 < data$node_2,]

readr::write_delim(data, '../data_processed/step1_dataprocessing/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows.csv', delim = ',')

rm(data, df, df_split, dyads, events, i, x) ; gc()

### remove impossible dyads -- dyads where one individual or the other is dead/unborn before end/at start of period ####
# import data for aggregated model (binomial)
#counts_df <- read_delim('../data_processed/step1_dataprocessing/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows.csv', delim = ',')
str(counts_df)

# check all dyads present
length(unique(counts_df$id_1)) ; length(unique(counts_df$id_2))
ids <- unique(c(counts_df$id_1, counts_df$id_2))
nrow(counts_df) == (cumsum(1:length(ids))[length(ids)-1])*length(unique(counts_df$period)) # number of dyads per time period should be the cumulative sum of n-1 where n is the number of elephants in the dataset. counts_df contains this for 36 time periods.

# estimate missing birth dates and standardise formatting
head(counts_df)
counts_df$birth_1 <- ifelse(counts_df$bmo_1 > 0 & counts_df$bmo_1 < 10,   # combine birth year and month to single value, set day to 1st of month
                            paste0(counts_df$byr_1,'-0',counts_df$bmo_1,'-01'),
                            paste0(counts_df$byr_1,'-',counts_df$bmo_1,'-01'))
counts_df$birth_1 <- ifelse(is.na(counts_df$bmo_1) == TRUE,   # if birth month unknown, make January of that year
                            paste0(counts_df$byr_1,'-01-01'),
                            counts_df$birth_1)
counts_df$birth_1 <- lubridate::as_date(counts_df$birth_1)    # convert to date format
length(which(is.na(counts_df$birth_1) == TRUE))
summary(counts_df$birth_1)

counts_df$birth_2 <- ifelse(counts_df$bmo_2 > 0 & counts_df$bmo_2 < 10,
                            paste0(counts_df$byr_2,'-0',counts_df$bmo_2,'-01'),
                            paste0(counts_df$byr_2,'-',counts_df$bmo_2,'-01'))
counts_df$birth_2 <- ifelse(is.na(counts_df$bmo_2) == TRUE,
                            paste0(counts_df$byr_2,'-01-01'),
                            counts_df$birth_2)
counts_df$birth_2 <- lubridate::as_date(counts_df$birth_2)
length(which(is.na(counts_df$birth_2) == TRUE))
summary(counts_df$birth_2)

counts_df$death_1 <- ifelse(counts_df$dmo_1 > 0 & counts_df$dmo_1 < 10 & counts_df$dyr_1 > 100,
                            paste0(counts_df$dyr_1,'-0',counts_df$dmo_1,'-01'),
                            ifelse(counts_df$dmo_1 > 0 & counts_df$dmo_1 >= 10 & counts_df$dyr_1 > 100,
                                   paste0(counts_df$dyr_1,'-',counts_df$dmo_1,'-01'),
                                   'living'))
counts_df$death_1 <- ifelse(is.na(counts_df$dmo_1) == TRUE & counts_df$dyr_1 > 100,
                            paste0(counts_df$dyr_1,'-01-01'),
                            counts_df$death_1)
counts_df$death_1[which(counts_df$death_1 == 'living')] <- '2022-08-08' # set date to when raw data were sent over, so definitely after last sighting.
sort(unique(counts_df$death_1))
counts_df$death_1 <- lubridate::as_date(counts_df$death_1)
length(which(is.na(counts_df$death_1) == TRUE))
summary(counts_df$death_1)

counts_df$death_2 <- ifelse(counts_df$dmo_2 > 0 & counts_df$dmo_2 < 10 & counts_df$dyr_2 > 100,
                            paste0(counts_df$dyr_2,'-0',counts_df$dmo_2,'-01'),
                            ifelse(counts_df$dmo_2 > 0 & counts_df$dmo_2 >= 10 & counts_df$dyr_2 > 100,
                                   paste0(counts_df$dyr_2,'-',counts_df$dmo_2,'-01'),
                                   'living'))
counts_df$death_2 <- ifelse(is.na(counts_df$dmo_2) == TRUE & counts_df$dyr_2 > 100,
                            paste0(counts_df$dyr_2,'-01-01'),
                            counts_df$death_2)
counts_df$death_2[which(counts_df$death_2 == 'living')] <- '2022-08-08' # set date to when raw data were sent over, so definitely after last sighting.
counts_df$death_2 <- lubridate::as_date(counts_df$death_2)
length(which(is.na(counts_df$death_2) == TRUE))
summary(counts_df$death_2)

counts_df$alive_1 <- ifelse(counts_df$birth_1 > counts_df$period_start, 'No',
                            ifelse(counts_df$death_1 < counts_df$period_end, 'No', 'Yes'))
length(which(is.na(counts_df$alive_1) == TRUE))
table(counts_df$alive_1)
counts_df$alive_2 <- ifelse(counts_df$birth_2 > counts_df$period_start, 'No',
                            ifelse(counts_df$death_2 < counts_df$period_end, 'No', 'Yes'))
length(which(is.na(counts_df$alive_2) == TRUE))
table(counts_df$alive_2)

counts_df <- counts_df[counts_df$alive_1 == 'Yes' & counts_df$alive_2 == 'Yes',1:24] ; counts_df <- counts_df[!is.na(counts_df$dyad),]

### add in count data ####
# add count data for both elephants to calculate times observed apart
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  janitor::clean_names()
males$id <- paste0('M',males$casename)
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA
males$p7 <- NA ; males$p8 <- NA ; males$p9 <- NA ; males$p10 <- NA ; males$p11 <- NA ; males$p12 <- NA
males$p13 <- NA ; males$p14 <- NA ; males$p15 <- NA ; males$p16 <- NA ; males$p17 <- NA ; males$p18 <- NA
males$p19 <- NA ; males$p20 <- NA ; males$p21 <- NA ; males$p22 <- NA ; males$p23 <- NA ; males$p24 <- NA
males$p25 <- NA ; males$p26 <- NA ; males$p27 <- NA ; males$p28 <- NA ; males$p29 <- NA ; males$p30 <- NA
males$p31 <- NA ; males$p32 <- NA ; males$p33 <- NA ; males$p34 <- NA ; males$p35 <- NA ; males$p36 <- NA

for(i in 1:nrow(males)){
  counts <- ate[ate$id == males$id[i],]
  males$count_all[i] <- nrow(counts)
  for(j in 1:nrow(periods)){
    males[i,(j+22)] <- length(which(counts$obs_date >= periods$period_start[j] & counts$obs_date < periods$period_end[j]))
  }
}
rm(counts) ; gc()

males_long <- pivot_longer(males, cols = 23:58, names_to = 'period', values_to = 'count_period') %>% 
  separate(period, into = c('p','period'), sep = 1, remove = T) %>% 
  select(casename, id, period, count_period) %>% 
  mutate(id_1 = id, id_2 = id,
         count_period_1 = count_period, count_period_2 = count_period,
         period = as.numeric(period))
cdf <- left_join(counts_df, males_long[,c(5,3,7)], by = c('period','id_1'))
cdf <- left_join(cdf, males_long[,c(6,3,8)], by = c('period','id_2'))

str(cdf)
summary(cdf$count_period_1)
summary(cdf$count_period_2)

head(cdf)
which(is.na(cdf$count_period_1))

cdf$count_period_dyad <- (cdf$count_period_1 + cdf$count_period_2) - cdf$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

length(which(cdf$event_count > cdf$count_period_dyad)) # 17 coming up with more events together than total of times individuals observed
#check <- cdf[which(cdf$event_count > cdf$count_period_dyad),
#             c('dyad_id',"id_1","id_2","node_1", "node_2","period","event_count",
#               "count_period_1","count_period_2","count_period_dyad")]
#table(check$id_1) ; table(check$id_2) # a few of the same ones

# delete dyads where one or other is not seen at any point during a time window -- if never seen, you don't know that they were actually within the population at all for that time, and therefore they should be treated as if you don't know they exist
table(cdf$count_period_1)
table(cdf$count_period_2)
cdf_non0 <- cdf[cdf$count_period_1 > 0,]           # 1866589 to 1198594
table(cdf_non0$count_period_2) # 762470 where one never seen
cdf_non0 <- cdf_non0[cdf_non0$count_period_2 > 0,] # 1198594 to 436124

# add column for total number of sightings per pair where they were NOT together
cdf_non0$apart <- cdf_non0$count_period_dyad - cdf_non0$event_count

# create variable for age difference (no point doing the same for sex difference as all are males)
cdf_non0$age_start_1 <- round(as.numeric(as.character((cdf_non0$period_start - cdf_non0$birth_1)/365.25)),0) # approx. age id_1 at start of period
cdf_non0$age_start_2 <- round(as.numeric(as.character((cdf_non0$period_start - cdf_non0$birth_2)/365.25)),0) # approx. age id_2 at start of period
cdf_non0$age_diff <- abs(cdf_non0$age_start_1 - cdf_non0$age_start_2)

### write new data frame to file including counts ####
write_delim(cdf_non0, file = '../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv', delim = ',')

## clean environment
rm(list = ls()[! ls() %in% c('ids', 'ate', 'sightings')])

############### produce data for 7 time windows (longer window) -- can't just combine 5 windows into 1 as already deleted dyads where an individual was not seen or died in the 500 day window  -- HAVE NOW CORRECTED EVERYTHING ABOVE SO RUNNING FROM HERE SHOULD WORK TO PRODUCE SOMETHING WHERE EVENT COUNT IS NEVER HIGHER THAN DURATION ####
## read in processed data files
anp <- readRDS(file = '../data_processed/anp_allpairwiseevents/anp_bayesian_allpairwiseevents.RDS')

## merge sightings information into dyad data
colnames(anp)[4] <- 'obs_id_intfact'
anp <- left_join(x = anp, y = sightings, by = 'obs_id_intfact') %>% distinct()

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
anp$dyad <- paste(anp$node_1, anp$node_2, sep = '_')
anp$dyad_id <- as.integer(as.factor(anp$dyad))
anp$location_id <- as.integer(as.factor(anp$grid_code))
head(anp)

max(sightings$obs_date) - min(sightings$obs_date) # 18197 days -- split into 7 time windows
num_periods <- 7
periods <- data.frame(period = 1:num_periods,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = (num_periods+1))[1:num_periods],
                      period_end = seq(from = min(sightings$obs_date),
                                       to = max(sightings$obs_date),
                                       length.out = (num_periods+1))[2:(num_periods+1)])
periods$period_start <- as.Date(periods$period_start)
periods$period_end <- as.Date(periods$period_end)

events <- anp[anp$social_event == 1,] ; events <- events[!is.na(events$node_1),]
events$period <- NA
for(i in 1:nrow(events)){
  x <- which(periods$period_start <= events$obs_date[i]) # vector of which periods started before observation
  events$period[i] <- x[length(x)] # take last value in vector (most recent = current)
  if(i %% 1000 == 0) {print(i)}
}
periods$period_start ; View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

### convert to Binomial model data format -- aggregate all sightings of each dyad together into a count ####
df_split <- events %>%
  group_by(node_1, node_2, period) %>%
  summarise(event_count = sum(social_event),
            dyad_id = cur_group_id())
head(df_split)

eles <- ate[,c('id','node_id')] %>% distinct()
colnames(eles) <- c('id_1','node_1')
df <- left_join(df_split, eles, by = 'node_1') %>% distinct()
colnames(eles) <- c('id_2','node_2')
df <- left_join(df, eles, by = 'node_2') %>% distinct()
head(df) ; tail(df)

df$dyad <- paste(df$id_1, df$id_2, sep = '_')
#df$dyad_id_period <- df$dyad_id              # every dyad has it's own ID number, including if same dyad in a different time window
#df$dyad_id <- as.integer(as.factor(df$dyad)) # every dyad has it's own ID number, but same dyad in different windows share ID number

### create dyad row for all pairs per period ####
colnames(eles) <- c('id','node') # unnecessary, but makes it clearer that these are not part of a dyad but rather all of the elephants in the population
dyads <- data.frame(id_1 = rep(sort(eles$id), each = nrow(eles)),
                    id_2 = rep(sort(eles$id), nrow(eles)),
                    node_1 = rep(sort(eles$node), each = nrow(eles)),
                    node_2 = rep(sort(eles$node), nrow(eles))
)

# duplicate data frame for every period
dyads <- data.frame(id_1 = rep(dyads$id_1, num_periods),
                    id_2 = rep(dyads$id_2, num_periods),
                    node_1 = rep(dyads$node_1, num_periods),
                    node_2 = rep(dyads$node_2, num_periods),
                    period = rep(sort(unique(df$period)), each = nrow(dyads)))
head(df) ; head(dyads)

#rm(events, df_split) # may need to clear environment a little before running next step
data <- left_join(x = dyads, y = df, by = c('id_1','id_2','period'))
data <- data[,c(1:5,8:10)]
colnames(data)[3:4] <- c('node_1','node_2')
data$event_count <- ifelse(is.na(data$event_count) == TRUE, 0, data$event_count)
table(data$event_count)
data$dyad <- paste(data$id_1, data$id_2, sep = '_')
data$dyad_id <- as.integer(as.factor(data$dyad))
head(data, 20)

data <- left_join(x = data, y = periods, by = 'period')
head(data)

### add data about nodes ####
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% 
  janitor::clean_names()
males$id <- paste('M',males$casename, sep = '')
males <- males[,c(21,2,5,6,8,9,14,18)]
colnames(males)
colnames(data)

colnames(males)[1] <- 'id_1'
data <- left_join(x = data, y = males, by = 'id_1')
colnames(data)[c(11:17)] <- c('id_no_1','bmo_1','byr_1','dmo_1','dyr_1','indyr_1','musthyr_1')
colnames(males)[1] <- 'id_2'
data <- left_join(x = data, y = males, by = 'id_2')
colnames(data)[c(18:24)] <- c('id_no_2','bmo_2','byr_2','dmo_2','dyr_2','indyr_2','musthyr_2')
data <- data[,c(8,7,1:5,9,10,6,11,18,12,19,13,20,14,21,15,22,16,23,17,24)]
head(data)

### write to file ####
data <- data[data$id_1 != data$id_2,]
data <- data[data$node_1 < data$node_2,]

readr::write_delim(data, '../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_longwindows.csv', delim = ',')

## clean environment
rm(list = ls() [ ! ls() %in% c('periods','ate','sightings')])

### remove impossible dyads -- dyads where one individual or the other is dead/unborn before end/at start of period ####
# add column for total number of sightings per pair
counts_df <- read_delim('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_longwindows.csv', delim = ',')
str(counts_df)

# check all dyads present
length(unique(counts_df$id_1)) ; length(unique(counts_df$id_2))
ids <- unique(c(counts_df$id_1, counts_df$id_2))
nrow(counts_df) == (cumsum(1:length(ids))[length(ids)-1])*length(unique(counts_df$period)) # number of dyads per time period should be the cumulative sum of n-1 where n is the number of elephants in the dataset. counts_df contains this for 7 time periods.

# estimate missing birth dates and standardise formatting
head(counts_df)
counts_df$birth_1 <- ifelse(counts_df$bmo_1 > 0 & counts_df$bmo_1 < 10,   # combine birth year and month to single value, set day to 1st of month
                            paste0(counts_df$byr_1,'-0',counts_df$bmo_1,'-01'),
                            paste0(counts_df$byr_1,'-',counts_df$bmo_1,'-01'))
counts_df$birth_1 <- ifelse(is.na(counts_df$bmo_1) == TRUE,   # if birth month unknown, make January of that year
                            paste0(counts_df$byr_1,'-01-01'),
                            counts_df$birth_1)
counts_df$birth_1 <- lubridate::as_date(counts_df$birth_1)    # convert to date format
length(which(is.na(counts_df$birth_1) == TRUE))
summary(counts_df$birth_1)

counts_df$birth_2 <- ifelse(counts_df$bmo_2 > 0 & counts_df$bmo_2 < 10,
                            paste0(counts_df$byr_2,'-0',counts_df$bmo_2,'-01'),
                            paste0(counts_df$byr_2,'-',counts_df$bmo_2,'-01'))
counts_df$birth_2 <- ifelse(is.na(counts_df$bmo_2) == TRUE,
                            paste0(counts_df$byr_2,'-01-01'),
                            counts_df$birth_2)
counts_df$birth_2 <- lubridate::as_date(counts_df$birth_2)
length(which(is.na(counts_df$birth_2) == TRUE))
summary(counts_df$birth_2)

counts_df$death_1 <- ifelse(counts_df$dmo_1 > 0 & counts_df$dmo_1 < 10 & counts_df$dyr_1 > 100,
                            paste0(counts_df$dyr_1,'-0',counts_df$dmo_1,'-01'),
                            ifelse(counts_df$dmo_1 > 0 & counts_df$dmo_1 >= 10 & counts_df$dyr_1 > 100,
                                   paste0(counts_df$dyr_1,'-',counts_df$dmo_1,'-01'),
                                   'living'))
counts_df$death_1 <- ifelse(is.na(counts_df$dmo_1) == TRUE & counts_df$dyr_1 > 100,
                            paste0(counts_df$dyr_1,'-01-01'),
                            counts_df$death_1)
counts_df$death_1[which(counts_df$death_1 == 'living')] <- '2022-08-08' # set date to when raw data were sent over, so definitely after last sighting.
sort(unique(counts_df$death_1))
counts_df$death_1 <- lubridate::as_date(counts_df$death_1)
length(which(is.na(counts_df$death_1) == TRUE))
summary(counts_df$death_1)

counts_df$death_2 <- ifelse(counts_df$dmo_2 > 0 & counts_df$dmo_2 < 10 & counts_df$dyr_2 > 100,
                            paste0(counts_df$dyr_2,'-0',counts_df$dmo_2,'-01'),
                            ifelse(counts_df$dmo_2 > 0 & counts_df$dmo_2 >= 10 & counts_df$dyr_2 > 100,
                                   paste0(counts_df$dyr_2,'-',counts_df$dmo_2,'-01'),
                                   'living'))
counts_df$death_2 <- ifelse(is.na(counts_df$dmo_2) == TRUE & counts_df$dyr_2 > 100,
                            paste0(counts_df$dyr_2,'-01-01'),
                            counts_df$death_2)
counts_df$death_2[which(counts_df$death_2 == 'living')] <- '2022-08-08' # set date to when raw data were sent over, so definitely after last sighting.
counts_df$death_2 <- lubridate::as_date(counts_df$death_2)
length(which(is.na(counts_df$death_2) == TRUE))
summary(counts_df$death_2)

counts_df$alive_1 <- ifelse(counts_df$birth_1 > counts_df$period_start, 'No',
                            ifelse(counts_df$death_1 < counts_df$period_end, 'No', 'Yes'))
length(which(is.na(counts_df$alive_1) == TRUE))
table(counts_df$alive_1)
counts_df$alive_2 <- ifelse(counts_df$birth_2 > counts_df$period_start, 'No',
                            ifelse(counts_df$death_2 < counts_df$period_end, 'No', 'Yes'))
length(which(is.na(counts_df$alive_2) == TRUE))
table(counts_df$alive_2)

counts_df <- counts_df[counts_df$alive_1 == 'Yes' & counts_df$alive_2 == 'Yes', 1:28] ; counts_df <- counts_df[!is.na(counts_df$dyad),]

### add in count data ####
# add count data for both elephants to calculate times observed apart
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA ; males$p7 <- NA
for(i in 1:nrow(males)){
  counts <- ate[ate$id == males$id[i],]
  males$count_all[i] <- nrow(counts)
  for(j in 1:nrow(periods)){
    males[i,(j+22)] <- length(which(counts$obs_date >= periods$period_start[j] & counts$obs_date < periods$period_end[j]))
  }
}
rm(counts) ; gc()

males_long <- pivot_longer(males, cols = 23:29, names_to = 'period', values_to = 'count_period') %>% 
  separate(period, into = c('p','period'), sep = 1, remove = T) %>% 
  select(casename, id, period, count_period) %>% 
  mutate(id_1 = id, id_2 = id,
         count_period_1 = count_period, count_period_2 = count_period,
         period = as.numeric(period))
cdf <- left_join(counts_df, males_long[,c(5,3,7)], by = c('period','id_1'))
cdf <- left_join(cdf, males_long[,c(6,3,8)], by = c('period','id_2'))

str(cdf)
summary(cdf$count_period_1)
summary(cdf$count_period_2)

head(cdf)
which(is.na(cdf$count_period_1))

cdf$count_period_dyad <- (cdf$count_period_1 + cdf$count_period_2) - cdf$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

length(which(cdf$event_count > cdf$count_period_dyad))

# delete dyads where one or other is not seen at any point during a time window -- if never seen, you don't know that they were actually within the population at all for that time, and therefore they should be treated as if you don't know they exist
table(cdf$count_period_1)
table(cdf$count_period_2)
cdf_non0 <- cdf[cdf$count_period_1 > 0,]
table(cdf_non0$count_period_2)
cdf_non0 <- cdf_non0[cdf_non0$count_period_2 > 0,]

# add column for total number of sightings per pair where they were NOT together
cdf_non0$apart <- cdf_non0$count_period_dyad - cdf_non0$event_count

# create variable for age difference (no point doing the same for sex difference as all are males)
cdf_non0$age_start_1 <- round(as.numeric(as.character((cdf_non0$period_start - cdf_non0$birth_1)/365.25)),0) # approx. age id_1 at start of period
cdf_non0$age_start_2 <- round(as.numeric(as.character((cdf_non0$period_start - cdf_non0$birth_2)/365.25)),0) # approx. age id_2 at start of period
cdf_non0$age_diff <- abs(cdf_non0$age_start_1 - cdf_non0$age_start_2)

### write new data frame to file including counts ####
write_delim(cdf_non0, file = '../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_longperiods_impossiblepairsremoved.csv', delim = ',')

## clean environment
rm(list = ls())
