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

load('ate_dataprocessing1_preparallelisation.RData')

############### produce data for 500 day time windows (match MOTNP) ####
### read in data ####
## read in processed data files
anp <- readRDS('../data_processed/anp_bayesian_allpairwiseevents_block1.RDS') %>% as.data.frame()

# check that data processing has worked properly
block_test <- data.frame(block = 1:488, check = NA)
for(i in 2:488){
  # read in block data RDS and append to rest of data
  new_block <- readRDS(paste0('../data_processed/anp_bayesian_allpairwiseevents_block',i,'.RDS')) %>% as.data.frame()
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
saveRDS(anp, '../data_processed/anp_bayesian_allpairwiseevents.RDS')
#anp <- readRDS(file = '../data_processed/anp_bayesian_allpairwiseevents.RDS')

## merge sightings information into dyad data
sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_old_std','grid_code')] %>% distinct()
sightings$obs_id_intfact <- as.integer(as.factor(sightings$obs_id))
nrow(sightings)              # 24386
length(unique(ate$obs_id))   # 24386

colnames(anp)[4] <- 'obs_id_intfact'
anp <- left_join(x = anp, y = sightings, by = 'obs_id_intfact') %>% distinct()

# clean environment
rm(ate_asnipe, block_table, block_test, factor_blocks, gbi_check, new_block, x, block_size, factor_levels, i, length_na)

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
anp$dyad <- paste(anp$node_1, anp$node_2, sep = '_')
anp$dyad_id <- as.integer(as.factor(anp$dyad))
anp$location_id <- as.integer(as.factor(anp$grid_code))
head(anp)

max(sightings$obs_date) - min(sightings$obs_date) # 18197 days. Time period for ALERT data = 504 days --> split into 36 time windows (so need 37 boundaries)
periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 37)
events <- anp[anp$social_event == 1,] ; events <- events[!is.na(events$node_1),] # cuts down to 114188 dyad pairs from 38098504 -- assigning to time window doesn't take forever
events$period <- NA
for(i in 1:nrow(events)){
  events$period[i] <- which(periods <= events$obs_date[i])[length(which(periods <= events$obs_date[i]))] # take last value in vector
  if(i %% 1000 == 0) {print(i)}
}
periods ; View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

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

rm(test3, test5, test18)

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

periods <- data.frame(period_start = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 37)[1:36],
                      period = 1:36)
data <- left_join(x = data, y = periods, by = 'period')
head(data)

### add data about nodes ####
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
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
readr::write_delim(data, '../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows.csv', delim = ',')

### remove impossible dyads -- dyads where one individual or the other is dead/unborn before end/at start of period ####
# add column for total number of sightings per pair
ate <- readxl::read_excel(path = '../data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1) %>%
  janitor::clean_names()
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

# import data for aggregated model (binomial)
#counts_df <- read_delim('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows.csv', delim = ',') %>% 
select(-c('indyr_1','indyr_2','musthyr_1','musthyr_2'))
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

periods <- data.frame(period = 1:36,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = 37)[1:36],
                      period_end = seq(from = min(sightings$obs_date),
                                       to = max(sightings$obs_date),
                                       length.out = 37)[2:37])

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

### add in count data ####
# add count data for both elephants to calculate times observed apart
#males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA ; males$p7 <- NA
males$p8 <- NA ; males$p9 <- NA ; males$p10 <- NA ; males$p11 <- NA ; males$p12 <- NA ; males$p13 <- NA
males$p14 <- NA ; males$p15 <- NA ; males$p16 <- NA ; males$p17 <- NA ; males$p18 <- NA ; males$p19 <- NA
males$p20 <- NA ; males$p21 <- NA ; males$p22 <- NA ; males$p23 <- NA ; males$p24 <- NA ; males$p25 <- NA
males$p26 <- NA ; males$p27 <- NA ; males$p28 <- NA ; males$p29 <- NA ; males$p30 <- NA ; males$p31 <- NA
males$p32 <- NA ; males$p33 <- NA ; males$p34 <- NA ; males$p35 <- NA ; males$p36 <- NA ; males$p37 <- NA
periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 37)
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
  males$p32[i] <- length(which(counts$obs_date >= periods[32] & counts$obs_date < periods[33]))
  males$p33[i] <- length(which(counts$obs_date >= periods[33] & counts$obs_date < periods[34]))
  males$p34[i] <- length(which(counts$obs_date >= periods[34] & counts$obs_date < periods[35]))
  males$p35[i] <- length(which(counts$obs_date >= periods[35] & counts$obs_date < periods[36]))
  males$p36[i] <- length(which(counts$obs_date >= periods[36] & counts$obs_date < periods[37]))
}

# counts per dyad per time period
periods <- data.frame(period = 1:36,
                      period_start = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 37)[1:36],
                      period_end   = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 37)[2:37])
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

head(counts_df)
which(is.na(counts_df$period_count_1))[1]

counts_df$period_count_dyad <- (counts_df$period_count_1 + counts_df$period_count_2) - counts_df$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

# delete dyads where one or other is not seen at any point during a time window -- if never seen, you don't know that they were actually within the population at all for that time, and therefore they should be treated as if you don't know they exist
table(counts_df$period_count_1)
table(counts_df$period_count_2)
counts_df_non0 <- counts_df[counts_df$period_count_1 > 0,]           # 1866589 to 1198594
table(counts_df_non0$period_count_2) # 762470 where one never seen
counts_df_non0 <- counts_df_non0[counts_df_non0$period_count_2 > 0,] # 1198594 to 436124

counts_df_non0$id_1[which(counts_df_non0$period_count_1 == 76)][1]   # M022
test <- counts_df_non0[counts_df_non0$period_count_1 == 76,]
test_raw <- ate[ate$obs_date > test$period_start[1] & ate$obs_date < test$period_end[1],]
length(which(test_raw$id == 'M022')) # 76
test$period_end[1] - test$period_start[1]

# add column for total number of sightings per pair where they were NOT together
counts_df_non0$apart <- counts_df_non0$period_count_dyad - counts_df_non0$event_count

# create variable for age difference (no point doing the same for sex difference as all are males)
counts_df_non0$age_start_1 <- round(as.numeric(as.character((counts_df_non0$period_start - counts_df_non0$birth_1)/365.25)),0) # approx. age id_1 at start of period
counts_df_non0$age_start_2 <- round(as.numeric(as.character((counts_df_non0$period_start - counts_df_non0$birth_2)/365.25)),0) # approx. age id_2 at start of period
counts_df_non0$age_diff <- abs(counts_df_non0$age_start_1 - counts_df_non0$age_start_2)

# add back in death2 variable as I seem to have deleted it at some point
counts_df_non0$death_2 <- ifelse(counts_df_non0$dmo_2 > 0 & counts_df_non0$dmo_2 < 10 & counts_df_non0$dyr_2 > 100,
                                 paste0(counts_df_non0$dyr_2,'-0',counts_df_non0$dmo_2,'-01'),
                                 ifelse(counts_df_non0$dmo_2 > 0 & counts_df_non0$dmo_2 >= 10 & counts_df_non0$dyr_2 > 100,
                                        paste0(counts_df_non0$dyr_2,'-',counts_df_non0$dmo_2,'-01'),
                                        'living'))
counts_df_non0$death_2 <- ifelse(is.na(counts_df_non0$dmo_2) == TRUE & counts_df_non0$dyr_2 > 100,
                                 paste0(counts_df_non0$dyr_2,'-01-01'),
                                 counts_df_non0$death_2)
counts_df_non0$death_2[which(counts_df_non0$death_2 == 'living')] <- '2022-08-08' # set date to when raw data were sent over, so definitely after last sighting.
counts_df_non0$death_2 <- lubridate::as_date(counts_df_non0$death_2)
length(which(is.na(counts_df_non0$death_2) == TRUE))
summary(counts_df_non0$death_2)

### write new data frame to file including counts ####
write_delim(counts_df_non0, file = '../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_shortwindows_impossiblepairsremoved.csv', delim = ',')

## clean environment
rm(list = ls())

############### produce data for 7 time windows (longer window) -- can't just combine 5 windows into 1 as already deleted dyads where an individual was not seen or died in the 500 day window ####
## read in processed data files
anp <- readRDS(file = '../data_processed/anp_bayesian_allpairwiseevents.RDS')

## merge sightings information into dyad data
sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_old_std','grid_code')] %>% distinct()
sightings$obs_id_intfact <- as.integer(as.factor(sightings$obs_id))
nrow(sightings)              # 24386
length(unique(ate$obs_id))   # 24386

colnames(anp)[4] <- 'obs_id_intfact'
anp <- left_join(x = anp, y = sightings, by = 'obs_id_intfact') %>% distinct()

# clean environment
rm(ate_asnipe, block_table, factor_blocks, block_size, factor_levels, i, length_na)

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
anp$dyad <- paste(anp$node_1, anp$node_2, sep = '_')
anp$dyad_id <- as.integer(as.factor(anp$dyad))
anp$location_id <- as.integer(as.factor(anp$grid_code))
head(anp)

max(sightings$obs_date) - min(sightings$obs_date) # 18197 days -- split into 7 time windows
periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 8)
events <- anp[anp$social_event == 1,] ; events <- events[!is.na(events$node_1),] # cuts down to 114188 dyad pairs from 38098504 -- assigning to time window doesn't take forever
events$period <- NA
for(i in 1:nrow(events)){
  events$period[i] <- which(periods <= events$obs_date[i])[length(which(periods <= events$obs_date[i]))] # take last value in vector
  if(i %% 1000 == 0) {print(i)}
}
periods ; View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

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

#rm(events, df_split) # may need to clear environment a little before running next step
data <- left_join(x = dyads, y = df, by = c('id_1','id_2','period'))
data <- data[,c(1:5,8:10)]
colnames(data)[3:4] <- c('node_1','node_2')
data$event_count <- ifelse(is.na(data$event_count) == TRUE, 0, data$event_count)
table(data$event_count)
data$dyad <- paste(data$id_1, data$id_2, sep = '_')
data$dyad_id <- as.integer(as.factor(data$dyad))
head(data, 20)

periods <- data.frame(period_start = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 8)[1:7],
                      period = 1:7)
data <- left_join(x = data, y = periods, by = 'period')
head(data)

### add data about nodes ####
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
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
readr::write_delim(data, '../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_longwindows.csv', delim = ',')

## clean environment
rm(list = ls())

### remove impossible dyads -- dyads where one individual or the other is dead/unborn before end/at start of period ####
# add column for total number of sightings per pair
ate <- readxl::read_excel(path = '../data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1) %>%
  janitor::clean_names()
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

# import data for aggregated model (binomial)
#counts_df <- read_delim('../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_longwindows.csv', delim = ',') %>% 
select(-c('indyr_1','indyr_2','musthyr_1','musthyr_2'))
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

periods <- data.frame(period = 1:7,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = 8)[1:7],
                      period_end = seq(from = min(sightings$obs_date),
                                       to = max(sightings$obs_date),
                                       length.out = 8)[2:8])

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

counts_df <- counts_df[counts_df$alive_1 == 'Yes' & counts_df$alive_2 == 'Yes', 1:24] ; counts_df <- counts_df[!is.na(counts_df$dyad),]

### add in count data ####
# add count data for both elephants to calculate times observed apart
males <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
males$id <- paste0('M',males$casename)
males <- males %>% dplyr::filter(id %in% ids)

# counts per male per time period
males$count_all <- NA
males$p1 <- NA ; males$p2 <- NA ; males$p3 <- NA ; males$p4 <- NA ; males$p5 <- NA ; males$p6 <- NA ; males$p7 <- NA
periods <- seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 8)
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
}

# counts per dyad per time period
periods <- data.frame(period = 1:7,
                      period_start = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 8)[1:7],
                      period_end   = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 8)[2:8])
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

head(counts_df)
which(is.na(counts_df$period_count_1))[1]

counts_df$period_count_dyad <- (counts_df$period_count_1 + counts_df$period_count_2) - counts_df$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

# delete dyads where one or other is not seen at any point during a time window -- if never seen, you don't know that they were actually within the population at all for that time, and therefore they should be treated as if you don't know they exist
table(counts_df$period_count_1)
table(counts_df$period_count_2)
counts_df_non0 <- counts_df[counts_df$period_count_1 > 0,]           # 271486 to 236533
table(counts_df_non0$period_count_2) # 762470 where one never seen
counts_df_non0 <- counts_df_non0[counts_df_non0$period_count_2 > 0,] # 236533 to 154807

# add column for total number of sightings per pair where they were NOT together
counts_df_non0$apart <- counts_df_non0$period_count_dyad - counts_df_non0$event_count

# create variable for age difference (no point doing the same for sex difference as all are males)
counts_df_non0$age_start_1 <- round(as.numeric(as.character((counts_df_non0$period_start - counts_df_non0$birth_1)/365.25)),0) # approx. age id_1 at start of period
counts_df_non0$age_start_2 <- round(as.numeric(as.character((counts_df_non0$period_start - counts_df_non0$birth_2)/365.25)),0) # approx. age id_2 at start of period
counts_df_non0$age_diff <- abs(counts_df_non0$age_start_1 - counts_df_non0$age_start_2)

### write new data frame to file including counts ####
write_delim(counts_df_non0, file = '../data_processed/anp_bayesian_pairwiseevents_aggregated_allperiods_longperiods_impossiblepairsremoved.csv', delim = ',')

## clean environment
rm(list = ls())
