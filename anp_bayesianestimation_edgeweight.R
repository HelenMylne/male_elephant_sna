#### Bayesian analysis of ATE data ####
# Script to process association data from Amboseli National Park, Kenya
# Data collected: 1972-2021 by Amboseli Trust for Elephants
# Data supplied by: Vicki Fishlock (March 2022) and Phyllis Lee (December 2021)
# Data input: raw data provided by ATE processed in scripts 22.03.20_anp_dataprocessing1.R and 22.03.22_anp_dataprocessing2.R

#### Set up ####
# load packages
library(tidyverse)
library(dplyr)
library(rstan)
library(rethinking)
library(igraph)
library(dagitty)
library(cmdstanr)

# information
sessionInfo()
R.Version()
rstan::stan_version()
#packageVersion('')
#citation('')

# set stan path
set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')

# load model
mod_2.2 <- cmdstan_model("models/edge_weight_estimation/simpleBetaNet_HKM_2.2_22.02.03.stan")
mod_2.2

# set seed
set.seed(12345)

################ 1) Draw DAGS ################
# plot with full names
binom <- dagitty::dagitty("dag{
                         age_dyad [exposure];
                         sex_dyad [exposure];
                         weight [outcome];
                         relatedness_dyad [unobserved];
                         age_dyad -> weight <- sex_dyad;
                         weight <- relatedness_dyad;
                         }")
dagitty::coordinates(binom) <- list(x = c(age_dyad = 0, sex_dyad = 1, weight = 2, relatedness = 3),
                                    y = c(age_dyad = 0, sex_dyad = 2, weight = 1, relatedness = 0))
drawdag(binom)

# plot with letters
binom <- dagitty::dagitty("dag{
                         A [exposure];
                         S [exposure];
                         W [outcome];
                         R [unobserved];
                         A -> W <- S;
                         W <- R;
                         }")
dagitty::coordinates(binom) <- list(x = c(A = 0.5, S = 1, W = 1.0, R = 1.5),
                                    y = c(A = 0.5, S = 2, W = 1.2, R = 0.5))
drawdag(binom, radius = 6, cex = 1.6)

# clear environment and reset plot window
rm(binom)
dev.off()
################ 2) Create data lists ################
# long version (put into data processing script): -- do not rerun unnecessarily, have commented out data read-in lines ####
### add column for total number of sightings per pair
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
head(ate)
rm(date_row, ate_nums, lu, i)

sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_std','grid_code')] %>% distinct()
sightings <- sightings[c(1:12,15:24176),]

### import data for aggregated model (binomial)
#counts_df <- read_delim('data_processed/anp_bayesian_pairwiseevents_22.03.28.csv', delim = ',') %>% 
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
#males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
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
write_csv(counts_df_non0, path = 'data_processed/anp_bayesian_pairwiseevents_22.04.06.csv')

# clean environment
rm(counts_1, counts_2, counts, counts_df, ate, periods, sightings, i, ids)

# short version: ####
counts_df_non0 <- read_csv('data_processed/anp_bayesian_pairwiseevents_22.04.06.csv')

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
head(ate)
rm(date_row, ate_nums, lu, i)

sightings <- ate[,c('obs_id','obs_date','correct_time_hms','obs_num_std','grid_code')] %>% distinct()
sightings <- sightings[c(1:12,15:24176),]

males <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()
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

# time windows
periods <- data.frame(period = 1:31,
                      period_start = seq(from = min(sightings$obs_date),
                                         to = max(sightings$obs_date),
                                         length.out = 32)[1:31],
                      period_end = seq(from = min(sightings$obs_date),
                                       to = max(sightings$obs_date),
                                       length.out = 32)[2:32])

# check out total sightings per individual/dyad per period and see if they are reasonable ####
table(counts_df_non0$period_count_1) ; hist(counts_df_non0$period_count_1)
table(counts_df_non0$period_count_2) ; hist(counts_df_non0$period_count_2)
table(counts_df_non0$period_count_dyad) ; hist(counts_df_non0$period_count_dyad, breaks = 30)
table(counts_df_non0$event_count) ; hist(counts_df_non0$event_count, breaks = 3)
# many are not seen much, but repeat sightings of pairs do seem to be genuinely rare

################ 6) Run model on real standardised data -- all ################
### create data list -- can contain no NA values in any column, even if column is not specified in model -- NEED TO SEPARATE BY TIME PERIOD OR SOMETHING SOMEHOW
counts_ls <- list(
  n_dyads  = nrow(counts_df_non0),          # total number of times one or other of the dyad was observed
  together = counts_df_non0$event_count,    # count number of sightings seen together
  apart    = counts_df_non0$apart,          # count number of sightings seen apart
  period   = counts_df_non0$period)         # which period it's within

### Fit model
weight_anp_2.2 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2
# variable        mean      median     sd    mad          q5         q95 rhat ess_bulk ess_tail
#lp__      -3125125.08 -3125130.00 480.56 489.26 -3125910.00 -3124320.00 1.00     1247     2132
#weight[1]        0.18        0.16   0.11   0.10        0.04        0.39 1.00     5295     2343
#weight[2]        0.28        0.26   0.16   0.16        0.07        0.57 1.00     5887     2562
#weight[3]        0.29        0.27   0.16   0.17        0.06        0.57 1.00     4645     2136
#weight[4]        0.17        0.15   0.10   0.10        0.03        0.36 1.00     4791     2023
#weight[5]        0.25        0.23   0.14   0.15        0.06        0.52 1.00     5203     2650
#weight[6]        0.22        0.20   0.13   0.13        0.05        0.46 1.00     3960     2120
#weight[7]        0.29        0.27   0.16   0.16        0.07        0.58 1.00     4833     1618
#weight[8]        0.28        0.26   0.16   0.17        0.06        0.58 1.00     4751     2500
#weight[9]        0.29        0.26   0.16   0.16        0.06        0.58 1.00     4743     2359
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2$output_files()[1])
draws1_anp2.2 <- as.data.frame(output1$post_warmup_draws)
rm(output1)
write_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_chain1_22.04.06.csv', x = draws1_anp2.2)

output2 <- read_cmdstan_csv(weight_anp_2.2$output_files()[2])[,1:50]
draws2_anp2.2 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2$output_files()[3])
draws3_anp2.2 <- as.data.frame(output3$post_warmup_draws)
r(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2$output_files()[4])
draws4_anp2.2 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp2.2 <- rbind(draws1_anp2.2, draws2_anp2.2, draws3_anp2.2, draws4_anp2.2)

colnames(draws_anp2.2)[2:402699] <- counts_df_non0$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:402699, size = 20, replace = F)

# tidy data -- vector memory usually exhausted and won't run
tidy_draws_2.2 <- pivot_longer(draws_anp2.2[,2:402699], cols = everything(), names_to = 'dyad', values_to = 'draw')
tidy_draws_2.2$chain <- rep(1:4, each = 402698000)
tidy_draws_2.2$index <- rep(rep(1:1000, each = 402698),4)
head(tidy_draws_2.2, 10)
tail(tidy_draws_2.2, 10)

### save data 
write_csv(draws_anp2.2, 'data_processed/anp_bayesian_edgedistributions_a2.b2_22.04.05.csv')

################ check single chain ################
draws1_anp2.2 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_chain1_22.04.06.csv')
plot_cols <- sample(2:402699, size = 30, replace = F)

### build traceplots -- some quite wide, generally not bad ####
plot(draws1_anp2.2[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1_anp2.2[,plot_cols[2]], col = 'tan')
lines(draws1_anp2.2[,plot_cols[3]], col = 'orange')
lines(draws1_anp2.2[,plot_cols[4]], col = 'green')
lines(draws1_anp2.2[,plot_cols[5]], col = 'chocolate')
lines(draws1_anp2.2[,plot_cols[6]], col = 'blue')
lines(draws1_anp2.2[,plot_cols[7]], col = 'red')
lines(draws1_anp2.2[,plot_cols[8]], col = 'seagreen')
lines(draws1_anp2.2[,plot_cols[9]], col = 'purple')
lines(draws1_anp2.2[,plot_cols[10]],col = 'magenta')
lines(draws1_anp2.2[,plot_cols[11]],col = 'black')
lines(draws1_anp2.2[,plot_cols[12]], col = 'tan')
lines(draws1_anp2.2[,plot_cols[13]], col = 'orange')
lines(draws1_anp2.2[,plot_cols[14]], col = 'green')
lines(draws1_anp2.2[,plot_cols[15]], col = 'chocolate')
lines(draws1_anp2.2[,plot_cols[16]], col = 'blue')
lines(draws1_anp2.2[,plot_cols[17]], col = 'red')
lines(draws1_anp2.2[,plot_cols[18]], col = 'seagreen')
lines(draws1_anp2.2[,plot_cols[19]], col = 'purple')
lines(draws1_anp2.2[,plot_cols[20]],col = 'magenta')
lines(draws1_anp2.2[,plot_cols[21]],col = 'black')
lines(draws1_anp2.2[,plot_cols[22]], col = 'tan')
lines(draws1_anp2.2[,plot_cols[23]], col = 'orange')
lines(draws1_anp2.2[,plot_cols[24]], col = 'green')
lines(draws1_anp2.2[,plot_cols[25]], col = 'chocolate')
lines(draws1_anp2.2[,plot_cols[26]], col = 'blue')
lines(draws1_anp2.2[,plot_cols[27]], col = 'red')
lines(draws1_anp2.2[,plot_cols[28]], col = 'seagreen')
lines(draws1_anp2.2[,plot_cols[29]], col = 'purple')
lines(draws1_anp2.2[,plot_cols[30]],col = 'magenta')

### density plots -- most are very wide, high uncertainty####
dens(draws1_anp2.2[,2], ylim = c(0,100),xlim = c(0,1), las = 1)
for(i in 3:8000){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws1_anp2.2[,i], col = col.alpha('blue', alpha = 0.01))
}



### summarise data ####
# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws1_anp2.2[,2:402699]),
                        min = rep(NA, ncol(draws1_anp2.2)-1),
                        max = rep(NA, ncol(draws1_anp2.2)-1),
                        mean = rep(NA, ncol(draws1_anp2.2)-1),
                        median = rep(NA, ncol(draws1_anp2.2)-1),
                        sd = rep(NA, ncol(draws1_anp2.2)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws1_anp2.2[,i+1])
  summaries$max[i]    <- max(draws1_anp2.2[,i+1])
  summaries$mean[i]   <- mean(draws1_anp2.2[,i+1])
  summaries$median[i] <- median(draws1_anp2.2[,i+1])
  summaries$sd[i]     <- sd(draws1_anp2.2[,i+1])
}

summary(summaries) # generally higher than MOTNP but SD much higher

# organise dem_class -- report age category based on age at start of period
plot_data_anp2.2 <- left_join(x = summaries, y = counts_df_non0, by = 'dyad')
head(plot_data_anp2.2)
plot_data_anp2.2$age_cat_1 <- ifelse(plot_data_anp2.2$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp2.2$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp2.2$age_start_1 > 19, 'A','P')))
plot_data_anp2.2$age_cat_2 <- ifelse(plot_data_anp2.2$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp2.2$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp2.2$age_start_2 > 19, 'A','P')))
plot_data_anp2.2$age_catnum_1 <- ifelse(plot_data_anp2.2$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp2.2$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp2.2$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp2.2$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp2.2$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp2.2$age_start_1 < 40, 6, 7))))))
plot_data_anp2.2$age_catnum_2 <- ifelse(plot_data_anp2.2$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp2.2$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp2.2$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp2.2$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp2.2$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp2.2$age_start_2 < 40, 6, 7))))))

plot_data_anp2.2$age_dyad <- ifelse(plot_data_anp2.2$age_catnum_1 >= plot_data_anp2.2$age_catnum_2,
                                    paste(plot_data_anp2.2$age_cat_1, plot_data_anp2.2$age_cat_2, sep = ''),
                                    paste(plot_data_anp2.2$age_cat_2, plot_data_anp2.2$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp2.2, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp2.2, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries)
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_anp2.2[plot_data_anp2.2$age_dyad == 'AA' | plot_data_anp2.2$age_dyad == 'AP' | plot_data_anp2.2$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(counts_ls)
rm(draws1_anp2.2, draws2_anp2.2, draws3_anp2.2, draws4_anp2.2, tidy_draws1_2.2)
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- single chain ################
head(summaries)
length(unique(plot_data_anp2.2$id_1))+1 # number of individuals = 55

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_non0$id_1))+1,
                         NROW(unique(counts_df_non0$id_2))+1,
                         NROW(draws1_anp2.2)),
                    dimnames = list(c(unique(counts_df_non0$id_1),'M119'),
                                    c('M2',unique(counts_df_non0$id_2)),
                                    NULL))
N <- nrow(counts_df_non0)

for (i in 1:N) {
  dyad_row <- counts_df_non0[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws1_anp2.2[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
males <- males[,c(21,6,9,22:53)]
View(males) # all ones born before start of study

# create variables for different degrees of node connectedness
males$degree_0.1 <- NA
males$degree_0.2 <- NA
males$degree_0.3 <- NA
males$degree_0.4 <- NA
males$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males)){
  rows <- summaries[summaries$id_1 == males$id[i] | summaries$id_2 == males$id[i],]
  males$degree_0.1[i] <- length(which(rows$median > 0.1))
  males$degree_0.2[i] <- length(which(rows$median > 0.2))
  males$degree_0.3[i] <- length(which(rows$median > 0.3))
  males$degree_0.4[i] <- length(which(rows$median > 0.4))
  males$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males$degree_0.1 < males$degree_0.2)
which(males$degree_0.2 < males$degree_0.3)
which(males$degree_0.3 < males$degree_0.4)
which(males$degree_0.4 < males$degree_0.5)

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = NA,
     vertex.size = 5,
     layout = coords)
plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     edge.color = 'black',
     vertex.size = 8,
     vertex.label = males$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     #vertex.color= ifelse(males$age_class == 'Adult','seagreen1',
     #                    ifelse(males$age_class == 'Pubescent','skyblue',
     #                           ifelse(males$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.5,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.5,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     #vertex.color= ifelse(males$age_class == 'Adult','seagreen1',
     #                     ifelse(males$age_class == 'Pubescent','skyblue',
     #                            ifelse(males$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males$id[which(males$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males$id[which(males$degree_0.3 == 0)])

set.seed(2)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     #vertex.color = ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Adult',
     #                      'seagreen1',
     #                      ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Pubescent','skyblue',
     #                             ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Juvenile','yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = log(males$count_all[which(males$degree_0.3 != 0)]),
     vertex.color = 'green',
     vertex.label = males$id[which(males$degree_0.3 != 0)],
     vertex.label.color = 'magenta',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     #vertex.color = ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Adult',
     #                      'seagreen1',
     #                      ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Pubescent','skyblue',
     #                             ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Juvenile','yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

################ repeat for all chains ################
# draws_anp2.2 <- read_csv('data_processed/anp_bayesian_edgedistributions_a2.b2_22.04.05.csv') %>% 
#   data.matrix()
plot_cols <- sample(2:402699, size = 20, replace = F)
### build traceplots -- most are very wide, high uncertainty ####
plot(draws_anp2.2[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp2.2[,plot_cols[2]], col = 'tan')
lines(draws_anp2.2[,plot_cols[3]], col = 'orange')
lines(draws_anp2.2[,plot_cols[4]], col = 'green')
lines(draws_anp2.2[,plot_cols[5]], col = 'chocolate')
lines(draws_anp2.2[,plot_cols[6]], col = 'blue')
lines(draws_anp2.2[,plot_cols[7]], col = 'red')
lines(draws_anp2.2[,plot_cols[8]], col = 'seagreen')
lines(draws_anp2.2[,plot_cols[9]], col = 'purple')
lines(draws_anp2.2[,plot_cols[10]],col = 'magenta')
lines(draws_anp2.2[,plot_cols[11]],col = 'black')
lines(draws_anp2.2[,plot_cols[12]], col = 'tan')
lines(draws_anp2.2[,plot_cols[13]], col = 'orange')
lines(draws_anp2.2[,plot_cols[14]], col = 'green')
lines(draws_anp2.2[,plot_cols[15]], col = 'chocolate')
lines(draws_anp2.2[,plot_cols[16]], col = 'blue')
lines(draws_anp2.2[,plot_cols[17]], col = 'red')
lines(draws_anp2.2[,plot_cols[18]], col = 'seagreen')
lines(draws_anp2.2[,plot_cols[19]], col = 'purple')
lines(draws_anp2.2[,plot_cols[20]],col = 'magenta')

### density plots -- most are very wide, high uncertainty ####
dens(draws_anp2.2[,2], ylim = c(0,10),xlim = c(0,1), las = 1)
for(i in 3:402699){
  dens(add = T, draws_anp2.2[,i], col = col.alpha('blue', alpha = 0.1))
}

### summarise data ####
# summarise -- look for any anomalies in draw values or chain variation
summaries <- data.frame(dyad = colnames(draws_anp2.2[,2:402699]),
                        min = rep(NA, ncol(draws_anp2.2)-1),
                        max = rep(NA, ncol(draws_anp2.2)-1),
                        mean = rep(NA, ncol(draws_anp2.2)-1),
                        median = rep(NA, ncol(draws_anp2.2)-1),
                        sd = rep(NA, ncol(draws_anp2.2)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp2.2[,i+1])
  summaries$max[i]    <- max(draws_anp2.2[,i+1])
  summaries$mean[i]   <- mean(draws_anp2.2[,i+1])
  summaries$median[i] <- median(draws_anp2.2[,i+1])
  summaries$sd[i]     <- sd(draws_anp2.2[,i+1])
}

summary(summaries) # generally higher than MOTNP but SD much higher -- higher mean/median doesn't really make sense given that no pair was ever seen more than twice together in any 2-year period

# organise dem_class -- report age category based on age at start of period
plot_data_anp2.2 <- left_join(x = summaries, y = counts_df_non0, by = 'dyad')
head(plot_data_anp2.2)
plot_data_anp2.2$age_cat_1 <- ifelse(plot_data_anp2.2$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp2.2$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp2.2$age_start_1 > 19, 'A','P')))
plot_data_anp2.2$age_cat_2 <- ifelse(plot_data_anp2.2$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp2.2$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp2.2$age_start_2 > 19, 'A','P')))
plot_data_anp2.2$age_catnum_1 <- ifelse(plot_data_anp2.2$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp2.2$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp2.2$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp2.2$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp2.2$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp2.2$age_start_1 < 40, 6, 7))))))
plot_data_anp2.2$age_catnum_2 <- ifelse(plot_data_anp2.2$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp2.2$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp2.2$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp2.2$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp2.2$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp2.2$age_start_2 < 40, 6, 7))))))

plot_data_anp2.2$age_dyad <- ifelse(plot_data_anp2.2$age_catnum_1 >= plot_data_anp2.2$age_catnum_2,
                                    paste(plot_data_anp2.2$age_cat_1, plot_data_anp2.2$age_cat_2, sep = ''),
                                    paste(plot_data_anp2.2$age_cat_2, plot_data_anp2.2$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp2.2, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp2.2, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries)
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_anp2.2[plot_data_anp2.2$age_dyad == 'AA' | plot_data_anp2.2$age_dyad == 'AP' | plot_data_anp2.2$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(counts_ls)
rm(draws1_anp2.2, draws2_anp2.2, draws3_anp2.2, draws4_anp2.2, tidy_draws_2.2)
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots ################
head(summaries)
length(unique(plot_data_anp2.2$id_1))+1 # number of individuals = 55

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_non0$id_1))+1,
                         NROW(unique(counts_df_non0$id_2))+1,
                         NROW(draws_anp2.2)),
                    dimnames = list(c(unique(counts_df_non0$id_1),'M119'),
                                    c('M2',unique(counts_df_non0$id_2)),
                                    NULL))
N <- nrow(counts_df_non0)

for (i in 1:N) {
  dyad_row <- counts_df_non0[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp2.2[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
males <- males[,c(21,6,9,22:53)]
View(males) # all ones born before start of study

# create variables for different degrees of node connectedness
males$degree_0.1 <- NA
males$degree_0.2 <- NA
males$degree_0.3 <- NA
males$degree_0.4 <- NA
males$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males)){
  rows <- summaries[summaries$id_1 == males$id[i] | summaries$id_2 == males$id[i],]
  males$degree_0.1[i] <- length(which(rows$median > 0.1))
  males$degree_0.2[i] <- length(which(rows$median > 0.2))
  males$degree_0.3[i] <- length(which(rows$median > 0.3))
  males$degree_0.4[i] <- length(which(rows$median > 0.4))
  males$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males$degree_0.1 < males$degree_0.2)
which(males$degree_0.2 < males$degree_0.3)
which(males$degree_0.3 < males$degree_0.4)
which(males$degree_0.4 < males$degree_0.5)

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = NA,
     vertex.size = 5,
     layout = coords)
plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     edge.color = 'black',
     vertex.size = 8,
     vertex.label = males$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     #vertex.color= ifelse(males$age_class == 'Adult','seagreen1',
     #                    ifelse(males$age_class == 'Pubescent','skyblue',
     #                           ifelse(males$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.5,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.5,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     #vertex.color= ifelse(males$age_class == 'Adult','seagreen1',
     #                     ifelse(males$age_class == 'Pubescent','skyblue',
     #                            ifelse(males$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

### only elephants degree > 0.3 ####
g_mid_0.3 <- delete.vertices(graph = g_mid, v = males$id[which(males$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = males$id[which(males$degree_0.3 == 0)])

set.seed(2)
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     #vertex.color = ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Adult',
     #                      'seagreen1',
     #                      ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Pubescent','skyblue',
     #                             ifelse(nodes[which(males$degree_0.3 != 0),]$age_class == 'Juvenile','yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = log(males$count_all[which(males$degree_0.3 != 0)]),
     vertex.color = 'green',
     vertex.label = males$id[which(males$degree_0.3 != 0)],
     vertex.label.color = 'magenta',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     #vertex.color = ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Adult',
     #                      'seagreen1',
     #                      ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Pubescent','skyblue',
     #                             ifelse(males[which(males$degree_0.3 != 0),]$age_class == 'Juvenile','yellow','magenta'))),
     layout = coords_0.3, add = TRUE)



################ 7.1) Run model on real standardised data -- period 1 ################
### create data list
counts_df_period1 <- counts_df_non0[counts_df_non0$period == 1,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period1),          # total number of times one or other of the dyad was observed
  together = counts_df_period1$event_count,    # count number of sightings seen together
  apart    = counts_df_period1$apart,          # count number of sightings seen apart
  period   = counts_df_period1$period)         # which period it's within

### Fit model
weight_anp_2.2_period1 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period1
#variable        mean     median    sd   mad         q5        q95 rhat ess_bulk ess_tail
#lp__      -132933.05 -132930.00 97.62 97.85 -133098.00 -132776.00 1.00     1394     1933
#weight[1]       0.05       0.05  0.04  0.03       0.01       0.12 1.00     3484     1794
#weight[2]       0.05       0.04  0.04  0.03       0.01       0.12 1.00     3514     1903
#weight[3]       0.08       0.07  0.05  0.05       0.01       0.17 1.00     4739     2340
#weight[4]       0.08       0.07  0.05  0.05       0.01       0.18 1.00     4905     2312
#weight[5]       0.18       0.17  0.09  0.09       0.05       0.35 1.00     4790     2435
#weight[6]       0.07       0.06  0.05  0.04       0.01       0.16 1.00     4968     2548
#weight[7]       0.20       0.19  0.09  0.09       0.08       0.36 1.00     5330     2555
#weight[8]       0.09       0.08  0.06  0.05       0.02       0.19 1.00     4911     2139
#weight[9]       0.06       0.05  0.04  0.04       0.01       0.14 1.00     4463     2009
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period1$output_files()[1])
draws1_anp_1 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period1$output_files()[2])
draws2_anp_1 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period1$output_files()[3])
draws3_anp_1 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period1$output_files()[4])
draws4_anp_1 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_1 <- rbind(draws1_anp_1, draws2_anp_1, draws3_anp_1, draws4_anp_1)

colnames(draws_anp_1)[2:ncol(draws_anp_1)] <- counts_df_period1$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_1), size = 30, replace = F)

### save data 
write_csv(draws_anp_1, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period1_22.04.15.csv')
rm(draws1_anp_1, draws2_anp_1, draws3_anp_1, draws4_anp_1)

### build traceplots -- period 1 -- very wide ####
plot(draws_anp_1[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_1[,plot_cols[2]], col = 'tan')
lines(draws_anp_1[,plot_cols[3]], col = 'orange')
lines(draws_anp_1[,plot_cols[4]], col = 'green')
lines(draws_anp_1[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_1[,plot_cols[6]], col = 'blue')
lines(draws_anp_1[,plot_cols[7]], col = 'red')
lines(draws_anp_1[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_1[,plot_cols[9]], col = 'purple')
lines(draws_anp_1[,plot_cols[10]],col = 'magenta')
lines(draws_anp_1[,plot_cols[11]],col = 'black')
lines(draws_anp_1[,plot_cols[12]], col = 'tan')
lines(draws_anp_1[,plot_cols[13]], col = 'orange')
lines(draws_anp_1[,plot_cols[14]], col = 'green')
lines(draws_anp_1[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_1[,plot_cols[16]], col = 'blue')
lines(draws_anp_1[,plot_cols[17]], col = 'red')
lines(draws_anp_1[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_1[,plot_cols[19]], col = 'purple')
lines(draws_anp_1[,plot_cols[20]],col = 'magenta')
lines(draws_anp_1[,plot_cols[21]],col = 'black')
lines(draws_anp_1[,plot_cols[22]], col = 'tan')
lines(draws_anp_1[,plot_cols[23]], col = 'orange')
lines(draws_anp_1[,plot_cols[24]], col = 'green')
lines(draws_anp_1[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_1[,plot_cols[26]], col = 'blue')
lines(draws_anp_1[,plot_cols[27]], col = 'red')
lines(draws_anp_1[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_1[,plot_cols[29]], col = 'purple')
lines(draws_anp_1[,plot_cols[30]],col = 'magenta')

### density plots -- period 1 ####
dens(draws_anp_1[,2], ylim = c(0,10),xlim = c(0,1), las = 1)
for(i in 1:30){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws_anp_1[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 1 ####
# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_1[,2:ncol(draws_anp_1)]),
                        min = rep(NA, ncol(draws_anp_1)-1),
                        max = rep(NA, ncol(draws_anp_1)-1),
                        mean = rep(NA, ncol(draws_anp_1)-1),
                        median = rep(NA, ncol(draws_anp_1)-1),
                        sd = rep(NA, ncol(draws_anp_1)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_1[,i+1])
  summaries$max[i]    <- max(draws_anp_1[,i+1])
  summaries$mean[i]   <- mean(draws_anp_1[,i+1])
  summaries$median[i] <- median(draws_anp_1[,i+1])
  summaries$sd[i]     <- sd(draws_anp_1[,i+1])
}

# organise dem_class -- report age category based on age at start of period
plot_data_anp_1 <- left_join(x = summaries, y = counts_df_period1, by = 'dyad')
head(plot_data_anp_1)
plot_data_anp_1$age_cat_1 <- ifelse(plot_data_anp_1$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_1$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_1$age_start_1 > 19, 'A','P')))
plot_data_anp_1$age_cat_2 <- ifelse(plot_data_anp_1$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_1$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_1$age_start_2 > 19, 'A','P')))
plot_data_anp_1$age_catnum_1 <- ifelse(plot_data_anp_1$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_1$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_1$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_1$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_1$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_1$age_start_1 < 40, 6, 7))))))
plot_data_anp_1$age_catnum_2 <- ifelse(plot_data_anp_1$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_1$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_1$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_1$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_1$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_1$age_start_2 < 40, 6, 7))))))

plot_data_anp_1$age_dyad <- ifelse(plot_data_anp_1$age_catnum_1 >= plot_data_anp_1$age_catnum_2,
                                    paste(plot_data_anp_1$age_cat_1, plot_data_anp_1$age_cat_2, sep = ''),
                                    paste(plot_data_anp_1$age_cat_2, plot_data_anp_1$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_1, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_1, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries) # generally about the same as MOTNP, but higher SD
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_anp_1[plot_data_anp_1$age_dyad == 'AA' | plot_data_anp_1$age_dyad == 'AP' | plot_data_anp_1$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 1 ################
head(summaries)
length(unique(plot_data_anp_1$id_1))+1 # number of individuals = 55

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period1$id_1))+1,
                         NROW(unique(counts_df_period1$id_2))+1,
                         NROW(draws_anp_1)),
                    dimnames = list(c(unique(counts_df_period1$id_1),
                                      unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))]),
                                    c(unique(counts_df_period1$id_1),
                                      unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))]),
                                    NULL))
N <- nrow(counts_df_period1)

for (i in 1:N) {
  dyad_row <- counts_df_period1[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_1[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids1 <- c(unique(counts_df_period1$id_1), unique(counts_df_period1$id_2)[length(unique(counts_df_period1$id_2))])
males1 <- males[,c(21,6,9,22:53)]
males1 <- males1 %>% dplyr::filter(id %in% ids1)
males1

# create variables for different degrees of node connectedness
males1$degree_0.1 <- NA
males1$degree_0.2 <- NA
males1$degree_0.3 <- NA
males1$degree_0.4 <- NA
males1$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males1)){
  rows <- summaries[summaries$id_1 == males1$id[i] | summaries$id_2 == males1$id[i],]
  males1$degree_0.1[i] <- length(which(rows$median > 0.1))
  males1$degree_0.2[i] <- length(which(rows$median > 0.2))
  males1$degree_0.3[i] <- length(which(rows$median > 0.3))
  males1$degree_0.4[i] <- length(which(rows$median > 0.4))
  males1$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males1$degree_0.1 < males1$degree_0.2)
which(males1$degree_0.2 < males1$degree_0.3)
which(males1$degree_0.3 < males1$degree_0.4)
which(males1$degree_0.4 < males1$degree_0.5)

# age variable
males1$age <- lubridate::year(periods$period_start[periods$period == 1]) - males1$byr
summary(males1$age)
males1$age_class <- ifelse(males1$age < 10, 2,
                            ifelse(males1$age < 15, 3,
                                   ifelse(males1$age < 20, 4,
                                          ifelse(males1$age < 25, 5,
                                                 ifelse(males1$age < 40, 6, 7)))))

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = NA,
     vertex.size = 5,
     layout = coords)
plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     edge.color = 'black',
     vertex.size = 8,
     vertex.label = males1$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males1$age_class == 7,'seagreen4',
                          ifelse(males1$age_class == 6,'seagreen3',
                                 ifelse(males1$age_class == 5,'seagreen2',
                                        ifelse(males1$age_class == 4,'steelblue3',
                                               ifelse(males1$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.2,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.2,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males1$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males1$age_class == 7,'seagreen4',
                          ifelse(males1$age_class == 6,'seagreen3',
                                 ifelse(males1$age_class == 5,'seagreen2',
                                        ifelse(males1$age_class == 4,'steelblue3',
                                               ifelse(males1$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 1 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males1$id[which(males1$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males1$id[which(males1$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)
plot(g_mid_0.2,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_rng_0.2)$weight*3, 0),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.2)
plot(g_mid_0.2,
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_mid_0.2)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 7,'seagreen4',
                          ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 6,'seagreen3',
                                 ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 5,'seagreen2',
                                        ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 4,'steelblue3',
                                               ifelse(males1[which(males1$degree_0.2 > 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

par(mai = c(0.1,0.4,0.1,0.2))
plot(g_mid_0.2,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_rng_0.2)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.2)
plot(g_mid_0.2,
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_mid_0.2)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males1$p1[which(males1$degree_0.2 != 0)],
     vertex.label = males1$id[which(males1$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males1[which(males1$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

### save summary data -- period 1 ####
dyad_period_weights <- counts_df_non0
summaries$period <- 1
summaries <- summaries[,c(1,4:9)]
dyad_period_weights <- left_join(x = dyad_period_weights, y = summaries,
                                 by = c('dyad','period'))
rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.2,
   counts_df_period1, draws_anp_1, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males1, plot_data_anp_1, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids1, N, plot_cols)

################ 7.15) Run model on real standardised data -- period 15 ################
### create data list
counts_df_period15 <- counts_df_non0[counts_df_non0$period == 15,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period15),          # total number of times one or other of the dyad was observed
  together = counts_df_period15$event_count,    # count number of sightings seen together
  apart    = counts_df_period15$apart,          # count number of sightings seen apart
  period   = counts_df_period15$period)         # which period it's within

### Fit model
weight_anp_2.2_period15 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period15
#variable        mean     median    sd   mad         q5        q95 rhat ess_bulk ess_tail
#lp__      -132933.05 -132930.00 97.62 97.85 -133098.00 -132776.00 1.00     1394     1933
#weight[1]       0.05       0.05  0.04  0.03       0.01       0.12 1.00     3484     1794
#weight[2]       0.05       0.04  0.04  0.03       0.01       0.12 1.00     3514     1903
#weight[3]       0.08       0.07  0.05  0.05       0.01       0.17 1.00     4739     2340
#weight[4]       0.08       0.07  0.05  0.05       0.01       0.18 1.00     4905     2312
#weight[5]       0.18       0.17  0.09  0.09       0.05       0.35 1.00     4790     2435
#weight[6]       0.07       0.06  0.05  0.04       0.01       0.16 1.00     4968     2548
#weight[7]       0.20       0.19  0.09  0.09       0.08       0.36 1.00     5330     2555
#weight[8]       0.09       0.08  0.06  0.05       0.02       0.19 1.00     4911     2139
#weight[9]       0.06       0.05  0.04  0.04       0.01       0.14 1.00     4463     2009
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period15$output_files()[1])
draws1_anp_15 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period15$output_files()[2])
draws2_anp_15 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period15$output_files()[3])
draws3_anp_15 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period15$output_files()[4])
draws4_anp_15 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_15 <- rbind(draws1_anp_15, draws2_anp_15, draws3_anp_15, draws4_anp_15)

colnames(draws_anp_15)[2:ncol(draws_anp_15)] <- counts_df_period15$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_15), size = 30, replace = F)

### save data 
write_csv(draws_anp_15, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period15_22.04.15.csv')
rm(draws1_anp_15, draws2_anp_15, draws3_anp_15, draws4_anp_15)

### build traceplots -- period 15 ####
# reset plot window
dev.off()

# traceplot
plot(draws_anp_15[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_15[,plot_cols[2]], col = 'tan')
lines(draws_anp_15[,plot_cols[3]], col = 'orange')
lines(draws_anp_15[,plot_cols[4]], col = 'green')
lines(draws_anp_15[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_15[,plot_cols[6]], col = 'blue')
lines(draws_anp_15[,plot_cols[7]], col = 'red')
lines(draws_anp_15[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_15[,plot_cols[9]], col = 'purple')
lines(draws_anp_15[,plot_cols[10]],col = 'magenta')
lines(draws_anp_15[,plot_cols[11]],col = 'black')
lines(draws_anp_15[,plot_cols[12]], col = 'tan')
lines(draws_anp_15[,plot_cols[13]], col = 'orange')
lines(draws_anp_15[,plot_cols[14]], col = 'green')
lines(draws_anp_15[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_15[,plot_cols[16]], col = 'blue')
lines(draws_anp_15[,plot_cols[17]], col = 'red')
lines(draws_anp_15[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_15[,plot_cols[19]], col = 'purple')
lines(draws_anp_15[,plot_cols[20]],col = 'magenta')
lines(draws_anp_15[,plot_cols[21]],col = 'black')
lines(draws_anp_15[,plot_cols[22]], col = 'tan')
lines(draws_anp_15[,plot_cols[23]], col = 'orange')
lines(draws_anp_15[,plot_cols[24]], col = 'green')
lines(draws_anp_15[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_15[,plot_cols[26]], col = 'blue')
lines(draws_anp_15[,plot_cols[27]], col = 'red')
lines(draws_anp_15[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_15[,plot_cols[29]], col = 'purple')
lines(draws_anp_15[,plot_cols[30]],col = 'magenta')

### density plots -- period 15 ####
dens(draws_anp_15[,2], ylim = c(0,30),xlim = c(0,1), las = 1)
for(i in 1:30){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws_anp_15[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 15 ####
# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_15[,2:ncol(draws_anp_15)]),
                        min = rep(NA, ncol(draws_anp_15)-1),
                        max = rep(NA, ncol(draws_anp_15)-1),
                        mean = rep(NA, ncol(draws_anp_15)-1),
                        median = rep(NA, ncol(draws_anp_15)-1),
                        sd = rep(NA, ncol(draws_anp_15)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_15[,i+1])
  summaries$max[i]    <- max(draws_anp_15[,i+1])
  summaries$mean[i]   <- mean(draws_anp_15[,i+1])
  summaries$median[i] <- median(draws_anp_15[,i+1])
  summaries$sd[i]     <- sd(draws_anp_15[,i+1])
}

# organise dem_class -- report age category based on age at start of period
plot_data_anp_15 <- left_join(x = summaries, y = counts_df_period15, by = 'dyad')
head(plot_data_anp_15)
plot_data_anp_15$age_cat_1 <- ifelse(plot_data_anp_15$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_15$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_15$age_start_1 > 19, 'A','P')))
plot_data_anp_15$age_cat_2 <- ifelse(plot_data_anp_15$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_15$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_15$age_start_2 > 19, 'A','P')))
plot_data_anp_15$age_catnum_1 <- ifelse(plot_data_anp_15$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_15$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_15$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_15$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_15$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_15$age_start_1 < 40, 6, 7))))))
plot_data_anp_15$age_catnum_2 <- ifelse(plot_data_anp_15$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_15$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_15$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_15$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_15$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_15$age_start_2 < 40, 6, 7))))))

plot_data_anp_15$age_dyad <- ifelse(plot_data_anp_15$age_catnum_1 >= plot_data_anp_15$age_catnum_2,
                                    paste(plot_data_anp_15$age_cat_1, plot_data_anp_15$age_cat_2, sep = ''),
                                    paste(plot_data_anp_15$age_cat_2, plot_data_anp_15$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_15, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_15, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries) # generally about the same as MOTNP
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_anp_15[plot_data_anp_15$age_dyad == 'AA' | plot_data_anp_15$age_dyad == 'AP' | plot_data_anp_15$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))
hist(m_sum$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 15 ################
head(summaries)
length(unique(plot_data_anp_15$id_1))+1 # number of individuals = 181

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period15$id_1))+1,
                         NROW(unique(counts_df_period15$id_2))+1,
                         NROW(draws_anp_15)),
                    dimnames = list(c(unique(counts_df_period15$id_1),
                                      unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))]),
                                    c(unique(counts_df_period15$id_1),
                                      unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))]),
                                    NULL))
N <- nrow(counts_df_period15)

for (i in 1:N) {
  dyad_row <- counts_df_period15[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_15[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids15 <- c(unique(counts_df_period15$id_1), unique(counts_df_period15$id_2)[length(unique(counts_df_period15$id_2))])
males15 <- males[,c(21,6,9,22:53)]
males15 <- males15 %>% dplyr::filter(id %in% ids15)
males15

# create variables for different degrees of node connectedness
males15$degree_0.1 <- NA
males15$degree_0.2 <- NA
males15$degree_0.3 <- NA
males15$degree_0.4 <- NA
males15$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males15)){
  rows <- summaries[summaries$id_1 == males15$id[i] | summaries$id_2 == males15$id[i],]
  males15$degree_0.1[i] <- length(which(rows$median > 0.1))
  males15$degree_0.2[i] <- length(which(rows$median > 0.2))
  males15$degree_0.3[i] <- length(which(rows$median > 0.3))
  males15$degree_0.4[i] <- length(which(rows$median > 0.4))
  males15$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males15$degree_0.1 < males15$degree_0.2)
which(males15$degree_0.2 < males15$degree_0.3)
which(males15$degree_0.3 < males15$degree_0.4)
which(males15$degree_0.4 < males15$degree_0.5)

# age variable
males15$age <- lubridate::year(periods$period_start[periods$period == 15]) - males15$byr
summary(males15$age)
males15$age_class <- ifelse(males15$age < 10, 2,
                            ifelse(males15$age < 15, 3,
                                   ifelse(males15$age < 20, 4,
                                          ifelse(males15$age < 25, 5,
                                                 ifelse(males15$age < 40, 6, 7)))))

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = NA,
     vertex.size = 5,
     layout = coords)
plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     edge.color = 'black',
     vertex.size = 8,
     vertex.label = males15$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males15$age_class == 7,'seagreen4',
                          ifelse(males15$age_class == 6,'seagreen3',
                                 ifelse(males15$age_class == 5,'seagreen2',
                                        ifelse(males15$age_class == 4,'steelblue3',
                                               ifelse(males15$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.15,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.15,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males15$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males15$age_class == 7,'seagreen4',
                          ifelse(males15$age_class == 6,'seagreen3',
                                 ifelse(males15$age_class == 5,'seagreen2',
                                        ifelse(males15$age_class == 4,'steelblue3',
                                               ifelse(males15$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 15 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males15$id[which(males15$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males15$id[which(males15$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)
plot(g_mid_0.2,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_rng_0.2)$weight*3, 0),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.2)
plot(g_mid_0.2,
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_mid_0.2)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

par(mai = c(0.1,0.4,0.1,0.2))
plot(g_mid_0.2,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_rng_0.2)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.2)
plot(g_mid_0.2,
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_mid_0.2)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males15$p15[which(males15$degree_0.2 != 0)],
     vertex.label = males15$id[which(males15$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males15[which(males15$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

### save summary data -- period 15 ####
summaries$period <- 15
summaries <- summaries[,c(1,4:9)]
dyad_period_weights <- left_join(x = dyad_period_weights, y = summaries,
                                 by = c('dyad','period'))
colnames(dyad_period_weights)[36:40] <- c('min','max','mean','median','sd')
for(i in 1:nrow(dyad_period_weights)){
  dyad_period_weights$min[i] <- ifelse(is.na(dyad_period_weights$min[i]) == TRUE,
                                       ifelse(is.na(dyad_period_weights$min.x[i]) == FALSE,
                                              dyad_period_weights$min.x[i], NA),
                                       dyad_period_weights$min[i])
  dyad_period_weights$max[i] <- ifelse(is.na(dyad_period_weights$max[i]) == TRUE,
                                       ifelse(is.na(dyad_period_weights$max.x[i]) == FALSE,
                                              dyad_period_weights$max.x[i], NA),
                                       dyad_period_weights$max[i])
  dyad_period_weights$mean[i] <- ifelse(is.na(dyad_period_weights$mean[i]) == TRUE,
                                        ifelse(is.na(dyad_period_weights$mean.x[i]) == FALSE,
                                               dyad_period_weights$mean.x[i], NA),
                                        dyad_period_weights$mean[i])
  dyad_period_weights$median[i] <- ifelse(is.na(dyad_period_weights$median[i]) == TRUE,
                                          ifelse(is.na(dyad_period_weights$median.x[i]) == FALSE,
                                                 dyad_period_weights$median.x[i], NA),
                                          dyad_period_weights$median[i])
  dyad_period_weights$sd[i] <- ifelse(is.na(dyad_period_weights$sd[i]) == TRUE,
                                      ifelse(is.na(dyad_period_weights$sd.x[i]) == FALSE,
                                             dyad_period_weights$sd.x[i], NA),
                                      dyad_period_weights$sd[i])
}
dyad_period_weights <- dyad_period_weights[,c(1:30,36:40)]

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.2,
   counts_df_period15, draws_anp_15, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males15, plot_data_anp_15, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids15, N, plot_cols)

################ 7.28) Run model on real standardised data -- period 28 ################
### create data list
counts_df_period28 <- counts_df_non0[counts_df_non0$period == 28,]
counts_ls <- list(
  n_dyads  = nrow(counts_df_period28),          # total number of times one or other of the dyad was observed
  together = counts_df_period28$event_count,    # count number of sightings seen together
  apart    = counts_df_period28$apart,          # count number of sightings seen apart
  period   = counts_df_period28$period)         # which period it's within

### Fit model
weight_anp_2.2_period28 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_anp_2.2_period28
# variable       mean     median     sd    mad         q5        q95 rhat ess_bulk ess_tail
#lp__      -124881.32 -124878.00 109.57 111.19 -125063.00 -124701.00 1.00     1361     2334
#weight[1]       0.14       0.13   0.09   0.09       0.03       0.33 1.00     3988     2130
#weight[2]       0.12       0.11   0.07   0.07       0.02       0.26 1.00     3463     1692
#weight[3]       0.11       0.10   0.06   0.06       0.03       0.22 1.00     4186     2399
#weight[4]       0.10       0.09   0.07   0.06       0.02       0.23 1.00     3647     2072
#weight[5]       0.15       0.13   0.09   0.09       0.03       0.32 1.00     3699     2219
#weight[6]       0.18       0.16   0.09   0.09       0.05       0.34 1.00     4891     2856
#weight[7]       0.22       0.20   0.11   0.11       0.06       0.42 1.00     3636     2532
#weight[8]       0.20       0.19   0.10   0.10       0.06       0.39 1.00     4324     2476
#weight[9]       0.09       0.08   0.06   0.06       0.02       0.21 1.00     3740     1912
rm(counts_ls)

# showing 10 of 402700 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
output1 <- read_cmdstan_csv(weight_anp_2.2_period28$output_files()[1])
draws1_anp_28 <- as.data.frame(output1$post_warmup_draws)
rm(output1)

output2 <- read_cmdstan_csv(weight_anp_2.2_period28$output_files()[2])
draws2_anp_28 <- as.data.frame(output2$post_warmup_draws)
rm(output2)

output3 <- read_cmdstan_csv(weight_anp_2.2_period28$output_files()[3])
draws3_anp_28 <- as.data.frame(output3$post_warmup_draws)
rm(output3)

output4 <- read_cmdstan_csv(weight_anp_2.2_period28$output_files()[4])
draws4_anp_28 <- as.data.frame(output4$post_warmup_draws)
rm(output4)

draws_anp_28 <- rbind(draws1_anp_28, draws2_anp_28, draws3_anp_28, draws4_anp_28)

colnames(draws_anp_28)[2:ncol(draws_anp_28)] <- counts_df_period28$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:ncol(draws_anp_28), size = 30, replace = F)

### save data 
write_csv(draws_anp_28, 'data_processed/anp_bayesian_edgedistributions_a2.b2_period28_22.04.15.csv')
rm(draws1_anp_28, draws2_anp_28, draws3_anp_28, draws4_anp_28)

### build traceplots  -- period 28 ####
# reset plot window
dev.off()

# traceplot
plot(draws_anp_28[,plot_cols[1]], type = 'l', ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_anp_28[,plot_cols[2]], col = 'tan')
lines(draws_anp_28[,plot_cols[3]], col = 'orange')
lines(draws_anp_28[,plot_cols[4]], col = 'green')
lines(draws_anp_28[,plot_cols[5]], col = 'chocolate')
lines(draws_anp_28[,plot_cols[6]], col = 'blue')
lines(draws_anp_28[,plot_cols[7]], col = 'red')
lines(draws_anp_28[,plot_cols[8]], col = 'seagreen')
lines(draws_anp_28[,plot_cols[9]], col = 'purple')
lines(draws_anp_28[,plot_cols[10]],col = 'magenta')
lines(draws_anp_28[,plot_cols[11]],col = 'black')
lines(draws_anp_28[,plot_cols[12]], col = 'tan')
lines(draws_anp_28[,plot_cols[13]], col = 'orange')
lines(draws_anp_28[,plot_cols[14]], col = 'green')
lines(draws_anp_28[,plot_cols[15]], col = 'chocolate')
lines(draws_anp_28[,plot_cols[16]], col = 'blue')
lines(draws_anp_28[,plot_cols[17]], col = 'red')
lines(draws_anp_28[,plot_cols[18]], col = 'seagreen')
lines(draws_anp_28[,plot_cols[19]], col = 'purple')
lines(draws_anp_28[,plot_cols[20]],col = 'magenta')
lines(draws_anp_28[,plot_cols[21]],col = 'black')
lines(draws_anp_28[,plot_cols[22]], col = 'tan')
lines(draws_anp_28[,plot_cols[23]], col = 'orange')
lines(draws_anp_28[,plot_cols[24]], col = 'green')
lines(draws_anp_28[,plot_cols[25]], col = 'chocolate')
lines(draws_anp_28[,plot_cols[26]], col = 'blue')
lines(draws_anp_28[,plot_cols[27]], col = 'red')
lines(draws_anp_28[,plot_cols[28]], col = 'seagreen')
lines(draws_anp_28[,plot_cols[29]], col = 'purple')
lines(draws_anp_28[,plot_cols[30]],col = 'magenta')

### density plots  -- period 28 ####
dens(draws_anp_28[,2], ylim = c(0,20),xlim = c(0,1), las = 1)
for(i in 1:30){  # looks like a poisson distribution - peaks just above 0 and then exponentiol decline. almost nothing above 0.5
  dens(add = T, draws_anp_28[,plot_cols[i]], col = col.alpha('blue', alpha = 0.3))
}

### summarise data -- period 28 ####
# summarise -- look for any anomalies in draw values
summaries <- data.frame(dyad = colnames(draws_anp_28[,2:ncol(draws_anp_28)]),
                        min = rep(NA, ncol(draws_anp_28)-1),
                        max = rep(NA, ncol(draws_anp_28)-1),
                        mean = rep(NA, ncol(draws_anp_28)-1),
                        median = rep(NA, ncol(draws_anp_28)-1),
                        sd = rep(NA, ncol(draws_anp_28)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_anp_28[,i+1])
  summaries$max[i]    <- max(draws_anp_28[,i+1])
  summaries$mean[i]   <- mean(draws_anp_28[,i+1])
  summaries$median[i] <- median(draws_anp_28[,i+1])
  summaries$sd[i]     <- sd(draws_anp_28[,i+1])
}

# organise dem_class -- report age category based on age at start of period
plot_data_anp_28 <- left_join(x = summaries, y = counts_df_period28, by = 'dyad')
head(plot_data_anp_28)
plot_data_anp_28$age_cat_1 <- ifelse(plot_data_anp_28$age_start_1 < 6, 'C',
                                     ifelse(plot_data_anp_28$age_start_1 < 11, 'J',
                                            ifelse(plot_data_anp_28$age_start_1 > 19, 'A','P')))
plot_data_anp_28$age_cat_2 <- ifelse(plot_data_anp_28$age_start_2 < 6, 'C',
                                     ifelse(plot_data_anp_28$age_start_2 < 11, 'J',
                                            ifelse(plot_data_anp_28$age_start_2 > 19, 'A','P')))
plot_data_anp_28$age_catnum_1 <- ifelse(plot_data_anp_28$age_start_1 < 6, 1,
                                        ifelse(plot_data_anp_28$age_start_1 < 11, 2,
                                               ifelse(plot_data_anp_28$age_start_1 < 16, 3,
                                                      ifelse(plot_data_anp_28$age_start_1 < 20, 4,
                                                             ifelse(plot_data_anp_28$age_start_1 < 25, 5,
                                                                    ifelse(plot_data_anp_28$age_start_1 < 40, 6, 7))))))
plot_data_anp_28$age_catnum_2 <- ifelse(plot_data_anp_28$age_start_2 < 6, 1,
                                        ifelse(plot_data_anp_28$age_start_2 < 11, 2,
                                               ifelse(plot_data_anp_28$age_start_2 < 16, 3,
                                                      ifelse(plot_data_anp_28$age_start_2 < 20, 4,
                                                             ifelse(plot_data_anp_28$age_start_2 < 25, 5,
                                                                    ifelse(plot_data_anp_28$age_start_2 < 40, 6, 7))))))

plot_data_anp_28$age_dyad <- ifelse(plot_data_anp_28$age_catnum_1 >= plot_data_anp_28$age_catnum_2,
                                    paste(plot_data_anp_28$age_cat_1, plot_data_anp_28$age_cat_2, sep = ''),
                                    paste(plot_data_anp_28$age_cat_2, plot_data_anp_28$age_cat_1, sep = ''))

# boxplot edge weights by demographic type of dyad -- all types
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_anp_28, aes(y = mean, x = age_dyad))+
    geom_boxplot(notch = T, outlier.shape = 4, fill = c('blue','purple','blue','grey','purple','blue'))+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown)',
         y = 'mean edge weight')+
    coord_flip())

# scatter plot of all adult edges by age difference
(edge_vs_agediff <- 
    ggplot(plot_data_anp_28, aes(x = age_diff, y = mean))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    scale_x_continuous('age difference between dyad members',
                       expand = c(0.02,0))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))

# values for reporting
summary(summaries) # generally about the same as MOTNP
quantile(summaries$median, seq(0,1,length.out = 101))
quantile(summaries$mean,   seq(0,1,length.out = 101))
hist(summaries$median, breaks = 100, xlim = c(0,1), las = 1,
     xlab = 'median estimated dyad strength', main = '', col = 'skyblue')

m_sum <- plot_data_anp_28[plot_data_anp_28$age_dyad == 'AA' | plot_data_anp_28$age_dyad == 'AP' | plot_data_anp_28$age_dyad == 'PP',] ; m_sum <- m_sum[!is.na(m_sum$dyad),]
quantile(m_sum$median, seq(0,1,length.out = 101))

# clean up
rm(edge_vs_agediff, edge_vs_demtype_all, m_sum)

### create network plots -- period 28 ################
head(summaries)
length(unique(plot_data_anp_28$id_1))+1 # number of individuals = 198

par(mai = c(0.1,0.1,0.1,0.1))

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df_period28$id_1))+1,
                         NROW(unique(counts_df_period28$id_2))+1,
                         NROW(draws_anp_28)),
                    dimnames = list(c(unique(counts_df_period28$id_1),
                                      unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))]),
                                    c(unique(counts_df_period28$id_1),
                                      unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))]),
                                    NULL))
N <- nrow(counts_df_period28)

for (i in 1:N) {
  dyad_row <- counts_df_period28[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_anp_28[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- (adj_upper - adj_lower) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ids28 <- c(unique(counts_df_period28$id_1), unique(counts_df_period28$id_2)[length(unique(counts_df_period28$id_2))])
males28 <- males[,c(21,6,9,22:53)]
males28 <- males28 %>% dplyr::filter(id %in% ids28)
males28

# create variables for different degrees of node connectedness
males28$degree_0.1 <- NA
males28$degree_0.2 <- NA
males28$degree_0.3 <- NA
males28$degree_0.4 <- NA
males28$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(males28)){
  rows <- summaries[summaries$id_1 == males28$id[i] | summaries$id_2 == males28$id[i],]
  males28$degree_0.1[i] <- length(which(rows$median > 0.1))
  males28$degree_0.2[i] <- length(which(rows$median > 0.2))
  males28$degree_0.3[i] <- length(which(rows$median > 0.3))
  males28$degree_0.4[i] <- length(which(rows$median > 0.4))
  males28$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(males28$degree_0.1 < males28$degree_0.2)
which(males28$degree_0.2 < males28$degree_0.3)
which(males28$degree_0.3 < males28$degree_0.4)
which(males28$degree_0.4 < males28$degree_0.5)

# age variable
males28$age <- lubridate::year(periods$period_start[periods$period == 28]) - males28$byr
summary(males28$age)
males28$age_class <- ifelse(males28$age < 10, 2,
                            ifelse(males28$age < 15, 3,
                                   ifelse(males28$age < 20, 4,
                                          ifelse(males28$age < 25, 5,
                                                 ifelse(males28$age < 40, 6, 7)))))

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = NA,
     vertex.size = 5,
     layout = coords)
plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     edge.color = 'black',
     vertex.size = 8,
     vertex.label = males28$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males28$age_class == 7,'seagreen4',
                          ifelse(males28$age_class == 6,'seagreen3',
                                 ifelse(males28$age_class == 5,'seagreen2',
                                        ifelse(males28$age_class == 4,'steelblue3',
                                               ifelse(males28$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.3,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 8,
     vertex.label = males28$id,
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(males28$age_class == 7,'seagreen4',
                          ifelse(males28$age_class == 6,'seagreen3',
                                 ifelse(males28$age_class == 5,'seagreen2',
                                        ifelse(males28$age_class == 4,'steelblue3',
                                               ifelse(males28$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords, add = TRUE)

### only elephants degree > 0.2 -- period 28 ####
g_mid_0.2 <- delete.vertices(graph = g_mid, v = males28$id[which(males28$degree_0.2 == 0)])
g_rng_0.2 <- delete.vertices(graph = g_rng, v = males28$id[which(males28$degree_0.2 == 0)])

set.seed(3)
coords_0.2 <- layout_nicely(g_mid_0.2)
plot(g_mid_0.2,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_rng_0.2)$weight*3, 0),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.2)
plot(g_mid_0.2,
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_mid_0.2)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

par(mai = c(0.1,0.4,0.1,0.2))
plot(g_mid_0.2,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_rng_0.2)$weight*3, 0),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_0.2)
plot(g_mid_0.2,
     edge.width = ifelse(E(g_mid_0.2)$weight > 0.2, E(g_mid_0.2)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males28$p28[which(males28$degree_0.2 != 0)],
     vertex.label = males28$id[which(males28$degree_0.2 != 0)],
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color= ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 7,'seagreen4',
                          ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 6,'seagreen3',
                                 ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 5,'seagreen2',
                                        ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 4,'steelblue3',
                                               ifelse(males28[which(males28$degree_0.2 != 0),]$age_class == 3,'steelblue1',
                                                      'yellow'))))),
     layout = coords_0.2, add = TRUE)

### save summary data -- period 28 ####
summaries$period <- 28
summaries <- summaries[,c(1,4:9)]
dyad_period_weights <- left_join(x = dyad_period_weights, y = summaries,
                                 by = c('dyad','period'))
colnames(dyad_period_weights)[36:40] <- c('min','max','mean','median','sd')
for(i in 1:nrow(dyad_period_weights)){
  dyad_period_weights$min[i] <- ifelse(is.na(dyad_period_weights$min[i]) == TRUE,
                                       ifelse(is.na(dyad_period_weights$min.x[i]) == FALSE,
                                              dyad_period_weights$min.x[i], NA),
                                       dyad_period_weights$min[i])
  dyad_period_weights$max[i] <- ifelse(is.na(dyad_period_weights$max[i]) == TRUE,
                                       ifelse(is.na(dyad_period_weights$max.x[i]) == FALSE,
                                              dyad_period_weights$max.x[i], NA),
                                       dyad_period_weights$max[i])
  dyad_period_weights$mean[i] <- ifelse(is.na(dyad_period_weights$mean[i]) == TRUE,
                                        ifelse(is.na(dyad_period_weights$mean.x[i]) == FALSE,
                                               dyad_period_weights$mean.x[i], NA),
                                        dyad_period_weights$mean[i])
  dyad_period_weights$median[i] <- ifelse(is.na(dyad_period_weights$median[i]) == TRUE,
                                          ifelse(is.na(dyad_period_weights$median.x[i]) == FALSE,
                                                 dyad_period_weights$median.x[i], NA),
                                          dyad_period_weights$median[i])
  dyad_period_weights$sd[i] <- ifelse(is.na(dyad_period_weights$sd[i]) == TRUE,
                                      ifelse(is.na(dyad_period_weights$sd.x[i]) == FALSE,
                                             dyad_period_weights$sd.x[i], NA),
                                      dyad_period_weights$sd[i])
}
dyad_period_weights <- dyad_period_weights[,c(1:30,36:40)]

rm(adj_lower, adj_mid, adj_range, adj_upper, coords, coords_0.2,
   counts_df_period28, draws_anp_28, dyad_row, g_mid,g_mid_0.2, g_rng, g_rng_0.2, males28, plot_data_anp_28, rows, summaries, adj_quantiles, adj_tensor, i, ids, ids28, N, plot_cols)
