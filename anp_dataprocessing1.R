# Amboseli data
### Information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022
### Set up ####
library(tidyverse)
library(rstan)
library(rethinking)
library(cmdstanr)
library(lubridate)

### Import sightings data ####
ate <- readxl::read_excel(path = 'data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx', sheet = 1)
ate <- janitor::clean_names(ate)
str(ate)
## obs_id = unique group number through study
summary(ate$obs_id)
length(unique(ate$obs_id))
## casename = MaleID number -- make character string obvious so clearly different from node_id later on
ate$id <- paste('M',ate$casename, sep = '')
sort(unique(ate$id))
ate$node_id <- as.integer(as.factor(ate$casename))
## musth = Musth status, true/false
table(ate$musth)
## obs_date = Date of observation -- convert to date value
summary(ate$obs_date)
#ate$obs_date[1]
#ate[ate$obs_date == '1972-09-02',]
#ate[ate$obs_date == '1972-09-02 UTC',]
#ate[ate$obs_date == '1972-09-02 00:00:00',]
#ate[ate$obs_date == 1972-09-02,]
ate$obs_date <- lubridate::as_date(ate$obs_date)
summary(ate$obs_date)
#ate[ate$obs_date == '1972-09-02',]
#ate[ate$obs_date == '1972-09-02 UTC',]
#ate[ate$obs_date == '1972-09-02 00:00:00',]

## obs_time = Time of observation -- convert to time value
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
#ate$correct_time <- lubridate::hour(ate$correct_time)*60*60 + lubridate::minute(ate$correct_time) + lubridate::second(ate$correct_time)
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
summary(ate$correct_time_hms)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)
summary(ate$corrected_time)
ate[ate$corrected_time == 0,]

## obs_num = Encounter number of the day -- standardise
test <- ate[,c('obs_id', 'obs_date', 'correct_time_hms', 'obs_num', 'grp_size')] %>% 
  distinct()
lu <- function(x) { length(unique(x)) }
ate_nums <- tapply(X = ate$obs_num, INDEX = ate$obs_date, FUN = lu )
test_nums <- tapply(X = test$obs_num, INDEX = test$obs_date, FUN = lu )
which(ate_nums != test_nums)
ate_nums[1:50] ; test_nums[1:50]

table(ate$obs_date[which(ate$obs_num == '0')])
table(ate$obs_date[which(ate$obs_num == '00')])
ate$obs_num <- ifelse(ate$obs_num == '0','00', ate$obs_num)
table(ate$obs_num)

table(ate$obs_date[which(ate$obs_num == '0a')])
table(ate$obs_date[which(ate$obs_num == '0A')])
ate$obs_num <- ifelse(ate$obs_num == '0a','0A', ate$obs_num)
table(ate$obs_num)

table(ate$obs_date[which(ate$obs_num == '0b')])
table(ate$obs_date[which(ate$obs_num == '0B')])
ate$obs_num <- ifelse(ate$obs_num == '0b','0B', ate$obs_num)
table(ate$obs_num)

table(test$obs_date[which(test$obs_num == '1')])
table(test$obs_date[which(test$obs_num == '01')])
unique(ate$obs_num[which(ate$obs_date == '2020-08-28')])
ate$obs_num <- ifelse(ate$obs_num == '1','01', ate$obs_num)
table(ate$obs_num)

int_fact <- function(x) { as.integer(as.factor(x)) }
test$obs_num_std <- tapply(X = test$obs_num, INDEX = test$obs_date, FUN = int_fact)

ate$obs_num_std <- NA
for(i in 1:length(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_std <- int_fact(sort(date_row$obs_num))
  ate$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}

## grid_code = 1km grid square, or larger location
table(ate$grid_code)

## utm_lat and utm_long = There is no mask applied to GPS at the moment, so this needs checking for outliers and impossible values, as well as note that short values are possible -- try plotting and see where they are
summary(ate$utm_lat)
summary(ate$utm_long)
no_gps <- unique(ate$obs_date[which(ate$utm_lat == 0)])
gps <- unique(ate$obs_date[which(ate$utm_lat != 0)])
min(gps)
no_gps[4100:5023] # many that are after the start of using actual GPS

which(ate$utm_lat == 0 & ate$utm_long != 0) # no latitude of 0 for high longitude 
which(ate$utm_lat != 0 & ate$utm_long == 0) # 5 have latitude value but no longitude

plot(utm_long ~ utm_lat, data = ate[ate$utm_lat != 0 & ate$utm_long != 0,])

summary(ate$utm_long[which(ate$utm_lat != 0 & ate$utm_long != 0)]) # should all be around 37
summary(ate$utm_lat[which(ate$utm_lat != 0 & ate$utm_long != 0)])  # should all be around -2

plot(utm_long ~ utm_lat, data = ate[ate$utm_lat != 0 & ate$utm_long != 0,],
     xlim = c(1,3), ylim = c(30,40))         # nothing in the range it should actually be

plot(utm_long ~ utm_lat, data = ate[ate$utm_lat != 0 & ate$utm_long != 0,],
     xlim = c(0,100000), ylim = c(0,400000)) # large mass of points but doesn't seem to be even an order of magnitude wrong to explain it

## hab_code_1 = habitat code for most of the group: 0 = Not described, 01 = short grass plain, 02 = Consimilis tall grasslands, 03 = Acacia tortilis woodland, 04 =	Acacia xanthophloea woodland, 05 = Salvadora/sueda, 06 = Palm woodland, 07 = Swamp edge woodland, 08 = Swamp edge, 09 =	Swamp, A =	Bushed grassland, B =	Open bushland north, C = Open bushland south, D	= Dense bushland north, E =	Acacia nubica, F = Dense bushland south, U = Unknown
## hab_code_2 = second habitat code if >1 category
sort(unique(ate$hab_code_1))
sort(unique(ate$hab_code_2))

## bull_q_r = quality of recognition for males (3 = knew all, 2 = knew > half, 1 = knew < half, U = unknown, blanks relate to data edits I haven't completed in the base data)
table(ate$bull_q_r)

## obs_type = Group Type: B = male only, M = males and females, U = unknown
table(ate$obs_type)

## grp_q_c = Group Quality of count: 3 = exact, 2 = good estimate, 1 = estimate, 0 = no estimate
table(ate$grp_q_c)

## grp_size = N elephants encountered in group; -1 is missing data. Zeros are typos that I haven't had chance to correct at present
table(ate$grp_size)
table(test$grp_size)

## bull_q_c = Male Quality of count: 3 = exact, 2 = good estimate, 1 = estimate, 0 = no estimate
table(ate$bull_q_c)

## num_bulls = N independent males present. Zeros are typos I haven't had chance to correct at present
table(ate$num_bulls)

## bulls_1_2 = Males aged 10-24, present/absent (1 = present,0 = absent, U = unknown)
table(ate$bulls_1_2)

## bulls_3_5 = Males aged 25+, present/absent (1 = present,0 = absent, U = unknown)
table(ate$bulls_3_5)

## musth_male = Musth male present/ absent (1 = present,0 = absent, U = unknown)
table(ate$musth_male)

## oestrus_fem = Oestrue female present/absent ((1 = present,0 = absent, U = unknown)
table(ate$oestrus_fem)

## act_code = Group activiity: 00 = not specified, 01 = walking while feeding, 02 = feeding, 03 = resting, 04 =	comfort (dusting, mudwallowing etc), 05 =	interacting, 06 =	drinking, 07 = walking, 08 = standing vigilant, 09 = more than one activity
table(ate$act_code)

## rearrange
ate <- ate[,c(1:3,25,26,4,27,28,8,29,9:24)]
colnames(ate)

## clean environment
rm(ate_nums, date_row, test, gps, i, no_gps, test_nums)

### Import nodes data ####
nodes <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx')
nodes <- janitor::clean_names(nodes)

## make character string of ID number
nodes$id <- paste('M', nodes$casename, sep = '')

## check ID numbers match in both data frames
nodes_id <- sort(unique(nodes$id))
ate_id <- sort(unique(ate$id))
length(nodes_id) ; length(ate_id) ## DO NOT MATCH, OR EVEN CLOSE!

### create group-by-individual matrix ####
ate$obs_id_std <- int_fact(ate$obs_id)
ate_asnipe <- ate[,c(4,27)]
colnames(ate_asnipe) <- c('ID','group')
ate_asnipe <- data.table::setDT(ate_asnipe) # just converts to a data table, no other change. 
gbi_matrix <- spatsoc::get_gbi(DT = ate_asnipe, group = 'group', id = 'ID')
# NOTE: THE ORDER OF SIGHTINGS IS ACCORDING TO ate$obs_id_std NOT ate$obs_id -- USE OBS_ID_STD TO MATCH UP SIGHTING INFORMATION TO INDIVIDUALS

gbi_matrix_test <- gbi_matrix[1:200,1:20]
#rownames(gbi_matrix_test) <- ate_asnipe$group

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), row_num = numeric())
for (obs_id in 1:nrow(gbi_matrix_test)) {
  for (i in which(gbi_matrix_test[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix_test)) {
      if (i != j) {
        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix_test[obs_id, i] == gbi_matrix_test[obs_id, j]),
                                           row_num = obs_id)
      }
    }
  }
}
gbi_df
table(gbi_df$row_num)
gbi_df$obs_id <- NA
for(i in 1:length(gbi_df$obs_id)){
  gbi_df$obs_id[i] <- unique(ate$obs_id[which(ate$obs_id_std == gbi_df$row_num[i])])
}

### code to convert gbi matrix format to dyadic data frame, shared by Prof Dan Franks and available from @JHart96 GitHub repository (https://github.com/JHart96/bison_examples/blob/main/examples/convert_gbi.md) -- NOTE: this step takes a long time to run
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1:nrow(gbi_matrix)) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_rows1.2219089.csv', delim = ',')

nrow(gbi_matrix) # 24174 observations -- do in chunks of 100 sightings (900 sightings takes ~2 days)
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1:100) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10) {print(obs_id)}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings1.100.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 101:300) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings101.300.csv', delim = ',')

compare_1 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings1.100.csv')
compare_2 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings101.300.csv')
compare_3 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_rows1.2219089.csv')

compare_2.1 <- compare_2[1:13167,]
compare_3.1 <- compare_3[1:13167,]
which(compare_3$node_1 == 1 & compare_3$node_2 == 5 & compare_3$obs_id == 101) # 141373
141373+514206 # 655579
compare_3.2 <- compare_3[141373:655578,]

which(compare_1$node_1 != compare_3.1$node_1) # running altogether vs in sections did not alter which ones it did
which(compare_1$node_2 != compare_3.1$node_2) # running altogether vs in sections did not alter which ones it did
which(compare_1$social_event != compare_3.1$social_event) # running altogether vs in sections did not alter which ones it did
which(compare_1$obs_id != compare_3.1$obs_id) # running altogether vs in sections did not alter which ones it did

which(compare_1$node_1 != compare_2.1$node_1) # definitely started over on a different part
which(compare_1$node_2 != compare_2.1$node_2) # definitely started over on a different part
which(compare_1$social_event != compare_2.1$social_event) # definitely started over on a different part
which(compare_1$obs_id != compare_2.1$obs_id) # definitely started over on a different part

which(compare_2$node_1 != compare_3.2$node_1) # running altogether vs in sections did not alter which ones it did
which(compare_2$node_2 != compare_3.2$node_2) # running altogether vs in sections did not alter which ones it did
which(compare_2$social_event != compare_3.2$social_event) # running altogether vs in sections did not alter which ones it did
which(compare_2$obs_id != compare_3.2$obs_id) # running altogether vs in sections did not alter which ones it did

### same script as run previously on MOTNP data, but broken down into blocks of 250 sightings else it will take about a year to run!
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1:250) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings1.250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 251:500) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings251.500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 501:750) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings501.750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 751:999) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings751.999.csv', delim = ',')



gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1000:2000) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      print(new) # print in nice format
      if (i != j) {
        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
write_csv(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings1000.2000.csv')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2001:2250) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings2001.2250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2251:2500) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings2251.2500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2501:2750) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings2501.2750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2751:2999) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings2751.2999.csv', delim = ',')

test1 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings1.250.csv')
test2 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings251.500.csv')
test3 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings501.750.csv')
test4 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings751.999.csv')
test5 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings1000.2000.csv')
test6 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings2001.2250.csv')
test7 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings2501.2750.csv')
test8 <- read_csv('data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings2751.2999.csv')

length(test1$node_1)
length(test2$node_1)
length(test3$node_1)
length(test4$node_1)
length(test5$node_1)
length(test6$node_1)
length(test7$node_1)
length(test8$node_1)

table(test1$obs_id)
table(test2$obs_id)
table(test3$obs_id)
table(test4$obs_id)
table(test5$obs_id)
table(test6$obs_id)
table(test7$obs_id)
table(test8$obs_id)

length(unique(test1$node_1)) ; length(unique(test1$node_2))
summary(test1$node_1)
cumsum(1:694)[693] # should be 240471 dyads per observation. test1$obs_id == 1 only compares elephant 408 to everyone else, never tries 1--2, 1--3, 2--3...  test1$obs_id == 2 only compares elephants 39&48 to everyone else.

gbi_matrix[,408]    # elephant 408 observed 148 times
sum(gbi_matrix[1,]) # only one elephant observed
gbi_matrix[1,408]   # only considers the elephant that was seen, not every dyad

test <- rbind(test1, test2, test3, test4, test5, test6, test7, test8)
test_distinct <- distinct(test)
nrow(test)
nrow(test_distinct)
which(test$node_1 != test_distinct$node_1)[1] # 1425
test[1424:1426,]
test_distinct[1424:1426,]
which(test$node_1 == 39 & test$node_2 == 48 & test$obs_id == 2) # 740 and 1425 -- 39 and 48 are the two elephants recorded as present so there is a duplicate for node_1 and node_2 
nrow(test) - nrow(test_distinct)
test1_distinct <- distinct(test1)
nrow(test1) - nrow(test1_distinct) # 1979 rows removed as duplicates
sum(gbi_matrix[1:250,])            # 727 elephants observed in first 250 sightings

sightings <- data.frame(encounter = sort(unique(ate$obs_id_std)),
                        total_id = NA,
                        mat_sum = NA)
for( i in 1:nrow(sightings)){
  obs_row <- ate[ate$obs_id_std == sightings$encounter[i],]
  sightings$total_id[i] <- nrow(obs_row)
  sightings$mat_sum[i] <- sum(gbi_matrix[i,])
}
summary(sightings$total_id)
summary(sightings$mat_sum)

sightings$encounter_id <- as.numeric(sightings$encounter)
sightings <- sightings[sightings$encounter_id < 251,]
sum(sightings$total_id) # 729 elephants in first 150 sightings??
sum(sightings$mat_sum)  # 727 elephants in first 150 sightings??
which(sightings$total_id != sightings$mat_sum) # 195
sightings[195,] # 4 identified, 2 in matrix

which(gbi_matrix[195,] == 1) # M83 and M96
ate[ate$obs_id_std == 195,]  # duplicated rows for some reason?

ate_distinct <- distinct(ate) # 55146 down to 55141
ate_distinct[ate_distinct$obs_id_std == 195,] # removed this weirdness!

sightings <- sightings[sightings$mat_sum > 1,] # 131 sightings where there will be duplication in gbi_df and therefore require removal -- should in total come to 1979 duplicate dyads
sum(sightings$mat_sum) # 608
sightings$dyads <- NA
for(i in 1:nrow(sightings)){
  j <- sightings$mat_sum[i]
  sightings$dyads[i] <- cumsum(1:j)[j-1]
}
sum(sightings$dyads) # 1979!!! YES!!
# all seems to be working as planned after all!


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3001:3250) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings3001.3250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3251:3500) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings3251.3500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3501:3750) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings3501.3750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3751:3999) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings3751.3999.csv', delim = ',')



gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3001:3250) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings3001.3250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3251:3500) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings3251.3500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3501:3750) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings3501.3750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3751:3999) {
  for (i in which(gbi_matrix[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi_matrix)) {
      if (i != j) {
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi_df[nrow(gbi_df) + 1, ] <- list(node_1 = node_1,
                                           node_2 = node_2,
                                           social_event = (gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j]),
                                           obs_id = obs_id)
      }
    }
  }
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.03.03_sightings3751.3999.csv', delim = ',')



