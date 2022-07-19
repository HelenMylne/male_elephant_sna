# Amboseli data
#### Information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022
#### Set up ####
library(tidyverse)
library(rstan)
library(rethinking)
library(cmdstanr)
library(lubridate)

#### Import sightings data ####
ate <- read_csv('data_processed/anp_sightings_updated_22.06.22.csv') %>% janitor::clean_names()
str(ate)

## casename = MaleID number -- make character string obvious so clearly different from node_id later on
ate$id <- paste('M',ate$casename, sep = '')
sort(unique(ate$id))
ate$node_id <- as.integer(as.factor(ate$casename))

## obs_date = Date of observation -- convert to date value
ate$obs_date <- lubridate::as_date(ate$obs_date)

## obs_time = Time of observation -- convert to time value
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
#ate$correct_time <- lubridate::hour(ate$correct_time)*60*60 + lubridate::minute(ate$correct_time) + lubridate::second(ate$correct_time)
ate$correct_time_hms <- hms::as_hms(ate$correct_time)
ate$corrected_time <- lubridate::hour(ate$correct_time_hms)*60*60 + lubridate::minute(ate$correct_time_hms) + lubridate::second(ate$correct_time_hms)
summary(ate$corrected_time)

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

ate$obs_num_std <- NA
for(i in 1:length(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_std <- as.integer(as.factor((sort(date_row$obs_num))))
  ate$obs_num_std[i] <- date_row$obs_num_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}

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

#### Import nodes data ####
nodes <- readxl::read_excel('data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()

## make character string of ID number
nodes$id <- paste('M', nodes$casename, sep = '')

## check ID numbers match in both data frames
nodes_id <- sort(unique(nodes$id))
ate_id <- sort(unique(ate$id))
length(nodes_id) ; length(ate_id) ## DO NOT MATCH, OR EVEN CLOSE!

#### create group-by-individual matrix ####
ate$obs_id_std <- as.integer(as.factor(ate$obs_id))
ate_asnipe <- ate[,c(4,27)]
colnames(ate_asnipe) <- c('ID','group')
ate_asnipe <- data.table::setDT(ate_asnipe) # just converts to a data table, no other change. 
gbi_matrix <- spatsoc::get_gbi(DT = ate_asnipe, group = 'group', id = 'ID')
# NOTE: THE ORDER OF SIGHTINGS IS ACCORDING TO ate$obs_id_std NOT ate$obs_id -- USE OBS_ID_STD TO MATCH UP SIGHTING INFORMATION TO INDIVIDUALS

#### convert group-by-individual matrix to dyadic data frame of sightings ####
### code to convert gbi matrix format to dyadic data frame, shared by Prof Dan Franks and available from @JHart96 GitHub repository (https://github.com/JHart96/bison_examples/blob/main/examples/convert_gbi.md) -- NOTE: this step takes at least 2 weeks to run
## sightings 1-2000 ####
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings1.250.csv', delim = ',')

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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings251.500.csv', delim = ',')

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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings501.750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 751:1000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings751.1000.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1001:1250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings1001.1250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1251:1500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings1251.1500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1501:1750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings1501.1750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1751:2000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings1751.2000.csv', delim = ',')

## sightings 2001-4000 ####
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings2001.2250.csv', delim = ',')

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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings2251.2500.csv', delim = ',')

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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings2501.2750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2751:3000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings2751.3000.csv', delim = ',')


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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings3001.3250.csv', delim = ',')

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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings3251.3500.csv', delim = ',')

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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings3501.3750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3751:4000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings3751.4000.csv', delim = ',')

## sightings 4001-6000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 4001:4250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings4001.4250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 4251:4500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings4251.4500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 4501:4750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings4501.4750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 4751:5000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings4751.5000.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 5001:5250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings5001.5250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 5251:5500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings5251.5500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 5501:5750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings5501.5750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 5751:6000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings5751.6000.csv', delim = ',')

## sightings 6001-8000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 6001:6250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings6001.6250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 6251:6500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings6251.6500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 6501:6750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings6501.6750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 6751:7000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings6751.7000.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 7001:7250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings7001.7250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 7251:7500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings7251.7500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 7501:7750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings7501.7750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 7751:8000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings7751.8000.csv', delim = ',')

## sightings 8001-10000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 8001:8250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings8001.8250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 8251:8500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings8251.8500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 8501:8750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings8501.8750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 8751:9000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings8751.9000.csv', delim = ',')




gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 9001:9250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings9001.9250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 9251:9500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings9251.9500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 9501:9750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings9501.9750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 9751:10000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings9751.10000.csv', delim = ',')

## sightings 10001-12000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 10001:10250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings10001.10250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 10251:10500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings10251.10500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 10501:10750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings10501.10750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 10751:11000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings10751.11000.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 11001:11250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings11001.11250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 11251:11500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings11251.11500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 11501:11750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings11501.11750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 11751:12000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings11751.12000.csv', delim = ',')

## sightings 12001-14000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 12001:12250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings12001.12250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 12251:12500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings12251.12500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 12501:12750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings12501.12750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 12751:13000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings12751.13000.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 13001:13250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings13001.13250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 13251:13500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings13251.13500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 13501:13750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings13501.13750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 13751:14000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings13751.14000.csv', delim = ',')

## sightings 14001-16000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 14001:14250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings14001.14250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 14251:14500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings14251.14500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 14501:14750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings14501.14750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 14751:15000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings14751.15000.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 15001:15250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings15001.15250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 15251:15500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings15251.15500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 15501:15750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings15501.15750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 15751:16000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings15751.16000.csv', delim = ',')

## sightings 16001-18000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 16001:16250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings16001.16250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 16251:16500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings16251.16500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 16501:16750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings16501.16750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 16751:17000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings16751.17000.csv', delim = ',')



gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 17001:17250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings17001.17250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 17251:17500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings17251.17500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 17501:17750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings17501.17750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 17751:18000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings17751.18000.csv', delim = ',')

## sightings 18001-20000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 18001:18250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings18001.18250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 18251:18500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings18251.18500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 18501:18750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings18501.18750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 18751:19000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings18751.19000.csv', delim = ',')



gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 19001:19250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings19001.19250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 19251:19500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings19251.19500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 19501:19750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings19501.19750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 19751:20000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings19751.20000.csv', delim = ',')

## sightings 20001-22000 ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 20001:20250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings20001.20250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 20251:20500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings20251.20500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 20501:20750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings20501.20750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 20751:21000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings20751.21000.csv', delim = ',')



gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 21001:21250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings21001.21250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 21251:21500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings21251.21500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 21501:21750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings21501.21750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 21751:22000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings21751.22000.csv', delim = ',')

## sightings 22001-end ####
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 22001:22250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings22001.22250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 22251:22500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings22251.22500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 22501:22750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings22501.22750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 22751:23000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings22751.23000.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 23001:23250) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings23001.23250.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 23251:23500) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings23251.23500.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 23501:23750) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings23501.23750.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 23751:24000) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings23751.24000.csv', delim = ',')

nrow(gbi_matrix) # 24174
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 24001:nrow(gbi_matrix)) {
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
write_delim(gbi_df, 'data_processed/anp_bayesian_allpairwiseevents_22.06.22_sightings24001.24174.csv', delim = ',')


















######
### add elephant ID numbers to assigned index factors
str(gbi_df)
gbi_id <- data.frame(id_1 = colnames(gbi_matrix), node_1 = as.numeric(1:472))
gbi_check <- left_join(x = gbi_df, y = gbi_id, by = 'node_1')
gbi_id <- data.frame(id_2 = colnames(gbi_matrix), node_2 = as.numeric(1:472))
gbi_check <- left_join(x = gbi_check, y = gbi_id, by = 'node_2')

### correct obs_id to match encounter numbers in other spreadsheets (some encounters have missing data: encounter numbers 1,2,3,5,6,8... where obs_id 1,2,3,4,5,6...). Warning: sometimes this step can exceed vector memory and not run.
eles$obs_id <- as.integer(as.factor(eles$encounter)) ; eles <- eles[,c(1,13,2:12)]
gbi_encounter <- left_join(x = gbi_check, y = eles, by = 'obs_id')
length(unique(gbi_encounter$encounter)) # 574 -- correct

### remove duplicate rows where a dyad is recording as being observed together twice during the same sighting
gbi_encounter$unique <- paste(gbi_encounter$node_1, gbi_encounter$node_2, gbi_encounter$social_event, gbi_encounter$obs_id, sep = '_')
gbi_distinct <- dplyr::distinct(gbi_encounter) # 40708781 obs
colnames(gbi_distinct) # "node_1","node_2","social_event","obs_id","id_1","id_2","encounter","elephant","date","time","location","gps_s","gps_e",'herd_type","total_elephants_numeric","total_elephants_uncert","total_id_hkm","perc_id_hkm","unique"
head(gbi_distinct, 10)
gbi_distinct <- gbi_distinct[,c(1:7,9:18)]
colnames(gbi_distinct)[c(7,16:17)] <- c('encounter_id','total_id','perc_id')

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power
gbi_distinct$dyad <- paste(gbi_distinct$id_1, gbi_distinct$id_2, sep = '_')
gbi_distinct$dyad_id <- as.integer(as.factor(gbi_distinct$dyad))
gbi_distinct$location_id <- as.integer(as.factor(gbi_distinct$location))
gbi_distinct <- dplyr::distinct(gbi_distinct)
head(gbi_distinct)

### convert to Binomial model data format -- aggregate all sightings of each dyad together into a count
df_agg <- gbi_distinct %>%
  group_by(id_1, id_2) %>%
  summarise(event_count=sum(social_event), dyad_id=cur_group_id()) %>%
  mutate(node_1_id=as.integer(as.factor(id_1)), node_2_id=as.integer(as.factor(id_2)))
length(df_agg$id_1) == cumsum(1:471)[471] # check have correct number of dyads -- number will be the (n-1)th value of the triangular number sequence in which n = total number of elephants in analysis (472). If TRUE, correct number of pairs.
head(df_agg) ; tail(df_agg)
##   id_1  id_2  event_count dyad_id node_1_id node_2_id
##   <chr> <chr>       <dbl>   <int>     <int>     <int>
## 1 F1    F10             0       1         1         1  -- both F1 and F10 have registered as elephant number 1
## 2 F1    F100            0       2         1         2
## 3 F1    F101            0       3         1         3
## 4 F1    F102            0       4         1         4
## 5 F1    F103            0       5         1         5
## 6 F1    F104            0       6         1         6  -- F10, F100-F104 never interacted with F1
## | |     |               |       |         |         |
## | |     |               |       |         |         |
## | |     |               |       |         |         |
## 1 U67   U7              0  111151         1         1  -- U67 and U7 both registering as elephant number 1 (like F1 and F10 above)
## 2 U67   U8              0  111152         1         2
## 3 U67   U9              0  111153         1         3
## 4 U7    U8              1  111154         1         1
## 5 U7    U9              1  111155         1         2
## 6 U8    U9              2  111156         1         1
# All good except node_1_id and node_2_id reset to 1 for every new value of id_1, so node_1_id contains nothing but "1" in every cell, and node_2_id counts all values when F1 is node_1, all-1 for F10, all-2 for F100..., only U8 and U9 when U7 is node_1, and only U9 when U8 is node_1

### correct values in node_1_id and node_2_id using factor values.
df_agg <- df_agg[,c(1:4)]
df_agg$node_1 <- as.integer(as.factor(df_agg$id_1))
df_agg$node_2 <- as.integer(as.factor(df_agg$id_2))+1 # add 1 so starts at 2 and F1 is "1"
head(df_agg,10) ; tail(df_agg,10)

### add data about nodes
colnames(nodes)
nodes <- nodes[,c(1,3:5,9,11:13)]
nodes$id_1 <- nodes$id ; nodes$id_2 <- nodes$id
colnames(nodes) ; colnames(df_agg)
dyads <- left_join(x = df_agg, y = nodes, by = 'id_1')
colnames(dyads)[c(2,7:15)] <- c('id_2','id_pad_1','name_1','age_class_1','age_category_1','sex_1','id1_deletecolumn','count_1','dem_class_1','deletecolumn1')
dyads <- left_join(x = dyads, y = nodes, by = 'id_2')
colnames(dyads)[c(1,16:24)] <- c('id_1','id_pad_2','name_2','age_class_2','age_category_2','sex_2','id2_deletecolumn','count_2','dem_class_2','deletecolumn2')
dyads <- dyads[,c(4,1,2,5,6,3,7,16,8,17,9,18,10,19,11,20,13,22,14,23)]
head(dyads)

### remove any elephants from whom their is disagreement in the different data frames regarding their names or ID numbers
unique(gbi_distinct$id_1) # 155 females, 250 males, 65 unknowns
unique(nodes$id_1)        # 151 females, 245 males, 66 unknowns
length(which(is.na(dyads$id_pad_1)))                 # 2639 entries where elephants have no information
unique(dyads$id_1[which(is.na(dyads$id_pad_1))])     # "F157" "F158" "F176" "M125" "M13"  "M138" "M21"  "M223" "M227"
length(which(is.na(dyads$id_pad_2)))                 # 1600 entries where elephants have no information
unique(dyads$node_2[which(is.na(dyads$id_pad_1))])   # 412 elephants -- all individuals that come after F157 in the sequence
# all of these should be removed -- these are sightings of elephants which were deleted previously from the ele_nodes and ele_links data frames due to very high uncertainty in their identity, and so their sightings are unreliable.

dyads <- dyads[dyads$id_1 != "F157" & dyads$id_1 != "F158" & dyads$id_1 != "F176" & 
                 dyads$id_1 != "M125" & dyads$id_1 != "M13" & dyads$id_1 !=  "M138" & 
                 dyads$id_1 != "M21" & dyads$id_1 != "M223" & dyads$id_1 != "M227", ]
dyads <- dyads[dyads$id_2 != "F157" & dyads$id_2 != "F158" & dyads$id_2 != "F176" & 
                 dyads$id_2 != "M125" & dyads$id_2 != "M13" & dyads$id_2 !=  "M138" & 
                 dyads$id_2 != "M21" & dyads$id_2 != "M223" & dyads$id_2 != "M227", ]
length(which(is.na(dyads$id_pad_1)))                 # 0 entries where elephants have no information
length(which(is.na(dyads$id_pad_2)))                 # 0 entries where elephants have no information

### write csv
readr::write_delim(dyads, 'data_processed/motnp_bayesian_trimmedpairwiseevents_22.01.10.csv', delim = ',')

