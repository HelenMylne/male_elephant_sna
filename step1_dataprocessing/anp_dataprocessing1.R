# Amboseli data
#### Information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022
#### Set up ####
#library(tidyverse, lib.loc = 'packages')
#library(lubridate, lib.loc = 'packages')
#library(janitor, lib.loc = 'packages')
#library(hms, lib.loc = 'packages')
#library(readxl, lib.loc = 'packages')
#library(data.table, lib.loc = 'packages')
#library(spatsoc, lib.loc = 'packages')

library(tidyverse)
library(lubridate)
library(janitor)
library(hms)
library(readxl)
library(data.table)
library(spatsoc)

#### Import sightings data -- remove final two columns and convert two to character just to ensure exactly the same as previously when it was working, other than the 13 removed rows ####
ate <- readxl::read_xlsx('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ATE_MaleSightingsCleaned_Fishlock220808.xlsx') %>% 
  janitor::clean_names()# %>% 
  #select(-obs_casename, -row_num)

old <- readxl::read_xlsx('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/old_anp/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% janitor::clean_names()

str(ate)
str(old)

which(ate$casename != ate$bull_id)

ate$bull_q_c <- as.character(ate$bull_q_c)
ate$grp_q_c <- as.character(ate$grp_q_c)

locations <- old[,c('obs_id','casename','obs_num','utm_lat','utm_long','bull_q_r','obs_type')]
ate <- left_join(ate, locations, by = c('obs_id','casename'))

ate <- ate[,c('obs_id','casename','musth','obs_date','obs_time','obs_num','grid_code','utm_lat','utm_long','hab_code_1','hab_code_2','bull_q_r','obs_type','grp_q_c','grp_size','bull_q_c','num_bulls','bulls_1_2','bulls_3_5','musth_male','oestrus_fem','act_code')]

colnames(ate) == colnames(old)
colnames(ate)[c(6,8,9,12,13)] <- c("obs_num_old","utm_lat_old","utm_long_old","bull_q_r_old","obs_type_old")

rm(locations,old)

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
test <- ate[,c('obs_id', 'obs_date', 'correct_time_hms', 'obs_num_old', 'grp_size')] %>% 
  distinct()
lu <- function(x) { length(unique(x)) }
ate_nums <- tapply(X = ate$obs_num_old, INDEX = ate$obs_date, FUN = lu )
test_nums <- tapply(X = test$obs_num_old, INDEX = test$obs_date, FUN = lu )
which(ate_nums != test_nums)
ate_nums[1:50] ; test_nums[1:50]

table(ate$obs_date[which(ate$obs_num_old == '0')])
table(ate$obs_date[which(ate$obs_num_old == '00')])
ate$obs_num_old <- ifelse(ate$obs_num_old == '0','00', ate$obs_num_old)
table(ate$obs_num_old)

table(ate$obs_date[which(ate$obs_num_old == '0a')])
table(ate$obs_date[which(ate$obs_num_old == '0A')])
ate$obs_num_old <- ifelse(ate$obs_num_old == '0a','0A', ate$obs_num_old)
table(ate$obs_num_old)

table(ate$obs_date[which(ate$obs_num_old == '0b')])
table(ate$obs_date[which(ate$obs_num_old == '0B')])
ate$obs_num_old <- ifelse(ate$obs_num_old == '0b','0B', ate$obs_num_old)
table(ate$obs_num_old)

table(test$obs_date[which(test$obs_num_old == '1')])
table(test$obs_date[which(test$obs_num_old == '01')])
unique(ate$obs_num_old[which(ate$obs_date == '2020-08-28')])
ate$obs_num_old <- ifelse(ate$obs_num_old == '1','01', ate$obs_num_old)
table(ate$obs_num_old)

ate$obs_num_old_std <- NA
for(i in 1:nrow(ate)){
  date_row <- ate[ate$obs_date == ate$obs_date[i],]
  date_row$obs_num_old_std <- as.integer(as.factor(date_row$obs_num_old))
  ate$obs_num_old_std[i] <- date_row$obs_num_old_std[which(date_row$obs_id == ate$obs_id[i])[1]]
}

table(ate$obs_num_old_std)

## utm_lat and utm_long = There is no mask applied to GPS at the moment, so this needs checking for outliers and impossible values, as well as note that short values are possible -- try plotting and see where they are
summary(ate$utm_lat_old)
summary(ate$utm_long_old)
no_gps <- unique(ate$obs_date[which(ate$utm_lat_old == 0)])
gps <- unique(ate$obs_date[which(ate$utm_lat_old != 0)])
min(gps)
no_gps[4100:5023] # many that are after the start of using actual GPS

which(ate$utm_lat_old == 0 & ate$utm_long_old != 0) # 0 have longitude value but no latitude
which(ate$utm_lat_old != 0 & ate$utm_long_old == 0) # 5 have latitude value but no longitude

plot(utm_long_old ~ utm_lat_old, data = ate[ate$utm_lat_old != 0 & ate$utm_long_old != 0,])

summary(ate$utm_long_old[which(ate$utm_lat_old != 0 & ate$utm_long_old != 0)]) # should all be around 37
summary(ate$utm_lat_old[which(ate$utm_lat_old != 0 & ate$utm_long_old != 0)])  # should all be around -2

plot(utm_long_old ~ utm_lat_old, data = ate[ate$utm_lat_old != 0 & ate$utm_long_old != 0,],
     xlim = c(1,3), ylim = c(30,40))         # nothing in the range it should actually be

plot(utm_long_old ~ utm_lat_old, data = ate[ate$utm_lat_old != 0 & ate$utm_long_old != 0,],
     xlim = c(0,100000), ylim = c(0,400000)) # large mass of points but doesn't seem to be even an order of magnitude wrong to explain it

## hab_code_1 = habitat code for most of the group: 0 = Not described, 01 = short grass plain, 02 = Consimilis tall grasslands, 03 = Acacia tortilis woodland, 04 =	Acacia xanthophloea woodland, 05 = Salvadora/sueda, 06 = Palm woodland, 07 = Swamp edge woodland, 08 = Swamp edge, 09 =	Swamp, A =	Bushed grassland, B =	Open bushland north, C = Open bushland south, D	= Dense bushland north, E =	Acacia nubica, F = Dense bushland south, U = Unknown
## hab_code_2 = second habitat code if >1 category
sort(unique(ate$hab_code_1))
sort(unique(ate$hab_code_2))

## bull_q_r = quality of recognition for males (3 = knew all, 2 = knew > half, 1 = knew < half, U = unknown, blanks relate to data edits I haven't completed in the base data)
table(ate$bull_q_r_old)

## obs_type = Group Type: B = male only, M = males and females, U = unknown
table(ate$obs_type_old)

## grp_q_c = Group Quality of count: 3 = exact, 2 = good estimate, 1 = estimate, 0 = no estimate
table(ate$grp_q_c)

## grp_size = N elephants encountered in group; -1 is missing data. Zeros are typos that I haven't had chance to correct at present
table(ate$grp_size)

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

## progress report
print(paste0('sightings data imported at ', Sys.time()))

#### Import nodes data ####
nodes_old <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/old_anp/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()

nodes <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ATE_LifeHistoryData_Fishlock220808.xlsx') %>% janitor::clean_names() %>% distinct()

colnames(nodes_old)
colnames(nodes)

## make character string of ID number
nodes$id <- paste('M', nodes$casename, sep = '')

## check ID numbers match in both data frames
nodes_id <- sort(unique(nodes$id))
ate_id <- sort(unique(ate$id))
length(nodes_id) ; length(ate_id) ## DO NOT MATCH, OR EVEN CLOSE!

## progress report
print(paste0('nodes data imported at ', Sys.time()))

#### create group-by-individual matrix ####
ate$obs_id_std <- as.integer(as.factor(ate$obs_id))
ate_asnipe <- ate[,c(4,27)]
colnames(ate_asnipe) <- c('ID','group')
ate_asnipe <- data.table::setDT(ate_asnipe) # just converts to a data table, no other change. 
gbi_matrix <- spatsoc::get_gbi(DT = ate_asnipe, group = 'group', id = 'ID')
# NOTE: THE ORDER OF SIGHTINGS IS ACCORDING TO ate$obs_id_std NOT ate$obs_id -- USE OBS_ID_STD TO MATCH UP SIGHTING INFORMATION TO INDIVIDUALS

## progress report
print(paste0('gbi_matrix created at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings1.250.csv', delim = ',')

## progress report
print(paste0('sightings 1-250 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings251.500.csv', delim = ',')

## progress report
print(paste0('sightings 251-500 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings501.750.csv', delim = ',')

## progress report
print(paste0('sightings 501-750 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings751.1000.csv', delim = ',')

## progress report
print(paste0('sightings 751-1000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings1001.1250.csv', delim = ',')

## progress report
print(paste0('sightings 1001-1250 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings1251.1500.csv', delim = ',')

## progress report
print(paste0('sightings 1251-1500 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings1501.1750.csv', delim = ',')

## progress report
print(paste0('sightings 1501-1750 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings1751.2000.csv', delim = ',')

## progress report
print(paste0('sightings 1751-2000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings2001.2250.csv', delim = ',')

## progress report
print(paste0('sightings 2001-2250 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings2251.2500.csv', delim = ',')

## progress report
print(paste0('sightings 2251-2500 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings2501.2750.csv', delim = ',')

## progress report
print(paste0('sightings 2501-2750 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings2751.3000.csv', delim = ',')

## progress report
print(paste0('sightings 2751-3000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings3001.3250.csv', delim = ',')

## progress report
print(paste0('sightings 3001-3250 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings3251.3500.csv', delim = ',')

## progress report
print(paste0('sightings 3251-3500 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings3501.3750.csv', delim = ',')

## progress report
print(paste0('sightings 3501-3750 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings3751.4000.csv', delim = ',')

## progress report
print(paste0('sightings 3751-4000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings4001.4250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings4251.4500.csv', delim = ',')

## progress report
print(paste0('sightings 4001-4500 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings4501.4750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings4751.5000.csv', delim = ',')

## progress report
print(paste0('sightings 4501-5000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings5001.5250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings5251.5500.csv', delim = ',')

## progress report
print(paste0('sightings 5001-5500 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings5501.5750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings5751.6000.csv', delim = ',')

## progress report
print(paste0('sightings 5501-6000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings6001.6250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings6251.6500.csv', delim = ',')

## progress report
print(paste0('sightings 6001-6500 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings6501.6750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings6751.7000.csv', delim = ',')

## progress report
print(paste0('sightings 6501-7000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings7001.7250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings7251.7500.csv', delim = ',')

## progress report
print(paste0('sightings 7001-7500 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings7501.7750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings7751.8000.csv', delim = ',')

## progress report
print(paste0('sightings 7501-8000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings8001.8250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings8251.8500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings8501.8750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings8751.9000.csv', delim = ',')

## progress report
print(paste0('sightings 8001-9000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings9001.9250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings9251.9500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings9501.9750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings9751.10000.csv', delim = ',')

## progress report
print(paste0('sightings 9001-10000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings10001.10250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings10251.10500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings10501.10750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings10751.11000.csv', delim = ',')

## progress report
print(paste0('sightings 10001-11000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings11001.11250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings11251.11500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings11501.11750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings11751.12000.csv', delim = ',')

## progress report
print(paste0('sightings 11001-12000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings12001.12250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings12251.12500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings12501.12750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings12751.13000.csv', delim = ',')

## progress report
print(paste0('sightings 12001-13000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings13001.13250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings13251.13500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings13501.13750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings13751.14000.csv', delim = ',')

## progress report
print(paste0('sightings 13001-14000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings14001.14250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings14251.14500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings14501.14750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings14751.15000.csv', delim = ',')

## progress report
print(paste0('sightings 14001-15000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings15001.15250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings15251.15500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings15501.15750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings15751.16000.csv', delim = ',')

## progress report
print(paste0('sightings 15001-16000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings16001.16250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings16251.16500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings16501.16750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings16751.17000.csv', delim = ',')

## progress report
print(paste0('sightings 16001-17000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings17001.17250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings17251.17500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings17501.17750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings17751.18000.csv', delim = ',')

## progress report
print(paste0('sightings 17001-18000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings18001.18250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings18251.18500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings18501.18750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings18751.19000.csv', delim = ',')

## progress report
print(paste0('sightings 18001-19000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings19001.19250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings19251.19500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings19501.19750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings19751.20000.csv', delim = ',')

## progress report
print(paste0('sightings 19001-20000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings20001.20250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings20251.20500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings20501.20750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings20751.21000.csv', delim = ',')

## progress report
print(paste0('sightings 20001-21000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings21001.21250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings21251.21500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings21501.21750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings21751.22000.csv', delim = ',')

## progress report
print(paste0('sightings 21001-22000 completed at ', Sys.time()))

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings22001.22250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings22251.22500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings22501.22750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings22751.23000.csv', delim = ',')

## progress report
print(paste0('sightings 22001-23000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings23001.23250.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings23251.23500.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings23501.23750.csv', delim = ',')

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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings23751.24000.csv', delim = ',')

## progress report
print(paste0('sightings 23001-24000 completed at ', Sys.time()))

## next section
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/anp_bayesian_allpairwiseevents_cleanedrawdata_sightings24001.24174.csv', delim = ',')

## progress report
print(paste0('sightings 24001-end completed at ', Sys.time()))
