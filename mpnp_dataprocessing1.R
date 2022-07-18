### Bayesian analysis of EfA data ####
# Script to process association data from Makgadikgadi Pans National Park, Botswana.
# Data collected: 1st November 2019 to 5th August 2021
# Collected by: Elephants for Africa (Dr Kate Evans)
# Data supplied by: EfA and Dr Kate Evans, 19th October 2021

### Set up ####
# load packages
library(tidyverse)
library(dplyr)
library(rstan)
library(rethinking)
library(igraph)
library(dagitty)
library(cmdstanr)
library(lubridate)

# information
sessionInfo()
R.Version()
rstan::stan_version()
#packageVersion('')
#citation('')

# set stan path
set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')

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
#### Sightings data frame ####
# sightings data
s <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211214.xlsx')
str(s)
colnames(s)[c(1:23,57)] <- s[2,c(1:23,57)]
colnames(s)[24:56] <- c('CM','CF','CU','CM','CF','CU','JM','JF','JU','YPM','YPF','YPU','OPM','OPF','OPU',
                          'YAM','YAF','YAU','MAM','MAF','MAU','OAM','OAF','OAU','UM','UF','UU','SM','SF','SU',
                          'AM','AF','AU')
s <- s[3:nrow(s),]
s <- janitor::clean_names(s)

check_na <- data.frame(col = 24:56,del = NA)
for(i in 1:nrow(check_na)) {
  check_na$del[i] <- nrow(s) - length(which(is.na(s[,i]) == 'FALSE'))
} # no columns to delete -- all have at least one recording
rm(check_na)

# individual data
efa <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')
str(efa)
efa$time_cat <- lubridate::hour(efa$Time)
efa <- separate(efa, Time, into = c('wrong_date','time'), sep = ' ')
efa$time <- hms::as_hms(efa$time)
efa <- efa[,c(1:5,7,33,8:32)]
efa$Age_Range_ID <- as.factor(efa$Age_Range_ID)
efa$Activity_ID <- as.factor(efa$Activity_ID)
efa$Distance_To_Observer <- as.numeric(efa$Distance_To_Observer)
efa$Physical_Condition_ID <- as.factor(efa$Physical_Condition_ID)
efa <- efa[,c(1:15,18:27,31:32)]
efa$Date <- lubridate::as_date(efa$Date)
efa <- janitor::clean_names(efa)

## taking a look ####
table(efa$elephant_id)
hist(table(efa$elephant_sighting_id), breaks = 30, las = 1)
length(which(table(efa$elephant_sighting_id) == 1)) # 2400 elephants seen singly
length(which(table(efa$elephant_sighting_id) == 2)) # 1017 elephants seen as pairs

table(efa$age_range_id)

summary(efa$reaction_indices)
efa$reaction_indices <- factor(efa$reaction_indices, levels = c(1:3))
table(efa$reaction_indices)
barplot(table(efa$reaction_indices), xlab = 'reaction index score', ylab = 'count',
        las = 1, ylim = c(0,500), xlim = c(0,4))
abline(h = 0, lwd = 2) ; axis(2, at = seq(0,500,50), las = 1)
text('No response\n= 349', x = 0.7, y = 480)
text('Calm response\n= 448', x = 1.9, y = 480)
text('Agitated, charge\nor flee = 0', x = 3.1, y = 480)

range(efa$date) # "2012-05-04" "2021-09-24"

prop.table(table(efa$sex_id))
#           1            2            3 
#0.9989531145 0.0006729978 0.0003738877   <-- 99.9% 1 = male, rest = BH/MX

## restructuring ####
str(efa)
#$ id_elephant_visual        : num [1:13429] 15993 15994 15995 15996 15997 ...
#$ elephant_sighting_id      : num [1:13429] 8552 8552 8552 8552 8552 ...
#$ elephant_id               : chr [1:13429] "T6409" "-" "-" "-" ...
#$ in_the_field_elephant_code: chr [1:13429] "A" "B" "C" "D" ...
#$ date                      : Date[1:13429], format: "2021-09-24" "2021-09-24" "2021-09-24" "2021-09-24" ...
#$ time                      : 'hms' num [1:13429] 10:15:00 10:15:00 10:15:00 10:15:00 ...
#$ time_cat                  : int [1:13429] 10 10 10 10 10 10 10 10 10 10 ...
#$ age_range_id              : Factor w/ 8 levels "10","2","3","4",..: 5 6 5 5 5 5 6 6 6 7 ...
#$ age                       : logi [1:13429] NA NA NA NA NA NA ...
#$ activity_id               : Factor w/ 13 levels "0","1","10","2",..: 11 11 11 11 11 11 11 11 11 11 ...
#$ in_musth                  : chr [1:13429] "N" "N" "N" "N" ...
#$ social_id                 : num [1:13429] 3 3 3 3 3 3 3 3 3 3 ...
#$ reaction_indices          : num [1:13429] 1 1 1 1 1 1 1 1 1 1 ...
#$ distance_to_observer      : num [1:13429] 193 193 193 193 193 193 193 193 193 193 ...
#$ habitat_id                : chr [1:13429] "1.4" "1.4" "1.4" "1.4" ...
#$ pictures_taken            : logi [1:13429] TRUE FALSE FALSE FALSE FALSE FALSE ...
#$ sick                      : chr [1:13429] "N" "N" "N" "N" ...
#$ physical_condition_id     : Factor w/ 9 levels "1.5","2","2.5",..: 5 4 4 4 4 4 4 4 4 4 ...
#$ l_tgs_length_id           : chr [1:13429] "UK" "UK" "UK" "UK" ...
#$ l_tgs_width_id            : chr [1:13429] "UK" "UK" "UK" "UK" ...
#$ l_tgs_age_id              : chr [1:13429] "UK" "UK" "UK" "UK" ...
#$ r_tgs_length_id           : chr [1:13429] "0" "UK" "UK" "0" ...
#$ r_tgs_width_id            : chr [1:13429] "0" "UK" "UK" "0" ...
#$ r_tgs_age_id              : chr [1:13429] "0" "UK" "UK" "0" ...
#$ tg_swelling_id            : chr [1:13429] "0" "0" "0" "0" ...
#$ sex_id                    : num [1:13429] 1 1 1 1 1 1 1 1 1 1 ...
#$ notes                     : chr [1:13429] "Collared elephant" NA NA NA ...

efa_long <- data.frame(encounter = efa$elephant_sighting_id,
                       date = efa$date,
                       time = efa$time,
                       gps_s = NA,
                       gps_e = NA,
                       total_elephants = NA,
                       total_id = NA,
                       perc_id = NA,
                       type = ifelse(efa$sex_id == 1, 'MO', 'BH/MX'),
                       elephant = ifelse(efa$elephant_id == '-', NA, efa$elephant_id),
                       sex = ifelse(efa$sex_id == 1, 'M','F/U'),
                       age_range = efa$age_range_id)

for(i in 1:length(efa_long$total_elephants)) {
  efa_long$total_elephants[i] <- length(which(efa$elephant_sighting_id == efa_long$encounter[i]))
}

efa_id <- efa_long[!is.na(efa_long$elephant),]
for(i in 1:length(efa_id$total_id)) {
  efa_id$total_id[i] <- length(which(efa_id$encounter == efa_long$encounter[i]))
}

efa_id$perc_id <- 100*(efa_id$total_id/efa_id$total_elephants)

colnames(s)[1] <- 'encounter'
s$encounter <- as.numeric(s$encounter)
efa_gps <- left_join(x = efa_id, y = s, by = 'encounter') # first few rows are blank because 2 datasets formed at slightly different times -- group sightings only goes as far as 29th January 2021, individual sightings go until 24th September 2021
efa_gps <- efa_gps[,c(1,13,10,2,3,18,19,6:9,11,12)]
colnames(efa_gps)[4] <- c('date')

efa_gps$location <- paste(efa_gps$latitude, efa_gps$longitude, sep = '_')

write_delim(efa_gps, 'data_processed/mpnp_eles_long_22.03.08.csv', delim = ',')

#### Nodes data frame ####
ele_nodes <- data.frame(id_no = sort(unique(efa$elephant_id)),
                        age_class = NA,
                        age_category = NA,
                        sex = NA,
                        count = NA,
                        dem_class = NA)

efa$age_range_id[range(which(efa$elephant_id == ele_nodes$id_no[2]))]                                   # 7 and 8, levels 2-10
efa_age <- efa[c(2,3,8)]
efa_age$age_range_id <- as.numeric(efa$age_range_id)
str(efa_age)
efa_age$age_range_id[range(which(efa_age$elephant_id == ele_nodes$id_no[2]))]                           # 7 and 8, no other levels
ele_nodes$age_class[2] <- mean(efa_age$age_range_id[which(efa_age$elephant_id == ele_nodes$id_no[2])])  # mean = 7.055556
test <- efa[efa$elephant_id == 'B1000',] ; mean(as.numeric(test$age_range_id))                          # mean = 7.055556

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

ele_nodes$age_class[2] <- getmode(efa_age$age_range_id[which(efa_age$elephant_id == ele_nodes$id_no[2])])
test <- efa[efa$elephant_id == 'B1000',] ; getmode(as.numeric(test$age_range_id))

ele_nodes <- ele_nodes[2:6646,]

for(i in 1:length(ele_nodes$age_class)) {
  ele_nodes$age_class[i] <- getmode(efa_age$age_range_id[which(efa_age$elephant_id == ele_nodes$id_no[i])])
  ele_nodes$sex[i] <- length(unique(efa$sex_id[which(efa$elephant_id == ele_nodes$id_no[i])]))
  ele_nodes$count[i] <- length(which(efa$elephant_id == ele_nodes$id_no[i]))
}
check <- ele_nodes$id_no[which(ele_nodes$sex > 1)]
for(i in 1:length(check)){
  print(efa$sex_id[which(efa$elephant_id == check[i])])
} # for all warnings, the sex is 1 and NA -- NA is causing apparent presence of extra sex_id
for(i in 1:length(ele_nodes$age_class)) {
  ele_nodes$sex[i] <- unique(efa$sex_id[which(efa$elephant_id == ele_nodes$id_no[i])])
}
ele_nodes$sex <- ifelse(ele_nodes$sex == 1,'M','F')

ele_nodes$age_category <- ifelse(ele_nodes$age_class < 3, 'Calf',
                                 ifelse(ele_nodes$age_class == 3, 'Juvenile',
                                        ifelse(ele_nodes$age_class > 3 & ele_nodes$age_class <= 5, 'Pubescent',
                                               ifelse(ele_nodes$age_class > 5 & ele_nodes$age_class <= 8, 'Adult', 'UNK'))))
ele_nodes$age_class <- as.factor(ele_nodes$age_class)

ele_nodes$age <- ifelse(ele_nodes$age_category == 'Adult', 'A',
                        ifelse(ele_nodes$age_category == 'Pubescent', 'P',
                               ifelse(ele_nodes$age_category == 'Juvenile', 'J', 'C')))
ele_nodes$dem_class <- paste(ele_nodes$age, ele_nodes$sex, sep = '')
ele_nodes <- ele_nodes[c(1:6)]

write_delim(ele_nodes, 'data_processed/mpnp_elenodes_22.03.08.csv')

### clean environment
rm(efa, efa_age, efa_id, efa_long, test, check, i)

#### Mapping ####
plot(efa_gps$longitude ~ efa_gps$latitude)
efa_gps$longitude[which(efa_gps$longitude < 20)] <- NA
efa_gps$longitude[which(efa_gps$longitude < 21)] # 20.512689999999
efa_gps$longitude[which(efa_gps$latitude < -20)] # 24.750029999999999 24.707709999999999 24.707709999999999 24.768260000000001 24.764379999999999 24.765647999999999 24.528880000000001 24.528880000000001
plot(efa_gps$longitude ~ efa_gps$latitude)
plot(efa_gps$longitude ~ efa_gps$latitude, ylim = c(24.2,24.9), xlim = c(-20.8,-20.2),
     pch = 19, col = rgb(0,0,1,0.05))

#### create dyadic data frame ####
## create gbi_matrix ####
colnames(efa_gps)

### create group-by-individual matrix
eles_asnipe <- efa_gps[,c(4,5,3,14)]  # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- hour(eles_asnipe$time)*60*60 + minute(eles_asnipe$time)*60 + second(eles_asnipe$time) # convert time values to seconds through day

eles_asnipe <- eles_asnipe[,c(5,6,3,4)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
str(eles_asnipe)

# get_gbi generates a group by individual matrix. The function accepts a data.table with individual identifiers and a group column. The group by individual matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_')
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter))
max(eles_asnipe$group) # 3765
eles_asnipe <- eles_asnipe[,c(3,7)]
eles_asnipe <- data.table::setDT(eles_asnipe)
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')

## convert gbi_matrix to series of dyadic data frames ####
### code to convert gbi matrix format to dyadic data frame, shared by Prof Dan Franks and available from @JHart96 GitHub repository (https://github.com/JHart96/bison_examples/blob/main/examples/convert_gbi.md) -- NOTE: this step takes a long time to run, and the for loops get exponentially slower with more and more observations. Running it in chunks rather than all at once reduces how much it slows down.
nrow(gbi_matrix) # 3765 -- run in 14 blocks of 250 sightings + 1 block of 265
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1:250) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings1.250_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 251:500) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings251.500_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 501:750) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings501.750_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 751:1000) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings751.1000_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1001:1250) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings1001.1250_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1251:1500) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings1251.1500_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1501:1750) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings1501.1750_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1751:2000) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings1751.2000_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2001:2250) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings2001.2250_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2251:2500) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings2251.2500_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2501:2750) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings2501.2750_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2751:3000) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings2751.3000_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3001:3250) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings3001.3250_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3251:3500) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings3251.3500_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3501:3765) {
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, 'data_processed/mpnp_bayesian_allpairwiseevents_sightings3501.3765_22.03.08.csv', delim = ',')

## merge dyadic data frames, create aggregated count data frame ####
### read in all files
s250 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings1.250_22.03.08.csv')
s500 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings251.500_22.03.08.csv')
s750 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings501.750_22.03.08.csv')
s1000 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings750.1000_22.03.08.csv')

s1250 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings1001.1250_22.03.08.csv')
s1500 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings1251.1500_22.03.08.csv')
s1750 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings1501.1750_22.03.08.csv')
s2000 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings1750.2000_22.03.08.csv')

s2250 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings2001.2250_22.03.08.csv')
s2500 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings2251.2500_22.03.08.csv')
s2750 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings2501.2750_22.03.08.csv')
s3000 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings2750.3000_22.03.08.csv')

s3250 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings3001.3250_22.03.08.csv')
s3500 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings3251.3500_22.03.08.csv')
s3765 <- read_csv('data_processed/mpnp_bayesian_allpairwiseevents_sightings3501.3765_22.03.08.csv')

gbi_df <- rbind(s250, s500, s750, s1000, s1250, s1500, s1750, s2000,
                s2250, s2500, s2750, s3000, s3250, s3500, s3765) %>% 
  distinct()
rm(s250, s500, s750, s1000, s1250, s1500, s1750, s2000, s2250, s2500, s2750, s3000, s3250, s3500, s3765)

### add elephant ID numbers to assigned index factors
str(gbi_df)
gbi_id <- data.frame(id_1 = colnames(gbi_matrix), node_1 = as.numeric(1:472))
gbi_check <- left_join(x = gbi_df, y = gbi_id, by = 'node_1')
gbi_id <- data.frame(id_2 = colnames(gbi_matrix), node_2 = as.numeric(1:472))
gbi_check <- left_join(x = gbi_check, y = gbi_id, by = 'node_2')

### correct obs_id to match encounter numbers in other spreadsheets (some encounters have missing data: encounter numbers 1,2,3,5,6,8... where obs_id 1,2,3,4,5,6...). Warning: sometimes this step can exceed vector memory and not run.
efa_gps$obs_id <- as.integer(as.factor(efa_gps$encounter)) ; efa_gps <- efa_gps[,c(1,13,2:12)]
gbi_encounter <- left_join(x = gbi_check, y = efa_gps, by = 'obs_id')
length(unique(gbi_encounter$encounter)) # 574 -- correct

### remove duplicate rows where an dyad is recording as being observed together twice during the same sighting
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
readr::write_delim(dyads, 'data_processed/mpnp_bayesian_trimmedpairwiseevents_22.01.13.csv', delim = ',')















#####
### import data for aggregated model (binomial)
counts_df <- read_delim('data_processed/mpnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',')

# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
counts_df$sex_1 <- as.character(counts_df$sex_1)
str(counts_df)  # sex_1 still comes up as logical?

# create variable for age difference
unique(counts_df$age_category_1) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
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
counts_df[is.na(counts_df$age_cat_id_1),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_cat_id_1[which(is.na(counts_df$age_cat_id_1))] <- 1

unique(counts_df$age_category_1[counts_df$age_class_1 == 'Calf'])      # shouldn't include any ages over 4-5
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Adult'])     # shouldn't include any ages under 20-25

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
counts_df[is.na(counts_df$age_cat_id_2),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_cat_id_2[which(is.na(counts_df$age_cat_id_2))] <- 1

counts_df$age_class_2 <- ifelse(counts_df$age_cat_id_2 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_2 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_2 > 4, 'Adult','Pubescent')))

unique(counts_df$age_category_2[counts_df$age_class_2 == 'Calf'])      # shouldn't include any ages over 4-5
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Adult'])     # shouldn't include any ages under 20-25

### correct dem_class with corrected age classes
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))

### add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$all_events

### create data list -- can contain no NA values in any column, even if column is not specified in model
counts_ls <- list(
  n_dyads  = nrow(counts_df),          # total number of times one or other of the dyad was observed
  together = counts_df$all_events,     # count number of sightings seen together
  apart    = counts_df$apart)          # count number of sightings seen apart

################ 3) Create simulated data set ################
### Population
# 120 individuals: 50 male (10 each of ages 3-7), 50 female (10 each of ages 3-7), 20 unknown (10 each of ages 1-2)
# males no preference for any other individual -- edge weight 0.2
# females in 8 cliques -- edge weight 0.8, otherwise 0.2
# unknowns attached to a single female -- edge weight 0.95

# assign females to family units
cliques_F <- sample(1:8,50,replace=T)
table(cliques_F)

# create data frame of all individuals
population <- data.frame(id = rep(NA,120),
                         num = c(1:50,1:50,1:20),
                         sex = c(rep('M',50), rep('F',50), rep('U',20)),
                         age = c(rep(c(3:7), each = 10), rep(c(3:7), each = 10), rep(c(1:2), each = 10)),
                         clique = c(rep(NA, 50), cliques_F, rep(NA, 20)))
population$id <- paste(population$sex, population$num, sep = '')

# randomly assign calves to females
mothers <- data.frame(id = population$id[101:120],
                      family = sample(population$id[61:100], replace = F, size = 20))
mothers
population <- left_join(x = population, y = mothers, by = 'id')

# need a clique for U elephants and a family for F elephants
females <- population[population$sex == 'F', c(1,5,6)]
unknown <- population[population$sex == 'U', c(1,5,6)]
females$mother <- females$id ; unknown$mother <- unknown$family   # create variable to join by
families <- left_join(x = unknown, y = females, by = 'mother')    # join to give unknowns a clique
unknown <- families[,c(1,4,6)]

families <- left_join(x = females, y = unknown, by = 'mother')    # join to give mothers their calf
females <- families[,c(1,5,2)]

colnames(unknown) <- c('id','mother','clique')
colnames(females) <- c('id','mother','clique')
families <- rbind(females, unknown)                               # put mother and calf data into a single data frame

# join to combine male and female/unknown data together
population <- left_join(x = population, y = families, by = 'id')
population <- population[,c(1,3,4,7,8)]
colnames(population)[c(4,5)] <- c('family','clique')

### Create dyadic data frame
dyads <- data.frame(id_1 = rep(population$id, each = 120), id_2 = rep(population$id, 120))               # create new data frame of all dyads
dyads$same <- ifelse(dyads$id_1 == dyads$id_2, 'yes', 'no') ; dyads <- dyads[dyads$same == 'no', c(1,2)] # remove self-dyads (e.g. M1-M1)
dyads$node_1 <- as.integer(as.factor(dyads$id_1)) ; dyads$node_2 <- as.integer(as.factor(dyads$id_2))    # create factor of  <- 
dyads$dyad <- ifelse(dyads$node_1 < dyads$node_2,                               # create variable of dyads
                     paste(dyads$id_1, dyads$id_2, sep = '_'),
                     paste(dyads$id_2, dyads$id_1, sep = '_'))

dyads <- data.frame(dyad = unique(dyads$dyad))                                  # create new data frame of only unique dyads
dyads <- separate(dyads, col = dyad, into = c('id_1','id_2'), remove = FALSE)   # nodes columns
dyads$node_1 <- as.integer(as.factor(dyads$id_1))                               # nodes columns
dyads$node_2 <- as.integer(as.factor(dyads$id_2))                               # nodes columns
dyads$dyad_id <- as.integer(as.factor(dyads$dyad))                              # dyad index variable

population$id_1 <- population$id ; population$id_2 <- population$id  # variable to join information on each individual
dyads <- left_join(x = dyads, y = population, by = 'id_1')           # add information about node 1
dyads <- dyads[,c(1:6,8:11)]
colnames(dyads)[c(3,7:10)] <- c('id_2','sex_1','age_1','family_1','clique_1')

dyads <- left_join(x = dyads, y = population, by = 'id_2')           # add information about node 2
dyads <- dyads[,c(1:10,12:15)]
colnames(dyads)[c(2,11:14)] <- c('id_1','sex_2','age_2','family_2','clique_2')

dyads <- dyads[,c(1:7,11,8,12,9,13,10,14)]

### Assign type of dyad
dyads$pair_type <- ifelse(is.na(dyads$clique_1) | is.na(dyads$clique_2), 'no_group',
                          ifelse(dyads$clique_1 == dyads$clique_2,
                                 ifelse(is.na(dyads$family_1) | is.na(dyads$family_2), 'group',
                                        ifelse(dyads$family_1 == dyads$id_2 | dyads$family_2 == dyads$id_1,
                                               'family', 'group')),
                                 'no_group'))

### Assign edge weight to dyad pairs
for(i in 1:nrow(dyads)){
  if(dyads$pair_type[i] == 'family')  {dyads$edge[i] <- rethinking::rbeta2(1,0.95,50)   # mother-calf: edge weight around 0.95
  } else 
    if(dyads$pair_type[i] == 'group') {dyads$edge[i] <- rethinking::rbeta2(1,0.8,50)    # same family, not mother-calf: edge ~ 0.8
    } else 
      dyads$edge[i] <- rethinking::rbeta2(1,0.20,50)                                        # different family: edge weight around 0.2
}
summary(dyads$edge)

boxplot(dyads$edge ~ dyads$pair_type)                                                   # check simulations are around what I expected
points(y = dyads$edge, x = rep(2, length(which(dyads$pair_type == 'group'))),           # had a few problems with this one, but looks good now
       col = col.alpha('black',0.2), pch = 19)

### clean up environment
rm(families, females, mothers, unknown, cliques_F)

### Sample observations
N <- 20
dyads$event_count <- rbinom(nrow(dyads), N, prob = dyads$edge)

### plot simulated observations against assigned edge weight to check model input is correlates to true value
plot(event_count ~ edge, data = dyads,
     las = 1, ylab = 'sightings together', xlab = 'assigned edge weight',
     pch = 16, col = col.alpha(rangi2, 0.2))

### add additional columns (dyadic sex, age and dem_class)
head(dyads)
# assign each node to age category (adult, pubescent, juvenile or calf) from age class
dyads$age_cat_1 <- ifelse(dyads$age_1 == 1, 'C', ifelse(dyads$age_1 == 2, 'J', ifelse(dyads$age_1 < 5, 'P', 'A')))
dyads$age_cat_2 <- ifelse(dyads$age_2 == 1, 'C', ifelse(dyads$age_2 == 2, 'J', ifelse(dyads$age_2 < 5, 'P', 'A')))

# combine age categories into single variable for dyad
dyads$age_dyad <- ifelse(dyads$age_cat_1 == 'A',
                         ifelse(dyads$age_cat_2 == 'A', 'A_A',
                                ifelse(dyads$age_cat_2 == 'P', 'A_P',
                                       ifelse(dyads$age_cat_2 == 'J', 'A_J', 'A_C'))),
                         ifelse(dyads$age_cat_1 == 'P',
                                ifelse(dyads$age_cat_2 == 'A', 'A_P',
                                       ifelse(dyads$age_cat_2 == 'P', 'P_P',
                                              ifelse(dyads$age_cat_2 == 'J', 'P_J', 'P_C'))),
                                ifelse(dyads$age_cat_1 == 'J',
                                       ifelse(dyads$age_cat_2 == 'A', 'A_J',
                                              ifelse(dyads$age_cat_2 == 'P', 'P_J',
                                                     ifelse(dyads$age_cat_2 == 'J', 'J_J', 'J_C'))),
                                       ifelse(dyads$age_cat_2 == 'A', 'A_C',
                                              ifelse(dyads$age_cat_2 == 'P', 'P_C',
                                                     ifelse(dyads$age_cat_2 == 'J', 'J_C', 'C_C'))))))
unique(dyads$age_dyad) # "P_P" "A_P" "P_C" "P_J" "A_A" "A_C" "A_J" "C_C" "J_C" "J_J"
dyads$age_diff <- abs(dyads$age_1 - dyads$age_2)

# create composite measure of demography per node (e.g. AM = adult male, CU = calf of unknown sex)
dyads$dem_class_1 <- paste(dyads$age_cat_1, dyads$sex_1, sep = '')
dyads$dem_class_2 <- paste(dyads$age_cat_2, dyads$sex_2, sep = '')

# combine demographies into single variable for dyad
dyads$demf_1 <- as.integer(as.factor(dyads$dem_class_1)) ; dyads$demf_2 <- as.integer(as.factor(dyads$dem_class_2))
dyads$dem_dyad <- ifelse(dyads$demf_1 < dyads$demf_2,
                         paste(dyads$dem_class_1, dyads$dem_class_2, sep = '_'),
                         paste(dyads$dem_class_2, dyads$dem_class_1, sep = '_'))
dyads$dem_dyad <- ifelse(dyads$dem_dyad == 'CU_JU', 'JU_CU',           # standardise dem_dyad so always in same order
                         ifelse(dyads$dem_dyad == 'CU_PF', 'PF_CU',
                                ifelse(dyads$dem_dyad == 'CU_PM','PM_CU',
                                       ifelse(dyads$dem_dyad == 'JU_PF','PF_JU',
                                              ifelse(dyads$dem_dyad == 'JU_PM','PM_JU',dyads$dem_dyad)))))
dyads$dem_dyad_f <- as.integer(as.factor(dyads$dem_dyad))              # make index variable

dyads$dem_diff <- ifelse(dyads$dem_class_1 == dyads$dem_class_2, 0, 1) # binary: same or different
dyads$sex_dyad <- paste(dyads$sex_1, dyads$sex_2, sep = '_')           # composite sex variable for dyad
dyads$sex_diff <- ifelse(dyads$sex_1 == dyads$sex_2, 0, 1)             # binary: same or different
dyads <- dyads[,c(1:23,26:30)]                                         # remove demf variables -- not needed

# together vs apart per dyad
summary(dyads$event_count)             # together, all out of 20
dyads$apart <- N - dyads$event_count   # apart = 20-together

### create data list -- can contain no NA values in any column, even if column is not specified in model
simdat_ls <- list(
  n_dyads  = nrow(dyads),       # total number of times one or other of the dyad was observed
  together = dyads$event_count, # count number of sightings seen together
  apart    = dyads$apart)       # count number of sightings seen apart

# clear environment
rm(population, i, N)

################ 4) Define likelihood distributions, model and parameters to estimate ################
# Binomial model using a beta distribution where shape 1 = times together and shape 2 = times apart

### Compile Stan model
mod_1.1 <- cmdstan_model("models/simpleBetaNet_DWF_22.01.23.stan")
mod_1.1

################ 5) Run model on simulated data ################
# create data list
simdat_ls <- list(
  n_dyads  = nrow(dyads),        # Number of dyads
  together = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart         # Number of sightings of each dyad apart
)

### Fit model
edge_weight_1.1 <- mod_1.1$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
edge_weight_1.1$summary()
#   variable        mean      median       sd     mad          q5        q95     rhat ess_bulk ess_tail
#   <chr>          <dbl>       <dbl>    <dbl>   <dbl>       <dbl>      <dbl>    <dbl>    <dbl>    <dbl>
# 1 lp__         -83508.      -83508     61.5    60.8      -83608.    -83407.    1.00     1347.    1966.
# 2 weight[1]      0.183       0.173   0.0806  0.0791      0.0680       0.330    1.00     8950.    2451.
# 3 weight[2]      0.272       0.267   0.0914  0.0896       0.134       0.434    1.00     9029.    2347.
# 4 weight[3]      0.272       0.265   0.0928  0.0974       0.133       0.436    1.00     8872.    2522.
# 5 weight[4]      0.408       0.406   0.0987   0.102       0.251       0.572   0.999    10311.    2992.
# 6 weight[5]     0.0457      0.0322   0.0450  0.0325     0.00221       0.138    1.00     5098.    2204.
# 7 weight[6]      0.137       0.125   0.0727  0.0689      0.0393       0.276    1.00     8291.    2548.
# 8 weight[7]      0.179       0.169   0.0778  0.0760      0.0669       0.325    1.01     9580.    3062.
# 9 weight[8]      0.272       0.266   0.0932  0.0934       0.132       0.437    1.00     8875.    2823.
#10 weight[9]     0.0446      0.0314   0.0434  0.0326     0.00219       0.132    1.00     4795.    1798.
# … with 7,131 more rows
output1 <- read_cmdstan_csv(edge_weight_1.1$output_files()[1])
output2 <- read_cmdstan_csv(edge_weight_1.1$output_files()[2])
output3 <- read_cmdstan_csv(edge_weight_1.1$output_files()[3])
output4 <- read_cmdstan_csv(edge_weight_1.1$output_files()[4])

draws1_1.1 <- as.data.frame(output1$post_warmup_draws)
draws2_1.1 <- as.data.frame(output2$post_warmup_draws)
draws3_1.1 <- as.data.frame(output3$post_warmup_draws)
draws4_1.1 <- as.data.frame(output4$post_warmup_draws)
draws_simdat_1.1 <- rbind(draws1_1.1, draws2_1.1, draws3_1.1, draws4_1.1)
rm(output1, output2, output3, output4)

# check chain mixing -- well mixed
mean(draws1_1.1$`1.weight[1]`) ; mean(draws2_1.1$`1.weight[1]`) ; mean(draws3_1.1$`1.weight[1]`) ; mean(draws4_1.1$`1.weight[1]`)
mean(draws1_1.1$`1.weight[2]`) ; mean(draws2_1.1$`1.weight[2]`) ; mean(draws3_1.1$`1.weight[2]`) ; mean(draws4_1.1$`1.weight[2]`)
mean(draws1_1.1$`1.weight[3]`) ; mean(draws2_1.1$`1.weight[3]`) ; mean(draws3_1.1$`1.weight[3]`) ; mean(draws4_1.1$`1.weight[3]`)

# plot each chain individually to check mixing -- all mixed well
plot(draws1_1.1$`1.weight[1]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_1.1$`1.weight[1]`, col = 'green')
lines(draws3_1.1$`1.weight[1]`, col = 'blue')
lines(draws4_1.1$`1.weight[1]`, col = 'magenta')
plot(draws1_1.1$`1.weight[2]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_1.1$`1.weight[2]`, col = 'green')
lines(draws3_1.1$`1.weight[2]`, col = 'blue')
lines(draws4_1.1$`1.weight[2]`, col = 'magenta')
plot(draws1_1.1$`1.weight[3]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_1.1$`1.weight[3]`, col = 'green')
lines(draws3_1.1$`1.weight[3]`, col = 'blue')
lines(draws4_1.1$`1.weight[3]`, col = 'magenta')

# build traceplots -- highly variable, but generally fits the inputs
plot(draws1_1.1$`1.weight[1]`, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1_1.1$`1.weight[2]`, col = 'tan')
lines(draws1_1.1$`1.weight[3]`, col = 'orange')
lines(draws1_1.1$`1.weight[4]`, col = 'green')
lines(draws1_1.1$`1.weight[5]`, col = 'chocolate')
lines(draws1_1.1$`1.weight[6]`, col = 'blue')
lines(draws1_1.1$`1.weight[7]`, col = 'red')
lines(draws1_1.1$`1.weight[8]`, col = 'seagreen')
lines(draws1_1.1$`1.weight[9]`, col = 'purple')
lines(draws1_1.1$`1.weight[10]`,col = 'magenta')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[1]+1], col = 'black')       # traceplot of mother-calf relationships
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[2]+1], col = 'tan')         # +1 because first column is not linked to a dyad
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[3]+1], col = 'orange')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[4]+1], col = 'green')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[6]+1], col = 'blue')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[7]+1], col = 'red')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[9]+1], col = 'purple')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[10]+1], col = 'magenta')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[1]+1], col = 'black')      # traceplots of elephants in the same group but which are not mother-calf pairs
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[2]+1], col = 'tan')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[3]+1], col = 'orange')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[4]+1], col = 'green')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[6]+1], col = 'blue')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[7]+1], col = 'red')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[9]+1], col = 'purple')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[10]+1], col = 'magenta')

# mean and credible interval for each dyad weight
means_1.1 <- data.frame(dyad = dyads$dyad,
                        mean = apply(draws_simdat_1.1[2:7141], 2, mean), 
                        median = apply(draws_simdat_1.1[2:7141], 2, median))
draws_1.1_pi <- apply(draws_simdat_1.1[2:7141], 2, rethinking::PI)
means_1.1$pi_lwr <- draws_1.1_pi[1,]
means_1.1$pi_upr <- draws_1.1_pi[2,]

# boxplot comparing types of dyads
plot_data_1.1 <- left_join(x = means_1.1, y = dyads, by = 'dyad')
ggplot(data = plot_data_1.1, aes(x = pair_type, y = mean, fill = pair_type))+
  geom_boxplot()+
  geom_jitter(aes(x = pair_type, y = edge, fill = pair_type),
              pch = 4, colour = col.alpha('red', 0.2))+
  scale_fill_viridis_d()+
  theme_classic()

# clear large objects from environment
rm(draws_1.1_pi, draws1_1.1, draws1_1.1_pi, draws2_1.1, draws3_1.1, draws4_1.1, means_1.1, plot_data_1.1, draws_simdat_1.1, dyads, simdat_ls)


################ 6) Run model on real standardised data ################
### run model
mod_1.1

### Fit model (slow) - just over an hour (1:06:30)
weight_mpnp_1.1 <- mod_1.1$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_mpnp_1.1
# variable       mean     median     sd    mad         q5        q95 rhat ess_bulk ess_tail
#lp__      -819804.08 -819806.00 248.32 249.82 -820201.05 -819385.90 1.00     1300     2048
#weight[1]       0.05       0.04   0.03   0.03       0.01       0.12 1.00     6812     2442
#weight[2]       0.05       0.04   0.03   0.03       0.01       0.12 1.00     6335     2736
#weight[3]       0.07       0.06   0.05   0.04       0.01       0.16 1.00     6288     2513
#weight[4]       0.05       0.04   0.03   0.03       0.01       0.12 1.00     7594     2252
#weight[5]       0.07       0.06   0.05   0.04       0.01       0.17 1.00     7396     2252
#weight[6]       0.05       0.04   0.03   0.03       0.01       0.12 1.00     7002     2412
#weight[7]       0.07       0.06   0.05   0.04       0.01       0.16 1.01     6784     2471
#weight[8]       0.05       0.04   0.03   0.03       0.01       0.11 1.00     7885     1953
#weight[9]       0.06       0.05   0.04   0.03       0.01       0.13 1.00     6779     2644
# showing 10 of 106954 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
#weight_mpnp_1.1$summary()  # exceeds vector memory on my laptop
output1 <- read_cmdstan_csv(weight_mpnp_1.1$output_files()[1])
output2 <- read_cmdstan_csv(weight_mpnp_1.1$output_files()[2])
output3 <- read_cmdstan_csv(weight_mpnp_1.1$output_files()[3])
output4 <- read_cmdstan_csv(weight_mpnp_1.1$output_files()[4])

draws1_mpnp1.1 <- as.data.frame(output1$post_warmup_draws)
draws2_mpnp1.1 <- as.data.frame(output2$post_warmup_draws)
draws3_mpnp1.1 <- as.data.frame(output3$post_warmup_draws)
draws4_mpnp1.1 <- as.data.frame(output4$post_warmup_draws)
draws_mpnp1.1 <- rbind(draws1_mpnp1.1, draws2_mpnp1.1, draws3_mpnp1.1, draws4_mpnp1.1)
rm(output1, output2, output3, output4)

colnames(draws_mpnp1.1)[2:106954] <- counts_df$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:106954, size = 30, replace = F)

# tidy data -- vector memory usually exhausted and won't run
tidy_draws_1.1 <- pivot_longer(draws_mpnp1.1[,2:106954], cols = everything(), names_to = 'dyad', values_to = 'draw')
tidy_draws_1.1$chain <- rep(1:4, each = 106953000)
tidy_draws_1.1$index <- rep(rep(1:1000, each = 106953),4)
head(tidy_draws_1.1, 10)
tail(tidy_draws_1.1, 10)

# subset for Sierra (F52) -- family = U17, F60, U21 and F98
tidy_sierra <- draws_mpnp1.1[,c('F52_M40','F52_M26','F52_F8','F52_M15','F52_F98','F52_U17','F52_U21','F52_F60')]
tidy_sierra <- pivot_longer(tidy_sierra, cols = everything(), names_to = 'dyad', values_to = 'draw')
tidy_sierra <- separate(tidy_sierra, dyad, sep = '_', remove = F, into = c('F52','id'))
tidy_sierra$label <- ifelse(tidy_sierra$id == 'U17', 'U17 (calf)',
                            ifelse(tidy_sierra$id == 'U21', 'U21 (relative)',
                                   ifelse(tidy_sierra$id == 'F60', 'F60 (relative)',
                                          ifelse(tidy_sierra$id == 'F98', 'F98 (relative)', 
                                                 ifelse(tidy_sierra$id == 'F8', 'F8 (unrelated female)',
                                                        ifelse(tidy_sierra$id == 'M26', 'M26 (adult male)',
                                                               ifelse(tidy_sierra$id == 'M15', 'M15 (unrelated pubescent male)', 'M40 (adult male)')))))))
tidy_sierra$label <- factor(tidy_sierra$label, levels = c('U17 (calf)', 'F60 (relative)', 'U21 (relative)', 'F98 (relative)', 'M40 (adult male)','M26 (adult male)', 'F8 (unrelated female)','M15 (unrelated pubescent male)'))
tidy_sierra$chain <- rep(1:4, each = length(tidy_sierra$dyad)/4)
tidy_sierra$index <- rep(rep(1:1000, each = length(unique(tidy_sierra$dyad))),4)

### save data 
write_csv(draws_mpnp1.1, 'data_processed/mpnp_bayesian_edgedistributions_a1.b1_22.03.03.csv')
#write_csv(tidy_draws, 'data_processed/mpnp_bayesian_edgedistributions_tidy_22.03.03.csv')

################ 7) Summarise and plot edge weights ################
# draws_mpnp1.1 <- read_csv('data_processed/mpnp_bayesian_edgedistributions_a1.b1_22.03.03.csv')
### build traceplots ####
# random dyads -- many are very wide, high uncertainty
plot(draws_mpnp1.1[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_mpnp1.1[,plot_cols[2]], col = 'tan')
lines(draws_mpnp1.1[,plot_cols[3]], col = 'orange')
lines(draws_mpnp1.1[,plot_cols[4]], col = 'green')
lines(draws_mpnp1.1[,plot_cols[5]], col = 'chocolate')
lines(draws_mpnp1.1[,plot_cols[6]], col = 'blue')
lines(draws_mpnp1.1[,plot_cols[7]], col = 'red')
lines(draws_mpnp1.1[,plot_cols[8]], col = 'seagreen')
lines(draws_mpnp1.1[,plot_cols[9]], col = 'purple')
lines(draws_mpnp1.1[,plot_cols[10]],col = 'magenta')
lines(draws_mpnp1.1[,plot_cols[11]],col = 'black')
lines(draws_mpnp1.1[,plot_cols[12]], col = 'tan')
lines(draws_mpnp1.1[,plot_cols[13]], col = 'orange')
lines(draws_mpnp1.1[,plot_cols[14]], col = 'green')
lines(draws_mpnp1.1[,plot_cols[15]], col = 'chocolate')
lines(draws_mpnp1.1[,plot_cols[16]], col = 'blue')
lines(draws_mpnp1.1[,plot_cols[17]], col = 'red')
lines(draws_mpnp1.1[,plot_cols[18]], col = 'seagreen')
lines(draws_mpnp1.1[,plot_cols[19]], col = 'purple')
lines(draws_mpnp1.1[,plot_cols[20]],col = 'magenta')
lines(draws_mpnp1.1[,plot_cols[21]],col = 'black')
lines(draws_mpnp1.1[,plot_cols[22]], col = 'tan')
lines(draws_mpnp1.1[,plot_cols[23]], col = 'orange')
lines(draws_mpnp1.1[,plot_cols[24]], col = 'green')
lines(draws_mpnp1.1[,plot_cols[25]], col = 'chocolate')
lines(draws_mpnp1.1[,plot_cols[26]], col = 'blue')
lines(draws_mpnp1.1[,plot_cols[27]], col = 'red')
lines(draws_mpnp1.1[,plot_cols[28]], col = 'seagreen')
lines(draws_mpnp1.1[,plot_cols[29]], col = 'purple')
lines(draws_mpnp1.1[,plot_cols[30]],col = 'magenta')

# Sierra = F52, herd members = F60+U21+F98, calf = U17. This looks really good, and not too wide (all seen a good number of times).
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,4000))
lines(draws_mpnp1.1$F52_M40, col = 'black')      # non-herd member, adult male
lines(draws_mpnp1.1$F52_M15, col = 'tan')        # non-herd member, pubescent male
lines(draws_mpnp1.1$F52_M203,col = 'orange')     # non-herd member, adult male
lines(draws_mpnp1.1$F52_M26, col = 'green')      # non-herd member, calf
lines(draws_mpnp1.1$F52_F8,  col = 'chocolate')  # non-herd member, adult female
lines(draws_mpnp1.1$F52_U9,  col = 'blue')       # non-herd member, calf
lines(draws_mpnp1.1$F52_F98, col = 'red')        # herd member most frequently absent from sightings
lines(draws_mpnp1.1$F52_U17, col = 'purple')     # calf
lines(draws_mpnp1.1$F52_U21, col = 'seagreen')   # sister
lines(draws_mpnp1.1$F52_F60, col = 'magenta')    # sister's calf

# ggplot version if tidying worked: plot for Sierra (F52)
(traceplot_sierra <- ggplot(tidy_sierra[tidy_sierra$chain == 1,], aes(x = index, y = draw, colour = label))+
    geom_line()+
    scale_color_viridis_d()+
    theme_classic()+
    scale_x_continuous('MCMC chain position',  expand = c(0,0))+
    scale_y_continuous('edge weight estimate', expand = c(0.02,0))+
    labs(colour = 'Partner (relationship)'))
ggsave(path = 'outputs/mpnp_22.03.03/', filename = 'mpnp_sierra_traceplot.pdf',
       plot = traceplot_sierra,
       width = 30, height = 24, units = 'cm')

(traceplot_sierra2 <- ggplot(tidy_sierra, aes(x = index, y = draw, colour = as.factor(chain)))+
    geom_line()+
    scale_color_viridis_d()+
    theme_light()+
    facet_wrap(.~label)+
    scale_x_continuous('MCMC chain position',  expand = c(0,0))+
    scale_y_continuous('edge weight estimate', expand = c(0.02,0))+
    theme(legend.position = c(0.85, 0.15), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))+
    labs(colour = 'Chain ID'))
ggsave(path = 'outputs/mpnp_22.02.01/', filename = 'mpnp_sierra_traceplot.pdf',
       plot = traceplot_sierra,
       width = 30, height = 24, units = 'cm')

# check chain mixing
plot(draws1_mpnp1.1$`1.weight[1]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_mpnp1.1$`1.weight[1]`, col = 'green')
lines(draws3_mpnp1.1$`1.weight[1]`, col = 'blue')
lines(draws4_mpnp1.1$`1.weight[1]`, col = 'magenta')
plot(draws1_mpnp1.1$`1.weight[2]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_mpnp1.1$`1.weight[2]`, col = 'green')
lines(draws3_mpnp1.1$`1.weight[2]`, col = 'blue')
lines(draws4_mpnp1.1$`1.weight[2]`, col = 'magenta')
plot(draws1_mpnp1.1$`1.weight[3]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_mpnp1.1$`1.weight[3]`, col = 'green')
lines(draws3_mpnp1.1$`1.weight[3]`, col = 'blue')
lines(draws4_mpnp1.1$`1.weight[3]`, col = 'magenta')
which(counts_df$dyad == 'F52_U17') # 41192
plot(draws1_mpnp1.1$`1.weight[41192]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_mpnp1.1$`1.weight[41192]`, col = 'green')
lines(draws3_mpnp1.1$`1.weight[41192]`, col = 'blue')
lines(draws4_mpnp1.1$`1.weight[41192]`, col = 'magenta')
which(counts_df$dyad == 'F52_F60') # 40896
plot(draws1_mpnp1.1$`1.weight[40896]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_mpnp1.1$`1.weight[40896]`, col = 'green')
lines(draws3_mpnp1.1$`1.weight[40896]`, col = 'blue')
lines(draws4_mpnp1.1$`1.weight[40896]`, col = 'magenta')
which(counts_df$dyad == 'F52_U21') # 41197
plot(draws1_mpnp1.1$`1.weight[41197]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_mpnp1.1$`1.weight[41197]`, col = 'green')
lines(draws3_mpnp1.1$`1.weight[41197]`, col = 'blue')
lines(draws4_mpnp1.1$`1.weight[41197]`, col = 'magenta')
which(counts_df$dyad == 'F52_F98') # 40937
plot(draws1_mpnp1.1$`1.weight[40937]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_mpnp1.1$`1.weight[40937]`, col = 'green')
lines(draws3_mpnp1.1$`1.weight[40937]`, col = 'blue')
lines(draws4_mpnp1.1$`1.weight[40937]`, col = 'magenta')

### density plots ####
# random columns
dens(draws_mpnp1.1[,plot_cols[1]], ylim = c(0,50), xlim = c(0,1), las = 1)
dens(add = T, draws_mpnp1.1[,plot_cols[2]], col = 'tan')
dens(add = T, draws_mpnp1.1[,plot_cols[3]], col = 'orange')
dens(add = T, draws_mpnp1.1[,plot_cols[4]], col = 'green')
dens(add = T, draws_mpnp1.1[,plot_cols[5]], col = 'chocolate')
dens(add = T, draws_mpnp1.1[,plot_cols[6]], col = 'blue')
dens(add = T, draws_mpnp1.1[,plot_cols[7]], col = 'red')
dens(add = T, draws_mpnp1.1[,plot_cols[8]], col = 'seagreen')
dens(add = T, draws_mpnp1.1[,plot_cols[9]], col = 'purple')
dens(add = T, draws_mpnp1.1[,plot_cols[10]],col = 'magenta')
dens(add = T, draws_mpnp1.1[,plot_cols[11]],col = 'black')
dens(add = T, draws_mpnp1.1[,plot_cols[12]],col = 'tan')
dens(add = T, draws_mpnp1.1[,plot_cols[13]],col = 'orange')
dens(add = T, draws_mpnp1.1[,plot_cols[14]],col = 'green')
dens(add = T, draws_mpnp1.1[,plot_cols[15]],col = 'chocolate')
dens(add = T, draws_mpnp1.1[,plot_cols[16]],col = 'blue')
dens(add = T, draws_mpnp1.1[,plot_cols[17]],col = 'red')
dens(add = T, draws_mpnp1.1[,plot_cols[18]],col = 'seagreen')
dens(add = T, draws_mpnp1.1[,plot_cols[19]],col = 'purple')
dens(add = T, draws_mpnp1.1[,plot_cols[20]],col = 'magenta')
dens(add = T, draws_mpnp1.1[,plot_cols[21]],col = 'black')
dens(add = T, draws_mpnp1.1[,plot_cols[22]],col = 'tan')
dens(add = T, draws_mpnp1.1[,plot_cols[23]],col = 'orange')
dens(add = T, draws_mpnp1.1[,plot_cols[24]],col = 'green')
dens(add = T, draws_mpnp1.1[,plot_cols[25]],col = 'chocolate')
dens(add = T, draws_mpnp1.1[,plot_cols[26]],col = 'blue')
dens(add = T, draws_mpnp1.1[,plot_cols[27]],col = 'red')
dens(add = T, draws_mpnp1.1[,plot_cols[28]],col = 'seagreen')
dens(add = T, draws_mpnp1.1[,plot_cols[29]],col = 'purple')
dens(add = T, draws_mpnp1.1[,plot_cols[30]],col = 'magenta')
dens(add = T, draws_mpnp1.1[,which(counts_df$dyad == 'F52_U17')+1],col = 'magenta')
dens(add = T, draws_mpnp1.1[,which(counts_df$dyad == 'F52_F60')+1],col = 'purple')
dens(add = T, draws_mpnp1.1[,which(counts_df$dyad == 'F52_U21')+1],col = 'blue')
dens(add = T, draws_mpnp1.1[,which(counts_df$dyad == 'F52_F98')+1],col = 'red')

# density plots -- only run if tidying data worked
dens(tidy_sierra[tidy_sierra$id == 'M40',]$draw, xlab = 'edge weight estimation',       # all chains, plotted together
     ylim = c(0,60), xlim = c(0,1), las = 1)                                 # non-herd member, adult male
dens(tidy_sierra[tidy_sierra$id == 'M15',]$draw, add = T, col = 'tan')       # non-herd member, pubescent male
dens(tidy_sierra[tidy_sierra$id == 'M26',]$draw, add = T, col = 'green')     # non-herd member, calf
dens(tidy_sierra[tidy_sierra$id == 'F8',]$draw,  add = T, col = 'chocolate') # non-herd member, adult female
dens(tidy_sierra[tidy_sierra$id == 'F98',]$draw, add = T, col = 'red')       # herd member
dens(tidy_sierra[tidy_sierra$id == 'U17',]$draw, add = T, col = 'purple')    # calf
dens(tidy_sierra[tidy_sierra$id == 'U21',]$draw, add = T, col = 'seagreen')  # sister
dens(tidy_sierra[tidy_sierra$id == 'F60',]$draw, add = T, col = 'magenta')   # sister's calf
text('M40', x = 0.15, y = 05, col = 'black')
text('M15', x = 0.05, y = 55, col = 'tan')
text('M26', x = 0.05, y = 22, col = 'green')
text('F8' , x = 0.06, y = 18, col = 'chocolate')
text('F98', x = 0.60, y = 08, col = 'red')
text('U17', x = 0.95, y = 16, col = 'purple')
text('U21', x = 0.80, y = 10, col = 'seagreen')
text('F60', x = 0.86, y = 11, col = 'magenta')

(sierra_density <-
    ggplot(data = tidy_sierra[tidy_sierra$chain == 1], aes(x = draw, fill = ID))+       # only 1 chain
    geom_density()+  # why different values to the rethinking::dens() version?
    scale_fill_viridis_d(alpha = 0.4)+
    scale_x_continuous('edge weight estimation', expand = c(0,0), limits = c(0,1.1))+
    scale_y_continuous(expand = c(0,0), limits = c(0,55))+
    theme_classic()+
    theme(legend.position = c(0.5, 0.6)))
ggsave(path = 'outputs/mpnp_22.01.31/', filename = 'mpnp_sierra_density.pdf',
       plot = sierra_density,
       width = 20, height = 24, units = 'cm')

(sierra_density <-
    ggplot(data = tidy_sierra, aes(x = draw, col = as.factor(chain)))+                 # all chains, plotted separately
    geom_density()+ 
    scale_fill_viridis_d(alpha = 0.2)+
    scale_x_continuous('edge weight estimation', expand = c(0,0), limits = c(0,1.1))+
    scale_y_continuous(expand = c(0,0), limits = c(0,55))+
    theme_light()+
    facet_wrap(. ~ ID)+
    theme(legend.position = c(0.85, 0.15), legend.title = element_text(size = 14))+
    labs(colour = 'Chain ID'))
ggsave(path = 'outputs/mpnp_22.02.01/', filename = 'mpnp_sierra_density.pdf',
       plot = sierra_density,
       width = 20, height = 24, units = 'cm')

### summarise data ####
# summarise -- look for any anomalies in draw values or chain variation
summaries <- data.frame(dyad = colnames(draws_mpnp1.1[2:106954]),
                        min = rep(NA, ncol(draws_mpnp1.1)-1),
                        max = rep(NA, ncol(draws_mpnp1.1)-1),
                        mean = rep(NA, ncol(draws_mpnp1.1)-1),
                        median = rep(NA, ncol(draws_mpnp1.1)-1),
                        sd = rep(NA, ncol(draws_mpnp1.1)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_mpnp1.1[,i+1])
  summaries$max[i]    <- max(draws_mpnp1.1[,i+1])
  summaries$mean[i]   <- mean(draws_mpnp1.1[,i+1])
  summaries$median[i] <- median(draws_mpnp1.1[,i+1])
  summaries$sd[i]     <- sd(draws_mpnp1.1[,i+1])
}

summary(summaries)
summaries$dyad[which(summaries$min == min(summaries$min))] # F145_M68 = lowest value drawn
summaries$dyad[which(summaries$min == max(summaries$min))] # F17_U1 = highest lowest value drawn
summaries$dyad[which(summaries$max == min(summaries$max))] # M14_M56 = lowest highest value drawn
summaries$dyad[which(summaries$max == max(summaries$max))] # F44_U15 M131_M178 = highest value drawn (1)
summaries$dyad[which(summaries$sd == min(summaries$sd))]   # M14_M56
summaries$dyad[which(summaries$sd == max(summaries$sd))]   # F113_F117

plot(draws_mpnp1.1$F145_M68, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')       # still shows a fair amount of variation
lines(draws_mpnp1.1$F17_U1, col = 'red')                # still shows a fair amount of variation
lines(draws_mpnp1.1$M14_M56, col = 'blue')              # very low variation, but enough to trust the chain
lines(draws_mpnp1.1$F44_U15, col = 'green')             # still shows a fair amount of variation
lines(draws_mpnp1.1$M131_M178, col = 'orange')          # looks pretty crazy! -- 2 young males, each only seen once
summaries$min[which(summaries$dyad == 'M131_M178')]
summaries$sd[which(summaries$dyad == 'M131_M178')]
summaries$mean[which(summaries$dyad == 'M131_M178')]
lines(draws_mpnp1.1$F113_F117, col = 'purple')          # OK this one DOES look a bit crazy...

# organise dem_class
plot_data_mpnp1.1 <- left_join(x = summaries, y = counts_df, by = 'dyad')
head(plot_data_mpnp1.1)
plot_data_mpnp1.1$dem_class_1_cat <- as.integer(as.factor(plot_data_mpnp1.1$dem_class_1))
plot_data_mpnp1.1$dem_class_2_cat <- as.integer(as.factor(plot_data_mpnp1.1$dem_class_2))
plot_data_mpnp1.1$dem_type_short <- ifelse(plot_data_mpnp1.1$dem_class_1_cat <= plot_data_mpnp1.1$dem_class_2_cat,
                                            paste(plot_data_mpnp1.1$dem_class_1, plot_data_mpnp1.1$dem_class_2, sep = '_'),
                                            paste(plot_data_mpnp1.1$dem_class_2, plot_data_mpnp1.1$dem_class_1, sep = '_'))
sort(unique(plot_data_mpnp1.1$dem_type_short))
types <- data.frame(all = sort(unique(plot_data_mpnp1.1$dem_type_short)), type1 = NA, type2 = NA)
types <- separate(types, col = all, into = c('type1', 'type2'), remove = F)
types <- separate(types, type1, into = c('class1','sex1'), sep = 1, remove = F)
types <- separate(types, type2, into = c('class2','sex2'), sep = 1, remove = F)
head(types)
types$join <- types$all
types$join <- case_when(types$class1 == 'C' & types$class2 == 'C' ~ 'CC',
                        types$class1 == 'J' & types$class2 == 'C' ~ 'JC',
                        types$class1 == 'C' & types$class2 == 'J' ~ 'JC',
                        types$class1 == 'J' & types$class2 == 'J' ~ 'JJ',
                        types$class1 == 'C' & types$class2 == 'P' ~ paste(types$type2,'C',sep='_'),
                        types$class1 == 'P' & types$class2 == 'C' ~ paste(types$type1,'C',sep='_'),
                        types$class1 == 'C' & types$class2 == 'A' ~ paste(types$type2,'C',sep='_'),
                        types$class1 == 'A' & types$class2 == 'C' ~ paste(types$type1,'C',sep='_'),
                        types$class1 == 'J' & types$class2 == 'P' ~ paste(types$type2,'J',sep='_'),
                        types$class1 == 'P' & types$class2 == 'J' ~ paste(types$type1,'J',sep='_'),
                        types$class1 == 'J' & types$class2 == 'A' ~ paste(types$type2,'J',sep='_'),
                        types$class1 == 'A' & types$class2 == 'J' ~ paste(types$type1,'J',sep='_'),
                        TRUE ~ types$all)
unique(types$join)
types <- types[,c(1,8)]
colnames(types)[1] <- 'dem_type_short'
plot_data_mpnp1.1 <- left_join(plot_data_mpnp1.1, types, by = 'dem_type_short')
plot_data_mpnp1.1$join <- paste(plot_data_mpnp1.1$join, ' ', sep = ' ')
which(is.na(plot_data_mpnp1.1$join))

# boxplot edge weights by demographic type of dyad -- all types
colours <- c('magenta','purple','grey','grey','magenta','purple','grey','blue','purple','purple','purple','blue','purple','grey','grey','grey','grey','grey','magenta','purple','grey','purple','purple','blue','purple','grey','grey')
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_mpnp1.1, aes(y = mean, x = join))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 fill = colours)+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown',
         y = 'mean edge weight')+
    coord_flip())
ggsave(path = 'outputs/mpnp_22.03.03/', filename = 'mpnp_all_edge_vs_demtype.pdf',
       plot = edge_vs_demtype_all,
       width = 20, height = 24, units = 'cm')  

# boxplot edge weights by demographic type of dyad -- only adults/pubescents
adults_mpnp1.1 <- plot_data_mpnp1.1[plot_data_mpnp1.1$age_cat_id_1 > 2 & plot_data_mpnp1.1$age_cat_id_2 > 2,]
adults_mpnp1.1 <- adults_mpnp1.1[!is.na(adults_mpnp1.1$dyad),]
colours <- c('magenta','purple','magenta','purple','grey','blue','purple','blue','purple','magenta','purple','grey','blue','purple')
(edge_vs_demtype <- 
    ggplot(data = adults_mpnp1.1, aes(y = mean, x = join))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 fill = colours)+
    theme_classic()+
    labs(x = 'dyad demography (A/P = adult/pubescent, M/F/U = male/female/unknown',
         y = 'mean edge weight')+
    coord_flip())
ggsave(path = 'outputs/mpnp_22.03.03/', filename = 'mpnp_adults_edge_vs_demtype.pdf',
       plot = edge_vs_demtype,
       width = 20, height = 24, units = 'cm')  

# scatter plot of all adult edges by age difference
adults_mpnp1.1$sex_type_names <- ifelse(adults_mpnp1.1$sex_type == 'F_F', 'female-female',
                                         ifelse(adults_mpnp1.1$sex_type == 'F_M', 'male-female', 'male-male'))
(edge_vs_agediff <- 
    ggplot(adults_mpnp1.1[adults_mpnp1.1$sex_1 != 'U' & adults_mpnp1.1$sex_2 != 'U',],
           aes(x = age_diff, y = mean, fill = sex_type))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    facet_grid(~sex_type_names)+
    scale_x_continuous('age category difference between dyad members',
                       expand = c(0.02,0),
                       breaks = c(0:5))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))
ggsave(path = 'outputs/mpnp_22.03.03/', filename = 'mpnp_adults_edge_vs_agediff.pdf',
       plot = edge_vs_agediff,
       width = 30, height = 24, units = 'cm')

rm(draws1_mpnp1.1, draws2_mpnp1.1, draws3_mpnp1.1, draws4_mpnp1.1)
rm(edge_vs_agediff, edge_vs_demtype, edge_vs_demtype_all, plot_data_mpnp1.1, types, colours, adults_mpnp1.1)

################ 8) Create network plots ################
head(summaries)
length(unique(plot_data_mpnp1.1$id_1))+1 # number of individuals = 463

### Create igraph object for all elephants and then just plot certain individuals ####
# read in draws data -- very slow!
#draws_mpnp1.1 <- read_csv('data_processed/mpnp_bayesian_edgedistributions_a1.b1_22.03.03.csv')

# reassign dyad numbers to remove gaps
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df$id_1))+1,
                         NROW(unique(counts_df$id_2))+1,
                         NROW(draws_mpnp1.1)),
                    dimnames = list(c(unique(counts_df$id_1),'U9'),
                                    c('F1',unique(counts_df$id_2)),
                                    NULL))
N <- nrow(counts_df)

for (i in 1:54000) {            # can cope with jumps of 50000 dyads at a time
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1.1[, i+1]
}
for (i in 54001:N) {
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1.1[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- ((adj_upper - adj_lower)/adj_mid) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ele_nodes <- read_csv('data_processed/mpnp_elenodes_22.01.13.csv')
nodes <- data.frame(id = sort(unique(ele_nodes$id)))
nodes <- left_join(nodes, ele_nodes, by = 'id')
nodes$sex       <- as.factor(nodes$sex)
nodes$age_class <- as.factor(nodes$age_class)
nodes$dem_class <- as.factor(nodes$dem_class)
str(nodes)

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_mid)$weight,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = 'black',
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.size = 8,
     vertex.label = c(sort(unique(counts_df$id_1)),'U9'),
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(nodes$age_class == 'Adult','seagreen1',
                          ifelse(nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(nodes$age_class == 'Juvenile', 'yellow','magenta'))),
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
     vertex.size = 5,
     vertex.label = c(sort(unique(counts_df$id_1)),'U9'),
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(nodes$age_class == 'Adult','seagreen1',
                          ifelse(nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

# create variables for different degrees of node connectedness
nodes$degree_0.1 <- NA
nodes$degree_0.2 <- NA
nodes$degree_0.3 <- NA
nodes$degree_0.4 <- NA
nodes$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(nodes)){
  rows <- summaries[summaries$id_1 == nodes$id[i] | summaries$id_2 == nodes$id[i],]
  nodes$degree_0.1[i] <- length(which(rows$median > 0.1))
  nodes$degree_0.2[i] <- length(which(rows$median > 0.2))
  nodes$degree_0.3[i] <- length(which(rows$median > 0.3))
  nodes$degree_0.4[i] <- length(which(rows$median > 0.4))
  nodes$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(nodes$degree_0.1 < nodes$degree_0.2)
which(nodes$degree_0.2 < nodes$degree_0.3)
which(nodes$degree_0.3 < nodes$degree_0.4)
which(nodes$degree_0.4 < nodes$degree_0.5)

# plot network with reduced nodes -- only those with degree values > 0.5
g_mid_0.5 <- delete.vertices(graph = g_mid, v = nodes$id[which(nodes$degree_0.5 == 0)])
g_rng_0.5 <- delete.vertices(graph = g_rng, v = nodes$id[which(nodes$degree_0.5 == 0)])

coords_0.5 <- layout_nicely(g_mid_0.5)
plot(g_mid_0.5,
     edge.width = ifelse(E(g_mid_0.5)$weight > 0.5, E(g_rng_0.5)$weight, 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.5)
plot(g_mid_0.5,
     edge.width = ifelse(E(g_mid_0.5)$weight > 0.5, E(g_mid_0.5)$weight, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(nodes[which(nodes$degree_0.5 != 0),]$age_class == 'Adult', 'seagreen1',
                           ifelse(nodes[which(nodes$degree_0.5 != 0),]$age_class == 'Pubescent','skyblue',
                                  ifelse(nodes[which(nodes$degree_0.5 != 0),]$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords_0.5,
     add = TRUE)

#  plot network with reduced nodes -- only those with degree values > 0.3
g_mid_0.3 <- delete.vertices(graph = g_mid, v = nodes$id[which(nodes$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = nodes$id[which(nodes$degree_0.3 == 0)])

coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(nodes[which(nodes$degree_0.3 != 0),]$age_class == 'Adult',
                           'seagreen1',
                           ifelse(nodes[which(nodes$degree_0.3 != 0),]$age_class == 'Pubescent',
                                  'skyblue',
                                  ifelse(nodes[which(nodes$degree_0.3 != 0),]$age_class == 'Juvenile',
                                         'yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

#plot males
g_mid_m <- delete.vertices(graph = g_mid,
                           v = nodes$id[which(nodes$dem_class != 'AM' &
                                                nodes$dem_class != 'PM')])
g_rng_m <- delete.vertices(graph = g_rng,
                           v = nodes$id[which(nodes$dem_class != 'AM' &
                                                nodes$dem_class != 'PM')])
males <- nodes[nodes$dem_class == 'AM' | nodes$dem_class == 'PM',]
coords_m <- layout_nicely(g_mid_m)
plot(g_mid_m,
     edge.width = E(g_mid_m)$weight,
     edge.color = 'black',
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m)
plot(g_mid_m,
     edge.color = rgb(0,0,0,0.25),
     edge.width = E(g_rng_m)$weight,
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males$age_class == 'Adult','seagreen1', 'skyblue'),
     layout = coords_m, add = TRUE)

# plot males >0.3
g_mid_m0.3 <- delete.vertices(graph = g_mid,
                              v = nodes$id[which(nodes$dem_class != 'AM' &
                                                   nodes$dem_class != 'PM' |
                                                   nodes$degree_0.3 == 0)])
g_rng_m0.3 <- delete.vertices(graph = g_rng,
                              v = nodes$id[which(nodes$dem_class != 'AM' &
                                                   nodes$dem_class != 'PM' |
                                                   nodes$degree_0.3 == 0)])
males0.3 <- males[males$degree_0.3 > 0,]
coords_m0.3 <- layout_nicely(g_mid_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_rng_m0.3)$weight*3, 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_mid_m0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males0.3$age_class == 'Adult', 'seagreen1', 'skyblue'),
     layout = coords_m0.3, add = TRUE)

m38 <- summaries[summaries$id_1 == 'M38' | summaries$id_2 == 'M38',]
m38 <- m38[m38$median > 0.3,] # leaves on elephants which are >0.3 connected to females/calves
m87 <- summaries[summaries$id_1 == 'M87' | summaries$id_2 == 'M87',]
m87 <- m87[m87$median > 0.3,] # leaves on elephants which are >0.3 connected to females/calves

g_mid_m0.3 <- delete_vertices(graph = g_mid_m0.3, v = c('M240','M71','M99','M233',
                                                        'M158','M148','M87','M38',
                                                        'M128','M10','M111','M95',
                                                        'M7','M70','M17','M117',
                                                        'M221','M156','M196','M80',
                                                        'M114','M177'))
g_rng_m0.3 <- delete_vertices(graph = g_rng_m0.3, v = c('M240','M71','M99','M233',
                                                        'M158','M148','M87','M38',
                                                        'M128','M10','M111','M95',
                                                        'M7','M70','M17','M117',
                                                        'M221','M156','M196','M80',
                                                        'M114','M177'))
coords_m0.3 <- layout_nicely(g_mid_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_rng_m0.3)$weight*3, 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_mid_m0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males0.3$age_class == 'Adult', 'seagreen1', 'skyblue'),
     layout = coords_m0.3, add = TRUE)

### clear environment
rm(coords, coords_0.3, coords_0.5, coords_m, coords_m0.3, dyad_row, edges, ele, ele_nodes, elephants, g_mid, g_mid_0.3, g_mid_0.5, g_mid_m, g_mid_m0.3, g_rng, g_rng_0.3, g_rng_0.5, g_rng_m, g_rng_m0.3, m38, m87, males, males0.3, rows, traceplot_sierra, traceplot_sierra2, i)

################ 10) Nodal regression ################
# The final common type of network analysis we will cover here is nodal regression, where a regression is performed to analyse the relationship between a nodal network metric(such as centrality) and nodal traits (such as age and sex).  These analyses are usually used to assess how network position depends on various biological factors, but can also be used where network position is a predictor.  Since node metrics are derivative measures of the network, uncertainty in edge weights should ideally propagate through to node metrics, and on to coefficients in regression analyses, giving an accurate estimate of the total uncertainty in inferred parameters.  The core of the nodal regression is similar to dyadic regression, taking the form:
#X(n)i j∼Poisson(λ(n)i jD(n)i j)
#logÄλ(n)i jä=ωi j+L(n)i j
#Mi(ω)∼Normal(β0+β1xi,σ2)
# where β0 is the intercept parameter, β1 is the slope parameter, and σ is the standard deviation of the residuals. Mi(ω) denotes the node metric estimates when applied to the edge weights estimated in the top part of the model. Calculating node metrics within the model may present a practical challenge when using standard model fitting software, as many node metrics can not be described purely in terms of closed-form equations (Freeman, 1978).  In this case, splitting up the model along the dashed line becomes an important step, as it allows the edge weights to be estimated using a standard piece of software, and leaves the regression part of the model to be fit using a sampler that supports numeric estimates of network centralities. As part of the code we have made available, we have written a custom Metropolis-Hastings sampler that allows numeric estimates of node metrics to be used when fitting the model. Importantly, our sampler samples from the joint distribution of edge weights, rather than sampling edge weights independently. This ensures that the effect of any structure in the observation data (such as location effects) is maintained and propagated through to the node metric estimates, and subsequently to regression coefficient estimates.

### create array of draws per dyad (distributions) -- same as above
adj_tensor <- array(0, c(NROW(unique(counts_df$id_1))+1,
                         NROW(unique(counts_df$id_2))+1,
                         NROW(draws_mpnp1.1)),
                    dimnames = list(c(unique(counts_df$id_1),'U9'),
                                    c('F1',unique(counts_df$id_2)),
                                    NULL))
N <- nrow(counts_df)

for (i in 1:54000) {            # can cope with jumps of 50000 dyads at a time
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1.1[, i+1]
}
for (i in 54001:N) {
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_mpnp1.1[, i+1]
}
adj_tensor[,,1]

### create centrality matrix
centrality_matrix <- matrix(0, nrow = NROW(draws_mpnp1.1), ncol = length(unique(counts_df$id_1))+1)
for (i in 1:NROW(draws_mpnp1.1)) {
  g <- graph_from_adjacency_matrix(adj_tensor[, , i], mode="undirected", weighted=TRUE)
  centrality_matrix[i, ] <- strength(g)
}

colnames(centrality_matrix) <- c(sort(unique(counts_df$id_1)),unique(counts_df$id_2)[462])
head(centrality_matrix)
write_csv(as.data.frame(centrality_matrix), 'data_processed/mpnp_bayesian_nodalregression_centralitymatrix_22.03.04.csv')

centrality_quantiles <- t(apply(centrality_matrix, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
centrality_quantiles

cq <- data.frame(id = rownames(centrality_quantiles),
                 strength_lwr = centrality_quantiles[,1],
                 strength_mid = centrality_quantiles[,2],
                 strength_upr = centrality_quantiles[,3])
cq$id[which(cq$strength_mid == max(cq$strength_mid))]
nodes <- left_join(nodes, cq, by = 'id')
head(nodes)

boxplot(nodes$strength_mid ~ nodes$sex, notch = T)
boxplot(nodes$strength_mid ~ nodes$age_class, notch = T)
boxplot(nodes$strength_mid ~ nodes$age_category, notch = T, horizontal = T, las = 1, ylab = '')

nodes$age <- ifelse(nodes$age_category == '0-3', 1,
                    ifelse(nodes$age_category == '1-2', 1,
                           ifelse(nodes$age_category == '3-4', 1,
                                  ifelse(nodes$age_category == '4-5', 1,
                                         ifelse(nodes$age_category == '5-6', 2,
                                                ifelse(nodes$age_category == '6-7', 2,
                                                       ifelse(nodes$age_category == '7-8', 2,
                                                              ifelse(nodes$age_category == '8-9', 2,
                                                                     ifelse(nodes$age_category == '9-10', 3,
                                                                            ifelse(nodes$age_category == '10-15', 3,
                                                                                   ifelse(nodes$age_category == '15-19', 4,
                                                                                          ifelse(nodes$age_category == '20-25', 5,
                                                                                                 ifelse(nodes$age_category == '20-35', 5,
                                                                                                        ifelse(nodes$age_category == '25-40', 6,
                                                                                                               ifelse(nodes$age_category == '35-50', 6,
                                                                                                                      ifelse(nodes$age_category == '40+', 7,
                                                                                                                             ifelse(nodes$age_category == '50+', 7, 'error')))))))))))))))))
nodes$age <- ifelse(nodes$id == 'U8', 1, nodes$age)
boxplot(nodes$strength_mid ~ nodes$age, notch = T, horizontal = T, las = 1, ylab = '')

nodes_tidy <- pivot_longer(data = nodes, cols = 14:21,
                           names_to = 'centrality_type', values_to = 'centrality')
nodes_tidy$centrality_type <- factor(nodes_tidy$centrality_type,
                                     levels = c('strength_lwr', 'strength_mid', 'strength_upr',
                                                'degree_0.1','degree_0.2','degree_0.3',
                                                'degree_0.4','degree_0.5'))
ggplot(nodes_tidy, aes(y = centrality, x = age, fill = centrality_type))+
  geom_boxplot()+
  facet_wrap(.~centrality_type, scales = 'free')+ # lwr, mid and upr = very similar patterns. degree = very variable depending on level
  theme_light()+
  theme(legend.position = 'none',
        strip.text = element_text(colour = 'black'))

ggplot(nodes_tidy[nodes_tidy$centrality_type == 'degree_0.2' |
                    nodes_tidy$centrality_type == 'strength_mid',],
       aes(y = centrality, x = age))+
  geom_boxplot()+
  facet_wrap( . ~ centrality_type )

### write nodes data frame to file
write_csv(nodes, 'data_processed/mpnp_elenodes_centrality_22.03.03.csv')