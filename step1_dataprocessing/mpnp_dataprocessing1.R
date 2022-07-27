### Bayesian analysis of EfA data ####
# Script to process association data from Makgadikgadi Pans National Park, Botswana.
# Data collected: 1st November 2019 to 5th August 2021
# Collected by: Elephants for Africa (Dr Kate Evans)
# Data supplied by: EfA and Dr Kate Evans, 19th October 2021

### Set up ####
# load packages
library(tidyverse)
library(dplyr)
library(dagitty)
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
s <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211214.xlsx')
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
efa <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')
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

write_delim(efa_gps, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long_22.03.08.csv', delim = ',')
# efa_gps <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long_22.03.08.csv')

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

write_delim(ele_nodes, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_elenodes_22.03.08.csv')
# ele_nodes <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_elenodes_22.03.08.csv')

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
nrow(gbi_matrix) # 3765 -- run in 74 blocks of 50 sightings + 1 block of 65
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1.250_22.03.10.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 251:387) {
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
head(gbi_df)
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings251.387_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 388:400) {
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings388.400_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 401:450) {
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings401.450_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 451:500) {
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings451.500_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 501:550) {
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings501.550_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 551:600) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings551.600_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 601:650) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings601.650_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 651:700) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings651.700_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 701:750) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings701.750_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 751:800) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings751.800_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 801:850) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings801.850_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 851:900) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings851.900_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 901:950) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings901.950_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 951:1000) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings951.1000_22.03.08.csv', delim = ',')









gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1001:1050) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1001.1050_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1051:1100) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1051.1100_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1101:1150) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1101.1150_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1151:1200) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1151.1200_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1201:1250) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1201.1250_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1251:1300) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1251.1300_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1301:1350) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1301.1350_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1351:1400) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1351.1400_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1401:1450) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1401.1450_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1451:1500) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1451.1500_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1501:1550) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1501.1550_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1551:1600) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1551.1600_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1601:1650) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1601.1650_22.03.08.csv', delim = ',')




gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1651:1700) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1651.1700_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1701:1750) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1701.1750_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1751:1800) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1751.1800_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1801:1850) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1801.1850_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1851:1900) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1851.1900_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1901:1950) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1901.1950_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1951:2000) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings1951.2000_22.03.08.csv', delim = ',')







gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2001:2050) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2001.2050_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2051:2100) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2051.2100_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2101:2150) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2101.2150_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2151:2200) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2151.2200_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2201:2250) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2201.2250_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2251:2300) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2251.2300_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2301:2350) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2301.2350_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2351:2400) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2351.2400_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2401:2450) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2401.2450_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2451:2500) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2451.2500_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2501:2550) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2501.2550_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2551:2600) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2551.2600_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2601:2650) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2601.2650_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2651:2700) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2651.2700_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2701:2750) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2701.2750_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2751:2800) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2751.2800_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2801:2850) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2801.2850_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2851:2900) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2851.2900_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2901:2950) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2901.2950_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 2951:3000) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings2951.3000_22.03.08.csv', delim = ',')





gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3001:3050) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3001.3050_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3051:3100) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3051.3100_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3101:3150) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3101.3150_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3151:3200) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3151.3200_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3201:3250) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3201.3250_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3251:3300) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3251.3300_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3301:3350) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3301.3350_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3351:3400) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3351.3400_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3401:3450) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3401.3450_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3451:3500) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3451.3500_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3501:3550) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3501.3550_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3551:3600) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3551.3600_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3601:3650) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3601.3650_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3651:3700) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3651.3700_22.03.08.csv', delim = ',')

gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3701:3750) {
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
  if(obs_id %% 5 == 0) {print(obs_id)}
  if(obs_id %% 5 == 0) {print(Sys.time())}
}
gbi_df
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3701.3750_22.03.08.csv', delim = ',')


gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 3751:3765) {
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
write_delim(gbi_df, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_allpairwiseevents_sightings3751.3765_22.03.08.csv', delim = ',')
