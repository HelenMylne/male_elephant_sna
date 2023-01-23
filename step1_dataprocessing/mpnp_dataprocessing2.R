# EfA data model prep
#### Information ####
# Data collected by Elephants for Africa (EfA) 2012-2021
# Data supplied by Dr Kate Evans
#### Set up ####
library(tidyverse, lib.loc = 'packages/')   # library(tidyverse)
library(janitor, lib.loc = 'packages/')     # library(janitor)
library(lubridate, lib.loc = 'packages/')   # library(lubridate)
library(hms, lib.loc = 'packages/')         # library(hms)
library(readxl, lib.loc = 'packages/')      # library(readxl)
library(dplyr, lib.loc = 'packages/')       # library(dplyr)

library(cmdstanr, lib.loc = 'packages/')    # library(cmdstanr)
library(rethinking, lib.loc = 'packages/')  # library(rethinking)
library(igraph, lib.loc = 'packages/')      # library(igraph)

########## Create initial data ###########
#### read in processed data files ####
data_files <- lapply(Sys.glob("../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_pairwiseevents/*.csv"), read.csv)
all <- as.data.frame(data_files[1])
for(i in 2:length(data_files)){
  new_data <- as.data.frame(data_files[i])
  all <- rbind(all, new_data)
  rm(new_data)
}
length(unique(all$obs_id))          # 3765 -- this was the correct number as produced in the GBI matrix

# clean environment
#rm(data_files)

#### recreate how data were produced to match up obs_id to actual encounters ####
# sightings data
#s <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211214.xlsx')  # import data
s <- readxl::read_excel('../data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211214.xlsx')  # import data
str(s)
colnames(s)[c(1:23,57)] <- s[2,c(1:23,57)]                     # set column names to values in second row of dataframe
colnames(s)[24:56] <- c('CM','CF','CU','CM','CF','CU','JM','JF','JU','YPM','YPF','YPU','OPM','OPF','OPU',
                        'YAM','YAF','YAU','MAM','MAF','MAU','OAM','OAF','OAU','UM','UF','UU','SM','SF','SU',
                        'AM','AF','AU')                        # set count column names to individual categories
s <- s[3:nrow(s),]                                             # remove first two rows
s <- janitor::clean_names(s)                                   # clean up

# individual data
#efa <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')  # read in
efa <- readxl::read_excel('../data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')  # read in
efa$time_cat <- hour(efa$Time)                                 # identify time category throughout the day
efa <- separate(efa, Time, into = c('wrong_date','time'), sep = ' ')  # remove wrong date in time column
efa$time <- as_hms(efa$time)                                   # set true time
efa$Date <- as_date(efa$Date)                                  # set true date
efa <- efa[,c(1:5,7,33,8:15,18:27,31:32)]                      # remove unnecessary/incorrect columns
efa <- janitor::clean_names(efa)                               # clean up

efa_long <- data.frame(encounter = efa$elephant_sighting_id,   # create a new dataframe for one row per individual but with fewer columns and some extras for counts per sighting etc.
                       date = efa$date, time = efa$time,       # time of sighting
                       gps_s = NA, gps_e = NA,                 # location of sighting
                       total_elephants = NA,                   # total number of individuals together
                       total_id = NA, perc_id = NA,            # count number identified and calculate proportion of total
                       type = ifelse(efa$sex_id == 1, 'MO', 'BH/MX'),  # herd type -- BH = Breeding herd, MO = male only, MX = Mixed (same codes as used for MOTNP)
                       elephant = ifelse(efa$elephant_id == '-', NA, efa$elephant_id), # ID of elephant (dash = unidentified)
                       sex = ifelse(efa$sex_id == 1, 'M','F/U'),       # M = Male, F = Female, U = Unknown
                       age_range = efa$age_range_id)           # estimated age category at time of sighting
for(i in 1:length(efa_long$total_elephants)) { # count total elephants in sighting
  efa_long$total_elephants[i] <- length(which(efa$elephant_sighting_id == efa_long$encounter[i]))
}
efa_id <- efa_long[!is.na(efa_long$elephant),]
for(i in 1:length(efa_id$total_id)) {          # count total elephants identified per encounter
  efa_id$total_id[i] <- length(which(efa_id$encounter == efa_long$encounter[i]))
}
efa_id$perc_id <- 100*(efa_id$total_id/efa_id$total_elephants) # calculate proportion of elephants identified from total

#efa_gps <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long.csv')
efa_gps <- read_csv('../data_processed/mpnp_eles_long.csv')

length(unique(efa_long$encounter)) # 5162
length(unique(efa_id$encounter))   # 3736
length(unique(efa_gps$encounter))  # 3736

efa_gps$location <- paste(efa_gps$latitude, efa_gps$longitude, sep = '_')  # create combined location variable that joins both columns for unique place
efa_gps$longitude[which(efa_gps$longitude < 24)] <- NA  # anything below 24 degrees is too far west
efa_gps$latitude[which(efa_gps$latitude < -20)]  <- NA  # anything closer to zero than -20 is too far north

### create group-by-individual matrix
eles_asnipe <- efa_gps[,c(1,4,5,3,14)]                              # encounter, date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                    # convert to integer
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)      # start from 1, not 1st January 1970
eles_asnipe$Time <- hour(eles_asnipe$time)*60*60 + minute(eles_asnipe$time)*60 + second(eles_asnipe$time) # convert time values to seconds through day

eles_asnipe <- eles_asnipe[,c('encounter','Date','Time','elephant','location')] # select desired columns
colnames(eles_asnipe)[4] <- 'ID'                                    # rename
eles_asnipe$ID <- as.character(eles_asnipe$ID)                      # make ID a character variable
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')        # pad out date
eles_asnipe$group <- with(eles_asnipe, paste(Date, Time, location, sep = '_')) # create unique variable per sighting
eles_asnipe$group_id <- as.integer(as.factor(eles_asnipe$group))    # unique value per sighting
length(unique(eles_asnipe$encounter))                               # 3736
length(unique(eles_asnipe$group))                                   # 3765

eles_asnipe2 <- eles_asnipe[,c(4,8)]                                # group and ID
colnames(eles_asnipe2)[2] <- 'group'                                # rename
eles_asnipe2 <- data.table::setDT(eles_asnipe2)                     # convert to data table
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe2, group = 'group', id = 'ID') # 3765 rows -- ROW 1 IN gbi_matrix = obs_id 1 IN FINAL OUTPUT

obs_id <- rownames(gbi_matrix)   # create new variable for sighting ID
asnipe <- unique(sort(as.integer(as.factor(eles_asnipe$group)))) ; which(obs_id != asnipe) # check if any groups don't match -- both 1:3765 

## clean up
rm(efa, efa_gps, efa_id, efa_long, eles_asnipe2, gbi_matrix, s, asnipe, i, obs_id)

#### create new variable to match up obs_id to groups ####
### SO: Create a new variable using eles_asnipe group -- do a for loop that adds 1 to the value every time a group changes and retains it whenever it is the same -- this will create a variable which is unique for 3765 encounters but which goes in the same order as the obs_id
### use observations dataframe which is unique by GROUP (NOT ENCOUNTER ID NUMBER) and then add number in order of group
eles_asnipe3 <- eles_asnipe %>% arrange(desc(encounter)) %>% arrange(desc(Time)) %>% arrange(desc(Date)) # convert to an order where group is genuinely going in order
which(eles_asnipe$encounter != eles_asnipe3$encounter)  # check to see which ones have changed -- row 9 = day 3421 = now in descending order as with date
eles_asnipe3$in_order_date <- c('yes', rep(NA,nrow(eles_asnipe3)-1)) # first is in order of date, rest NA
eles_asnipe3$in_order_time <- c('yes', rep(NA,nrow(eles_asnipe3)-1)) # first is in order of time, rest NA
for(i in 2:nrow(eles_asnipe3)){ # go through each and see if it is in order relative to the previous
  eles_asnipe3$in_order_date[i] <- ifelse(eles_asnipe3$Date[i] > eles_asnipe3$Date[i-1], 'no', 'yes')
  eles_asnipe3$in_order_time[i] <- ifelse(eles_asnipe3$Date[i] == eles_asnipe3$Date[i-1], 
                                          ifelse(eles_asnipe3$Time[i] <= eles_asnipe3$Time[i-1], 'yes', 'no'),
                                          'yes')
}
table(eles_asnipe3$in_order_date) ; table(eles_asnipe3$in_order_time) # sightings now in reverse chronological order, both by date and time

(N <- length(unique(eles_asnipe3$group)))   # 3765
length(unique(eles_asnipe3$encounter)) - N  # -29 -- encounters were not reliably numbered throughout dataset
eles_asnipe3$group_id_new <- c(N, rep(NA, nrow(eles_asnipe3)-1)) # first one in data is last sighting (highest number)
for(i in 2:nrow(eles_asnipe3)){             # new value for each sighting (start at 3765, down to 1)
  eles_asnipe3$group_id_new[i] <- ifelse(eles_asnipe3$group[i] == eles_asnipe3$group[i-1],
                                         eles_asnipe3$group_id_new[i-1],
                                         eles_asnipe3$group_id_new[i-1]-1)
}
summary(eles_asnipe3$group_id_new)          # min = 1, max = 3765 -- this has worked

observations <- eles_asnipe3[,c('Date','Time','location','group_id_new')] %>% distinct() # 3765 observations

## SO... to identify which sightings in main data frame are actually represented in binomial data (all): merge 'observations' data with Binomial social interactions
rm(eles_asnipe, eles_asnipe3, N)

#### merge sightings information into dyad data #####
colnames(observations) <- c('date', 'time', 'location', 'obs_id') # rename columns
all <- left_join(x = all, y = observations, by = 'obs_id')        # merge sightings info per dyad
head(all,10) ; tail(all,10)                                       # check structure
rm(observations)

#### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
all$dyad <- paste(all$node_1, all$node_2, sep = '_')              # create unique value per dyad
all$dyad_id <- as.integer(as.factor(all$dyad))                    # convert to integer
all$location_id <- as.integer(as.factor(all$location))            # create integer value for GPS locations
all <- all %>% distinct()                                         # remove duplicates wherever there is an event
head(all)                                                         # check structure

max(all$date) - min(all$date) # 3430 days. Time period for ALERT data = 504 days --> split into 7 time windows (8 boundaries)
periods <- seq(from = min(all$date), to = max(all$date), length.out = 8) # create date sequence split into chunks of similar length to MOTNP dataset
periods[8] <- periods[8]+1 # set to be one higher than the final higher otherwise it takes the last date and creates a whole new period
events <- all[all$social_event == 1,] ; events <- events[!is.na(events$node_1),] # cuts down to 26086 pairs from 55503976 dyad pairs -- assigning to time window doesn't take forever
events$period <- NA          # create new variable for ID of time window
for(i in 1:nrow(events)){    # fill time window variable
  events$period[i] <- which(periods <= events$date[i])[length(which(periods <= events$date[i]))] # fill cell with last value in vector
  if(i %% 1000 == 0) {print(i)} # progress marker
}
range(events$obs_id)    # check output
table(events$period)    # check output

periods ; View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

# check elephants all match up
length(unique(all$node_1)) ; length(unique(all$node_2))         # check contain the same number of elephants
length(unique(events$node_1)) ; length(unique(events$node_2))   # check contain the same number of elephants -- don't need to
#rm(all)

#### convert to Binomial model data format ####
df_split <- events %>%
  group_by(node_1, node_2, period) %>%          # aggregate all sightings of each dyad
  summarise(event_count = sum(social_event),    # convert to a count of sightings
            dyad_id = cur_group_id())
head(df_split)                                  # check output
rm(events)

### add ID numbers
#eles <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long.csv') %>% 
#  select(elephant) %>%                                              # take only elephant ID column
#  distinct()                                                        # select unique values only
eles <- read_csv('../data_processed/mpnp_eles_long.csv') %>% 
  select(elephant) %>%                                              # take only elephant ID column
  distinct()                                                        # select unique values only
eles$node_1 <- as.integer(as.factor(eles$elephant))                 # create integer variable
colnames(eles)[1] <- c('id_1')                                      # rename variable
df <- left_join(df_split, eles, by = 'node_1') %>% distinct()       # add node ID to data frame about elephant 1
colnames(eles) <- c('id_2','node_2')                                # rename variables
df <- left_join(df, eles, by = 'node_2') %>% distinct()             # add node ID to data frame about elephant 2
head(df) ; tail(df)                                                 # check output
table(df$event_count)                                               # seems sensible
rm(df_split)

df$dyad <- paste(df$id_1, df$id_2, sep = '_')  # unique value for every dyad
df$dyad_id_period <- df$dyad_id                # every dyad has it's own ID number, including if same dyad in a different time window
df$dyad_id <- as.integer(as.factor(df$dyad))   # every dyad has it's own ID number, but same dyad in different windows share ID number

#### create dyad row for all pairs per period ####
dyads <- data.frame(id_1 = rep(sort(eles$id_2), each = nrow(eles)),   # new data frame containing 1 row per dyad (duplicated so id1_id2 and id2_id1 are both present, as are dyads where id1 = id2)
                    id_2 = rep(sort(eles$id_2), nrow(eles)))

colnames(eles) <- c('id_1','node_1')            # rename variables for merging
dyads <- left_join(dyads, eles, by = 'id_1')    # merge elephant information into dyad data frame
colnames(eles) <- c('id_2','node_2')            # rename variables for merging
dyads <- left_join(dyads, eles, by = 'id_2')    # merge elephant information into dyad data frame
dyads <- dyads[dyads$node_1 < dyads$node_2,]    # remove any rows where id1 = id2 and duplicates (id2_id1)
rm(all, eles, i)

dyads <- data.frame(id_1 = rep(dyads$id_1, length(unique(df$period))),        # reform data frame to the structure desired
                    id_2 = rep(dyads$id_2, length(unique(df$period))),
                    node_1 = rep(dyads$node_1, length(unique(df$period))),
                    node_2 = rep(dyads$node_2, length(unique(df$period))),
                    period = rep(sort(unique(df$period)), each = nrow(dyads)))

head(df) ; head(dyads)                                                # check structures of data frames
data <- left_join(x = dyads, y = df, by = c('id_1','id_2','period'))  # merge by dyad per time window
head(data)                                                            # check structure
data <- data[,c(1:5,8:10)]                                            # select desired variables
colnames(data)[3:4] <- c('node_1','node_2')                           # rename variables
data$event_count <- ifelse(is.na(data$event_count) == TRUE, 0, data$event_count)   # set NA counts (dyad not observed) to 0
table(data$event_count)                                               # check numbers of times dyads were observed
data$dyad <- paste(data$id_1, data$id_2, sep = '_')                   # create character variable with unique value per dyad
data$dyad_id <- as.integer(as.factor(data$dyad))                      # create numeric variable with unique value per dyad
head(data, 20)                                                        # check structure

efa <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx') %>% janitor::clean_names()  # read in data
efa$date <- as_date(efa$date)  # obtain date values

windows <- data.frame(period_start = seq(from = min(efa$date), to = max(efa$date), length.out = 8)[1:7],
                      period = 1:7,
                      period_days = periods[1:7])  # identify time windows

## check all time windows same length
windows$period_start[1]+windows$period_days[2]-windows$period_days[1]
windows$period_start[2]+windows$period_days[3]-windows$period_days[2]
windows$period_start[3]+windows$period_days[4]-windows$period_days[3]
windows$period_start[4]+windows$period_days[5]-windows$period_days[4]
windows$period_start[5]+windows$period_days[6]-windows$period_days[5]
windows$period_start[6]+windows$period_days[7]-windows$period_days[6]

data <- left_join(x = data, y = windows, by = 'period')  # merge time window ID with event data
head(data)

#### write to file ####
write_delim(data, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_pairwiseevents_22.05.30.csv', delim = ',')

## clean environment
rm(df, dyads, efa)

########## add additional information and remove impossible sightings ##########
#data <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_bayesian_pairwiseevents_22.05.30.csv')

#### break down into periods so easier to manipulate ####
table(data$period)
mpnp1 <- data[data$period == 1,]  # 1st time window
mpnp2 <- data[data$period == 2,]  # 2nd time window
mpnp3 <- data[data$period == 3,]  # 3rd time window
mpnp4 <- data[data$period == 4,]  # 4th time window
mpnp5 <- data[data$period == 5,]  # 5th time window
mpnp6 <- data[data$period == 6,]  # 6th time window
mpnp7 <- data[data$period == 7,]  # 6th time window
rm(data)

#### add count data per elephant and remove sightings where elephant not seen at all in period ####
# insert column for sighting period
sightings <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long.csv')
sightings <- sightings[,c(3,4)]                            # select only elephant ID and date columns
sightings$day <- as.numeric(sightings$date)                # numeric date
sightings$day <- sightings$day - (min(sightings$day)-1)    # convert date to start at 1 on first day of study
summary(sightings$day)

sightings$period <- NA                                     # create new variable for the time window
for(i in 1:nrow(sightings)){
  sightings$period[i] <- which(periods <= sightings$day[i])[length(which(periods <= sightings$day[i]))] # set time window value as last value in vector
}
periods ; sightings[sample(1:nrow(sightings),20),]         # visual check that individual dates match time windows

# create data frame of counts per dyad per time window
counts <- data.frame(id = rep(unique(sightings$elephant), 7),
                     period = rep(1:7, each = length(unique(sightings$elephant))),
                     count_all = NA,
                     count_period = NA)
for(i in 1:nrow(counts)){
  individual <- sightings[sightings$elephant == counts$id[i],]    # create data frame of all sightings per individual
  counts$count_all[i] <- nrow(individual)                         # count all sightings = number of rows in produced data frame
  counts$count_period[i] <- nrow(individual[individual$period == counts$period[i],])  # count sighitngs within time window = number of rows where time window matches
  if(i %% 1000 == 0) {print(i)}
}
rm(individual, i)

summary(counts$count_all)      # check output
summary(counts$count_period)   # check output

# join counts data frame with dyad data
colnames(counts)[1] <- 'id_1'                                                   # rename ID column for merging
data1 <- left_join(x = mpnp1, y = counts, by = c('id_1','period')) ; rm(mpnp1)  # add id1 count data for 1st time window
data2 <- left_join(x = mpnp2, y = counts, by = c('id_1','period')) ; rm(mpnp2)  # add id1 count data for 2nd time window
data3 <- left_join(x = mpnp3, y = counts, by = c('id_1','period')) ; rm(mpnp3)  # add id1 count data for 3rd time window
data4 <- left_join(x = mpnp4, y = counts, by = c('id_1','period')) ; rm(mpnp4)  # add id1 count data for 4th time window
data5 <- left_join(x = mpnp5, y = counts, by = c('id_1','period')) ; rm(mpnp5)  # add id1 count data for 5th time window
data6 <- left_join(x = mpnp6, y = counts, by = c('id_1','period')) ; rm(mpnp6)  # add id1 count data for 6th time window
data7 <- left_join(x = mpnp7, y = counts, by = c('id_1','period')) ; rm(mpnp7)  # add id1 count data for 7th time window

data1 <- data1[data1$count_period > 0,]   # remove all dyads id1 was never seen in 1st time window
data2 <- data2[data2$count_period > 0,]   # remove all dyads id1 was never seen in 2nd time window
data3 <- data3[data3$count_period > 0,]   # remove all dyads id1 was never seen in 3rd time window
data4 <- data4[data4$count_period > 0,]   # remove all dyads id1 was never seen in 4th time window
data5 <- data5[data5$count_period > 0,]   # remove all dyads id1 was never seen in 5th time window
data6 <- data6[data6$count_period > 0,]   # remove all dyads id1 was never seen in 6th time window
data7 <- data7[data7$count_period > 0,]   # remove all dyads id1 was never seen in 7th time window

colnames(data1)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data2)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data3)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data4)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data5)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data6)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data7)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging

colnames(counts)[1] <- 'id_2'                                                   # rename ID column for merging
data1 <- left_join(x = data1, y = counts, by = c('id_2','period'))              # add id2 count data for 1st time window
data2 <- left_join(x = data2, y = counts, by = c('id_2','period'))              # add id2 count data for 2nd time window
data3 <- left_join(x = data3, y = counts, by = c('id_2','period'))              # add id2 count data for 3rd time window
data4 <- left_join(x = data4, y = counts, by = c('id_2','period'))              # add id2 count data for 4th time window
data5 <- left_join(x = data5, y = counts, by = c('id_2','period'))              # add id2 count data for 5th time window
data6 <- left_join(x = data6, y = counts, by = c('id_2','period'))              # add id2 count data for 6th time window
data7 <- left_join(x = data7, y = counts, by = c('id_2','period'))              # add id2 count data for 7th time window

data1 <- data1[data1$count_period > 0,]   # remove all dyads id2 was never seen in 1st time window
data2 <- data2[data2$count_period > 0,]   # remove all dyads id2 was never seen in 2nd time window
data3 <- data3[data3$count_period > 0,]   # remove all dyads id2 was never seen in 3rd time window
data4 <- data4[data4$count_period > 0,]   # remove all dyads id2 was never seen in 4th time window
data5 <- data5[data5$count_period > 0,]   # remove all dyads id2 was never seen in 5th time window
data6 <- data6[data6$count_period > 0,]   # remove all dyads id2 was never seen in 6th time window
data7 <- data7[data7$count_period > 0,]   # remove all dyads id2 was never seen in 7th time window

colnames(data1)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data2)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data3)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data4)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data5)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data6)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data7)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging

rm(counts,sightings)

table(data1$count_period_1)  # check id1 counts for 1st time window
table(data1$count_period_2)  # check id2 counts for 1st time window

table(data2$count_period_1)  # check id1 counts for 2nd time window
table(data2$count_period_2)  # check id2 counts for 2nd time window

table(data3$count_period_1)  # check id1 counts for 3rd time window
table(data3$count_period_2)  # check id2 counts for 3rd time window

table(data4$count_period_1)  # check id1 counts for 4th time window
table(data4$count_period_2)  # check id2 counts for 4th time window

table(data5$count_period_1)  # check id1 counts for 5th time window
table(data5$count_period_2)  # check id2 counts for 5th time window

table(data6$count_period_1)  # check id1 counts for 6th time window
table(data6$count_period_2)  # check id2 counts for 6th time window

table(data7$count_period_1)  # check id1 counts for 6th time window
table(data7$count_period_2)  # check id2 counts for 6th time window

## check ID types -- B numbers are identified, T numbers are yet to be identified, F number is a female
unique(data1$id_1) # 695 B numbers,    1 F number
unique(data2$id_1) # 
unique(data3$id_1) # 
unique(data4$id_1) # 
unique(data5$id_1) # 
unique(data6$id_1) # 
unique(data7$id_1) # 

## remove unidentified elephants
data1_id <- data1 %>%
  separate(id_1, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num) %>%                                               # remove split variables
  separate(id_2, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num)                                                   # remove split variables
data2_id <- data2 %>%
  separate(id_1, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num) %>%                                               # remove split variables
  separate(id_2, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num)                                                   # remove split variables
data3_id <- data3 %>%
  separate(id_1, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number 
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num) %>%                                               # remove split variables
  separate(id_2, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num)                                                   # remove split variables
data4_id <- data4 %>%
  separate(id_1, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number 
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num) %>%                                               # remove split variables
  separate(id_2, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num)                                                   # remove split variables
data5_id <- data5 %>%
  separate(id_1, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number 
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num) %>%                                               # remove split variables
  separate(id_2, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num)                                                   # remove split variables
data6_id <- data6 %>%
  separate(id_1, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number 
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num) %>%                                               # remove split variables
  separate(id_2, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num)                                                   # remove split variables
data7_id <- data7 %>%
  separate(id_1, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number 
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num) %>%                                               # remove split variables
  separate(id_2, into = c('BT','num'), sep = 1, remove = FALSE) %>%   # split ID1 number into B/T/F and number
  filter(BT == 'B') %>%                                               # select only identified (B) numbers
  select(-BT, -num)                                                   # remove split variables

#rm(data1, data2, data3, data4, data5, data6)

#### add data about individual elephant ages ####
# read in data
eles <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long.csv') %>% 
  select(elephant, sex, age_range, date) %>%   # read in age data, remove unnecessary variables
  distinct()                                   # select unique rows (NOTE -- this is NOT just 1 row per elephant as some have been estimated to be in multiple age categories within a time window)

# identify period of each sighting
eles$period <- NA                                      # create new variable for time window per elephant sighting
for(i in 1:nrow(eles)){
  eles$period[i] <- which(windows$period_start <= eles$date[i])[length(which(windows$period_start <= eles$date[i]))] # set time window ID as last value in vector
}
windows$period_start ; eles[sample(1:nrow(eles),20),]  # visual check of window against date

eles <- eles %>% separate(elephant, into = c('BT','id_number'), sep = 1, remove = FALSE) %>% # split ID to B/T/F and number
  filter(BT == 'B') %>%                                                                      # remove unidentified elephants
  select(-date, -BT, -id_number) %>%                                                         # remove unnecessary columns
  distinct()                                                                                 # remove duplicates
unique(eles$elephant)                                                                        # check ID numbers

# identify elephants whose age category changes within time window
table(eles$age_range) # 0x1 (sub 1 year old), 0x9 (females >20), 120x10 (unknown)
eles$age_range_NA <- ifelse(eles$age_range == 10, NA, eles$age_range)  # set unknown ages to NA

eles$age_unsure <- NA    # new variable: are there multiple estimated age categories
eles$age_min <- NA       # new variable: minimum estimated age category
eles$age_max <- NA       # new variable: maximum estimated age category
eles$age_maxmin <- NA    # new variable: number of estimated age categories between max and min
eles$age_median <- NA    # new variable: median of estimated age categories
for (i in 1:nrow(eles)) {
  individual <- eles[eles$elephant == eles$elephant[i] & eles$period == eles$period[i],] # data frame for every sighting of elephant in that time window in which a different age estimate was made
  summary <- summary(individual$age_range_NA, na.rm = T)   # summary values for elephant age estimates in time window
  eles$age_unsure[i] <- nrow(individual)                   # number of sightings with different age estimates
  eles$age_min[i]    <- summary[1]                         # minimum age estimate
  eles$age_max[i]    <- summary[6]                         # maximum age estimate
  eles$age_maxmin[i] <- summary[6] - summary[1]            # range of age estimates
  eles$age_median[i] <- summary[3]                         # number of sightings with different age estimates
}
summary(eles$age_median)  # select median estimate

eles <- eles %>% select(-sex, -age_range, -age_range_NA, -age_unsure) %>% distinct()  # remove unnecessary variables

# add age data
colnames(eles)[c(1,3:6)] <- c('id_1','age_min_1','age_max_1','age_range_1','age_median_1') # rename variables for merging
data1_id <- left_join(x = data1_id, y = eles, by = c('id_1','period'))  # add ID1 age estimates for 1st time window
data2_id <- left_join(x = data2_id, y = eles, by = c('id_1','period'))  # add ID1 age estimates for 2nd time window
data3_id <- left_join(x = data3_id, y = eles, by = c('id_1','period'))  # add ID1 age estimates for 3rd time window
data4_id <- left_join(x = data4_id, y = eles, by = c('id_1','period'))  # add ID1 age estimates for 4th time window
data5_id <- left_join(x = data5_id, y = eles, by = c('id_1','period'))  # add ID1 age estimates for 5th time window
data6_id <- left_join(x = data6_id, y = eles, by = c('id_1','period'))  # add ID1 age estimates for 6th time window
data7_id <- left_join(x = data7_id, y = eles, by = c('id_1','period'))  # add ID1 age estimates for 7th time window

colnames(eles)[c(1,3:6)] <- c('id_2','age_min_2','age_max_2','age_range_2','age_median_2') # rename variables for merging
data1_id <- left_join(x = data1_id, y = eles, by = c('id_2','period'))  # add ID2 age estimates for 1st time window
data2_id <- left_join(x = data2_id, y = eles, by = c('id_2','period'))  # add ID2 age estimates for 2nd time window
data3_id <- left_join(x = data3_id, y = eles, by = c('id_2','period'))  # add ID2 age estimates for 3rd time window
data4_id <- left_join(x = data4_id, y = eles, by = c('id_2','period'))  # add ID2 age estimates for 4th time window
data5_id <- left_join(x = data5_id, y = eles, by = c('id_2','period'))  # add ID2 age estimates for 5th time window
data6_id <- left_join(x = data6_id, y = eles, by = c('id_2','period'))  # add ID2 age estimates for 6th time window
data7_id <- left_join(x = data7_id, y = eles, by = c('id_2','period'))  # add ID2 age estimates for 7th time window

rm(eles, individual)

### add variable for sightings apart ####
data1_id$count_dyad <- (data1_id$count_period_1 + data1_id$count_period_2) - data1_id$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data1_id$apart <- data1_id$count_dyad - data1_id$event_count               # seen apart = total sightings - seen together
colnames(data1_id)[6] <- 'together'                                        # rename variable
table(data1_id$together)                                                   # check values
table(data1_id$apart)                                                      # check values

data2_id$count_dyad <- (data2_id$count_period_1 + data2_id$count_period_2) - data2_id$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data2_id$apart <- data2_id$count_dyad - data2_id$event_count               # seen apart = total sightings - seen together
colnames(data2_id)[6] <- 'together'                                        # rename variable
table(data2_id$together)                                                   # check values
table(data2_id$apart)                                                      # check values

data3_id$count_dyad <- (data3_id$count_period_1 + data3_id$count_period_2) - data3_id$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data3_id$apart <- data3_id$count_dyad - data3_id$event_count               # seen apart = total sightings - seen together
colnames(data3_id)[6] <- 'together'                                        # rename variable
table(data3_id$together)                                                   # check values
table(data3_id$apart)                                                      # check values

data4_id$count_dyad <- (data4_id$count_period_1 + data4_id$count_period_2) - data4_id$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data4_id$apart <- data4_id$count_dyad - data4_id$event_count               # seen apart = total sightings - seen together
colnames(data4_id)[6] <- 'together'                                        # rename variable
table(data4_id$together)                                                   # check values
table(data4_id$apart)                                                      # check values

data5_id$count_dyad <- (data5_id$count_period_1 + data5_id$count_period_2) - data5_id$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data5_id$apart <- data5_id$count_dyad - data5_id$event_count               # seen apart = total sightings - seen together
colnames(data5_id)[6] <- 'together'                                        # rename variable
table(data5_id$together)                                                   # check values
table(data5_id$apart)                                                      # check values

data6_id$count_dyad <- (data6_id$count_period_1 + data6_id$count_period_2) - data6_id$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data6_id$apart <- data6_id$count_dyad - data6_id$event_count               # seen apart = total sightings - seen together
colnames(data6_id)[6] <- 'together'                                        # rename variable
table(data6_id$together)                                                   # check values
table(data6_id$apart)                                                      # check values

data7_id$count_dyad <- (data7_id$count_period_1 + data7_id$count_period_2) - data7_id$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data7_id$apart <- data7_id$count_dyad - data7_id$event_count               # seen apart = total sightings - seen together
colnames(data7_id)[7] <- 'together'                                        # rename variable
table(data7_id$together)                                                   # check values
table(data7_id$apart)                                                      # check values

## save data #####
data1_id <- data1_id[,c(8,7,1:5,9,6,24,12,14,23,15,19,16,20,17,21,18,22)]  # select only necessary variables for 1st time window
data2_id <- data2_id[,c(8,7,1:5,9,6,24,12,14,23,15,19,16,20,17,21,18,22)]  # select only necessary variables for 2nd time window
data3_id <- data3_id[,c(8,7,1:5,9,6,24,12,14,23,15,19,16,20,17,21,18,22)]  # select only necessary variables for 3rd time window
data4_id <- data4_id[,c(8,7,1:5,9,6,24,12,14,23,15,19,16,20,17,21,18,22)]  # select only necessary variables for 4th time window
data5_id <- data5_id[,c(8,7,1:5,9,6,24,12,14,23,15,19,16,20,17,21,18,22)]  # select only necessary variables for 5th time window
data6_id <- data6_id[,c(8,7,1:5,9,6,24,12,14,23,15,19,16,20,17,21,18,22)]  # select only necessary variables for 6th time window
data7_id <- data7_id[,c(8,7,1:5,9,6,24,12,14,23,15,19,16,20,17,21,18,22)]  # select only necessary variables for 6th time window

# write data files
write_csv(data1_id,'../data_processed/mpnp_period1_pairwiseevents.csv')
write_csv(data2_id,'../data_processed/mpnp_period2_pairwiseevents.csv')
write_csv(data3_id,'../data_processed/mpnp_period3_pairwiseevents.csv')
write_csv(data4_id,'../data_processed/mpnp_period4_pairwiseevents.csv')
write_csv(data5_id,'../data_processed/mpnp_period5_pairwiseevents.csv')

length(unique(data6_id$id_1))   # not usable data -- too few identified individuals
#write_csv(data6_id,'../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period6_pairwiseevents.csv')

length(unique(data7_id$id_1))   # not usable data -- too few identified individuals
#write_csv(data7_id,'../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period7_pairwiseevents.csv')
