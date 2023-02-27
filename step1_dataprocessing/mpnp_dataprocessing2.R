# EfA data model prep
#### Information ####
# Data collected by Elephants for Africa (EfA) 2012-2021
# Data supplied by Dr Kate Evans
#### Set up ####
library(tidyverse, lib.loc = '../packages/')   # library(tidyverse)
library(janitor, lib.loc = '../packages/')     # library(janitor)
library(hms, lib.loc = '../packages/')         # library(hms)
library(dplyr, lib.loc = '../packages/')       # library(dplyr)
library(lubridate, lib.loc = '../packages/')   # library(lubridate)
library(readxl, lib.loc = '../packages/')      # library(readxl)

#library(cmdstanr, lib.loc = '../packages/')    # library(cmdstanr)
#set_cmdstan_path('../packages/')
#library(rethinking, lib.loc = '../packages/')  # library(rethinking)
#library(igraph, lib.loc = '../packages/')      # library(igraph)

########## Create initial data ###########
#### read in processed data files ####
data_files <- lapply(Sys.glob("../data_processed/mpnp_pairwiseevents/*.RDS"), readRDS)
all <- as.data.frame(data_files[1])
for(i in 2:length(data_files)){
  new_data <- as.data.frame(data_files[i])
  all <- rbind(all, new_data)
  rm(new_data)
}
length(unique(all$obs_id))          # 1438 -- this was the correct number as produced in the GBI matrix

# clean environment
rm(data_files) ; gc()

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
#efa <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')        # read in
efa <- readxl::read_excel('../data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')  # read in
efa$time_cat <- lubridate::hour(efa$Time)                            # create variable which is just hour of the day
efa <- separate(efa, Time, into = c('wrong_date','time'), sep = ' ') # time column registers right time but wrong date
efa$time <- hms::as_hms(efa$time)                                    # extract time value only from time column
efa$Date <- as_date(efa$Date)                                        # set true date
efa <- efa[,c(1:5,7,33,8:32)]                                        # remove incorrect date column
efa$Age_Range_ID <- as.factor(efa$Age_Range_ID)                      # make factor not character
efa$Activity_ID <- as.factor(efa$Activity_ID)                        # make factor not character
efa$Distance_To_Observer <- as.numeric(efa$Distance_To_Observer)     # make number not character
efa$Physical_Condition_ID <- as.factor(efa$Physical_Condition_ID)    # make factor not character
efa <- efa[,c(1:15,18:27,31:32)]                                     # rearrange
efa$Date <- lubridate::as_date(efa$Date)                             # make date
efa <- janitor::clean_names(efa)                                     # clean up

# correct IDs
efa <- efa %>% separate(col = elephant_id, into = c('letter','number'), sep = 1, remove = F)
unique(efa$letter)
efa$elephant_id <- ifelse(efa$letter == 'b', paste0('B',efa$number), efa$elephant_id)

# add additional encounter numbers for those that are missing
efa$encounter <- paste(efa$elephant_sighting_id, efa$date, efa$time, sep = '_')
for(i in 1:nrow(efa)){
  if(is.na(efa$elephant_sighting_id[i]) == TRUE){
    x <- efa[efa$encounter == efa$encounter[i],]
    if(nrow(x) == 1){
      efa$elephant_sighting_id[i] <- max(efa$elephant_sighting_id, na.rm = TRUE)+1
    }
    else{
      if(length(which(is.na(x$elephant_sighting_id) == FALSE)) > 0 ){
        sighting_id <- x$elephant_sighting_id[which(is.na(x$elephant_sighting_id) == FALSE)[1]]
        efa$elephant_sighting_id[which(efa$encounter == efa$encounter[i])] <- sighting_id
      } else{
        sighting_id <- max(efa$elephant_sighting_id, na.rm = TRUE)+1
        efa$elephant_sighting_id[which(efa$encounter == efa$encounter[i])] <- sighting_id
      }
    }
  }
}

# trim down unnecessary columns and calculate proportions of identified individuals
efa_long <- data.frame(encounter = efa$elephant_sighting_id,   # create a new dataframe for one row per individual but with fewer columns and some extras for counts per sighting etc.
                       date = efa$date, time = efa$time,       # time of sighting
                       gps_s = NA, gps_e = NA,                 # location of sighting
                       total_elephants = NA,                   # total number of individuals together
                       total_id = NA, perc_id = NA,            # count number identified and calculate proportion of total
                       elephant = efa$elephant_id,             # ID of elephant (dash/T = unidentified)
                       letter = efa$letter,                    # so can filter out un-ID'd
                       age_range = efa$age_range_id)           # estimated age category at time of sighting

for(i in sort(unique(efa_long$encounter))) {                   # count total elephants in sighting
  encounter <- efa_long[efa_long$encounter == i,]
  efa_long$total_elephants[efa_long$encounter == i] <- nrow(encounter)                  # total eles
  efa_long$total_id[efa_long$encounter == i] <- length(which(encounter$letter == 'B'))  # total identified
}
rm(encounter, x, i) ; gc()
efa_long$perc_id <- 100*(efa_long$total_id/efa_long$total_elephants) # calculate proportion of elephants identified from total
summary(efa_long$perc_id)

# remove unidentified individuals
efa_long <- efa_long %>% filter(elephant != '-') %>% 
  filter(letter != 'T') %>% filter(letter != 'F')
colnames(efa_long)
efa_long <- efa_long[,c(1,9,2:8,11)]

# add gps data
colnames(s)[1] <- 'encounter'
s$encounter <- as.numeric(s$encounter)
efa_gps <- left_join(x = efa_long, y = s, by = 'encounter') # join ID dataframe to sightings data so that all information included in single dataframe. first 100-ish rows are blank past column 12 because 2 datasets formed at slightly different times -- group sightings only goes as far as 29th January 2021, individual sightings go until 24th September 2021
which(as.numeric(efa_gps$number_of_elephants) != efa_gps$total_elephants) # many -- go with total_elephants as based on input data
efa_gps <- efa_gps[,c(1:11,16:17)]       # remove unnecessary columns
colnames(efa_gps)[3] <- c('date')        # rename from "date.x"
efa_gps$location <- paste(efa_gps$latitude, efa_gps$longitude, sep = '_')  # create combined location variable that joins both columns for unique place
efa_gps$longitude[which(efa_gps$longitude < 24)] <- NA  # anything below 24 degrees is too far west
efa_gps$latitude[which(efa_gps$latitude < -20)]  <- NA  # anything closer to zero than -20 is too far north

length(unique(efa_long$encounter)) # 1407
length(unique(efa_gps$encounter))  # 1407

efa <- efa_gps
rm(efa_gps, efa_long, s, sighting_id) ; gc()

### create group-by-individual matrix
eles_asnipe <- efa[,c(3,4,2,14)]                                         # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                         # convert to integer
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)           # start from 1, not 1st January 1970
eles_asnipe$Time <- hour(eles_asnipe$time)*60*60 + minute(eles_asnipe$time)*60 + second(eles_asnipe$time) # convert time values to seconds through day
eles_asnipe <- eles_asnipe[,c(5,6,3,4)]                                  # select desired variables
colnames(eles_asnipe) <- c('Date','Time','ID','Location')                # rename for asnipe package
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 4, pad = '0')             # ensure dates go in the right order as character
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value per sighting
eles_asnipe$observation <- as.integer(as.factor(eles_asnipe$encounter))  # unique count value per encounter
max(eles_asnipe$observation)                                             # 1437 (3765 if un-ID'd included)
eles_asnipe_gbi <- eles_asnipe[,c(3,7)]                                  # ID and observation
eles_asnipe_gbi <- data.table::setDT(eles_asnipe_gbi)                    # convert to data table
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe_gbi, group = 'observation', id = 'ID') # create group-by-individual matrix

# check gbi_matrix
obs_id <- rownames(gbi_matrix)   # create new variable for sighting ID
asnipe <- unique(sort(as.integer(as.factor(eles_asnipe$group)))) ; which(obs_id != asnipe) # check if any groups don't match -- both 1:3765 

## clean up
#rm(efa, efa_gps, efa_id, efa_long, eles_asnipe2, gbi_matrix, s, asnipe, i, obs_id) ; gc()

#### create new variable to match up obs_id to groups ####
### SO: Create a new variable using eles_asnipe group -- do a for loop that adds 1 to the value every time a group changes and retains it whenever it is the same -- this will create a variable which is unique for 3765 encounters but which goes in the same order as the obs_id
### use observations dataframe which is unique by GROUP (NOT ENCOUNTER ID NUMBER) and then add number in order of group
eles_asnipe2 <- eles_asnipe %>% arrange(desc(observation)) %>% arrange(desc(Time)) %>% arrange(desc(Date)) # convert to an order where group is genuinely going in order
which(eles_asnipe$encounter != eles_asnipe2$encounter)  # check to see which ones have changed -- row 3 = day 2856 = now in descending order as with date
eles_asnipe2$in_order_date <- c('yes', rep(NA,nrow(eles_asnipe2)-1)) # first is in order of date, rest NA
eles_asnipe2$in_order_time <- c('yes', rep(NA,nrow(eles_asnipe2)-1)) # first is in order of time, rest NA
for(i in 2:nrow(eles_asnipe2)){ # go through each and see if it is in order relative to the previous
  eles_asnipe2$in_order_date[i] <- ifelse(eles_asnipe2$Date[i] > eles_asnipe2$Date[i-1], 'no', 'yes')
  eles_asnipe2$in_order_time[i] <- ifelse(eles_asnipe2$Date[i] == eles_asnipe2$Date[i-1], 
                                          ifelse(eles_asnipe2$Time[i] <= eles_asnipe2$Time[i-1], 'yes', 'no'),
                                          'yes')
}
table(eles_asnipe2$in_order_date) ; table(eles_asnipe2$in_order_time) # sightings now in reverse chronological order, both by date and time

(N <- length(unique(eles_asnipe2$observation)))   # 1437
length(unique(efa$encounter)) - N                 # -30 -- encounters were not reliably numbered throughout dataset
eles_asnipe2$group_id_new <- c(N, rep(NA, nrow(eles_asnipe2)-1)) # first one in data is last sighting (highest number)
for(i in 2:nrow(eles_asnipe2)){                   # new value for each sighting (start at 1437, down to 1)
  eles_asnipe2$group_id_new[i] <- ifelse(eles_asnipe2$observation[i] == eles_asnipe2$observation[i-1],
                                         eles_asnipe2$group_id_new[i-1],
                                         eles_asnipe2$group_id_new[i-1]-1)
}
summary(eles_asnipe2$group_id_new)                # min = 1, max = 1437 -- this has worked
length(unique(eles_asnipe2$group_id_new))

observations <- eles_asnipe2[,c('Date','Time','Location','group_id_new')] %>% distinct() # 1437 observations

## SO... to identify which sightings in main data frame are actually represented in binomial data (all): merge 'observations' data with Binomial social interactions
rm(eles_asnipe, eles_asnipe2, N, asnipe, i, obs_id) ; gc()

#### merge sightings information into dyad data #####
colnames(observations) <- c('date', 'time', 'location', 'obs_id') # rename columns
all <- left_join(x = all, y = observations, by = 'obs_id')        # merge sightings info per dyad
head(all,10) ; tail(all,10)                                       # check structure
rm(observations) ; gc()

#### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
#all$dyad <- paste(all$node_1, all$node_2, sep = '_')              # create unique value per dyad
#all$dyad_id <- as.integer(as.factor(all$dyad))                    # convert to integer
all$location_id <- as.integer(as.factor(all$location))            # create integer value for GPS locations
all <- all %>% distinct()                                         # remove duplicates wherever there is an event
head(all)                                                         # check structure

all <- all[which(is.na(all$date) == FALSE),]
max(all$date) - min(all$date) # 2966 days. Time period for ALERT data = 504 days --> split into 6 time windows (7 boundaries)
periods <- seq(from = min(all$date), to = max(all$date), length.out = 7) # create date sequence split into chunks of similar length to MOTNP dataset
periods[7] <- periods[7]+1 # set to be one higher than the final higher otherwise it takes the last date and creates a whole new period
events <- all[all$social_event == 1,] ; events <- events[!is.na(events$node_1),] # cuts down to 5047 pairs from 3248447 dyad pairs -- assigning to time window doesn't take forever
events$period <- NA          # create new variable for ID of time window
for(i in 1:nrow(events)){    # fill time window variable
  events$period[i] <- which(periods <= events$date[i])[length(which(periods <= events$date[i]))] # fill cell with last value in vector
}
range(events$obs_id)    # check output
table(events$period)    # check output

periods ; #View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

# check elephants all match up
length(unique(all$node_1)) ; length(unique(all$node_2))         # check contain the same number of elephants
length(unique(events$node_1)) ; length(unique(events$node_2))   # check contain the same number of elephants -- don't need to

#### convert to Binomial model data format ####
df_split <- events %>%
  group_by(node_1, node_2, period) %>%          # aggregate all sightings of each dyad
  summarise(event_count = sum(social_event))
head(df_split)                                  # check output
rm(events) ; gc()

### add ID numbers
#eles <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long.csv') %>% 
#  select(elephant) %>%                                              # take only elephant ID column
#  distinct()                                                        # select unique values only
eles <- read_csv('../data_processed/mpnp_eles_long_Jan23check.csv') %>% 
  select(elephant) %>%                                              # take only elephant ID column
  distinct()                                                        # select unique values only
eles$node_1 <- as.integer(as.factor(eles$elephant))                 # create integer variable
colnames(eles)[1] <- c('id_1')                                      # rename variable
df <- left_join(df_split, eles, by = 'node_1') %>% distinct()       # add node ID to data frame about elephant 1
colnames(eles) <- c('id_2','node_2')                                # rename variables
df <- left_join(df, eles, by = 'node_2') %>% distinct()             # add node ID to data frame about elephant 2
head(df) ; tail(df)                                                 # check output
table(df$event_count)                                               # seems sensible
rm(df_split) ; gc()

#### create dyad row for all pairs per period ####
dyads <- data.frame(id_1 = rep(sort(eles$id_2), each = nrow(eles)),   # new data frame containing 1 row per dyad (duplicated so id1_id2 and id2_id1 are both present, as are dyads where id1 = id2)
                    id_2 = rep(sort(eles$id_2), nrow(eles)))

colnames(eles) <- c('id_1','node_1')            # rename variables for merging
dyads <- left_join(dyads, eles, by = 'id_1')    # merge elephant information into dyad data frame
colnames(eles) <- c('id_2','node_2')            # rename variables for merging
dyads <- left_join(dyads, eles, by = 'id_2')    # merge elephant information into dyad data frame
dyads <- dyads[dyads$node_1 < dyads$node_2,]    # remove any rows where id1 = id2 and duplicates (id2_id1)
rm(eles, i) ; gc()

dyads <- data.frame(id_1 = rep(dyads$id_1, length(unique(df$period))),        # reform data frame to the structure desired
                    id_2 = rep(dyads$id_2, length(unique(df$period))),
                    node_1 = rep(dyads$node_1, length(unique(df$period))),
                    node_2 = rep(dyads$node_2, length(unique(df$period))),
                    period = rep(sort(unique(df$period)), each = nrow(dyads)))

head(df) ; head(dyads)                                                # check structures of data frames
data <- left_join(x = dyads, y = df, by = c('id_1','id_2','period'))  # merge by dyad per time window
head(data)                                                            # check structure
data <- data[,c(1:5,8)]                                               # select desired variables
colnames(data)[3:4] <- c('node_1','node_2')                           # rename variables
data$event_count <- ifelse(is.na(data$event_count) == TRUE, 0, data$event_count)   # set NA counts (dyad not observed) to 0
table(data$event_count)                                               # check numbers of times dyads were observed
table(df$event_count)                                                 # check data against df -- match
data$dyad <- paste(data$id_1, data$id_2, sep = '_')                   # create character variable with unique value per dyad
data$dyad_id <- as.integer(as.factor(data$dyad))                      # create numeric variable with unique value per dyad
head(data, 20)                                                        # check structure

# identify start and end dates for each time window
efa <- readxl::read_excel('../data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx') %>% janitor::clean_names()  # read in sightings data
efa$date <- as_date(efa$date)  # obtain date values
efa <- efa %>% separate(elephant_id, into = c('letter','number'), sep = 1, remove = FALSE) %>% 
  filter(letter != 'T') %>% filter(letter != 'F') %>% filter(elephant_id != '-') # remove unidentified individuals
max(all$date)-min(all$date) ; max(efa$date)-min(efa$date)
windows <- data.frame(period_start = seq(from = min(efa$date), to = max(efa$date)+1, length.out = 7)[1:6],
                      period = 1:6,
                      period_days = periods[1:6])  # identify time windows

## check all time windows same length
windows$period_start[1]+(windows$period_days[2]-windows$period_days[1])
windows$period_start[2]+(windows$period_days[3]-windows$period_days[2])
windows$period_start[3]+(windows$period_days[4]-windows$period_days[3])
windows$period_start[4]+(windows$period_days[5]-windows$period_days[4])
windows$period_start[5]+(windows$period_days[6]-windows$period_days[5])

data <- left_join(x = data, y = windows, by = 'period')  # merge time window ID with event data
head(data)

## clean environment
rm(df, dyads, efa, periods) ; gc()

#### write to file ####
saveRDS(data, '../data_processed/mpnp_bayesian_pairwiseevents.RDS')

## clean up
rm(all, eles_asnipe_gbi) ; gc()

########## Short time windows ##########
#### add additional information and remove impossible sightings ##########
#data <- readRDS('../data_processed/mpnp_bayesian_pairwiseevents.RDS')

#### break down into periods so easier to manipulate
table(data$period)
mpnp1 <- data[data$period == 1,]  # 1st time window
mpnp2 <- data[data$period == 2,]  # 2nd time window
mpnp3 <- data[data$period == 3,]  # 3rd time window
mpnp4 <- data[data$period == 4,]  # 4th time window
mpnp5 <- data[data$period == 5,]  # 5th time window
mpnp6 <- data[data$period == 6,]  # 6th time window
#mpnp7 <- data[data$period == 7,]  # 7th time window
rm(data) ; gc()

## insert column for sighting period
#sightings <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long.csv')
sightings <- read_csv('../data_processed/mpnp_eles_long_Jan23check.csv')
sightings <- sightings[!is.na(sightings$elephant),c(2,3)]  # select only elephant ID and date columns
sightings$date <- as.Date(sightings$date, '%d/%m/%Y')      # numeric date
sightings$day <- as.numeric(sightings$date)                # numeric date
sightings$day <- sightings$day - (min(sightings$day)-1)    # convert date to start at 1 on first day of study
summary(sightings$day)

periods <- windows$period_days

sightings$period <- NA                                     # create new variable for the time window
for(i in 1:nrow(sightings)){
  sightings$period[i] <- which(periods <= sightings$day[i])[length(which(periods <= sightings$day[i]))] # set time window value as last value in vector
}
periods ; #View(sightings[sample(1:nrow(sightings),20),])         # visual check that individual dates match time windows

# create data frame of counts per dyad per time window
counts <- data.frame(id = rep(unique(sightings$elephant), 6),
                     period = rep(1:6, each = length(unique(sightings$elephant))),
                     count_all = NA,
                     count_period = NA)
for(i in 1:nrow(counts)){
  individual <- sightings[sightings$elephant == counts$id[i],]    # create data frame of all sightings per individual
  counts$count_all[i] <- nrow(individual)                         # count all sightings = number of rows in produced data frame
  counts$count_period[i] <- nrow(individual[individual$period == counts$period[i],])  # count sighitngs within time window = number of rows where time window matches
}
rm(individual, i) ; gc()

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
#data7 <- left_join(x = mpnp7, y = counts, by = c('id_1','period')) ; rm(mpnp7)  # add id1 count data for 7th time window
gc()

data1 <- data1[data1$count_period > 0,]   # remove all dyads id1 was never seen in 1st time window
data2 <- data2[data2$count_period > 0,]   # remove all dyads id1 was never seen in 2nd time window
data3 <- data3[data3$count_period > 0,]   # remove all dyads id1 was never seen in 3rd time window
data4 <- data4[data4$count_period > 0,]   # remove all dyads id1 was never seen in 4th time window
data5 <- data5[data5$count_period > 0,]   # remove all dyads id1 was never seen in 5th time window
data6 <- data6[data6$count_period > 0,]   # remove all dyads id1 was never seen in 6th time window
#data7 <- data7[data7$count_period > 0,]   # remove all dyads id1 was never seen in 7th time window

colnames(data1)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data2)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data3)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data4)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data5)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
colnames(data6)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging
#colnames(data7)[11:12] <- c('count_all_1','count_period_1')                     # rename post merging

colnames(counts)[1] <- 'id_2'                                                   # rename ID column for merging
data1 <- left_join(x = data1, y = counts, by = c('id_2','period'))              # add id2 count data for 1st time window
data2 <- left_join(x = data2, y = counts, by = c('id_2','period'))              # add id2 count data for 2nd time window
data3 <- left_join(x = data3, y = counts, by = c('id_2','period'))              # add id2 count data for 3rd time window
data4 <- left_join(x = data4, y = counts, by = c('id_2','period'))              # add id2 count data for 4th time window
data5 <- left_join(x = data5, y = counts, by = c('id_2','period'))              # add id2 count data for 5th time window
data6 <- left_join(x = data6, y = counts, by = c('id_2','period'))              # add id2 count data for 6th time window
#data7 <- left_join(x = data7, y = counts, by = c('id_2','period'))              # add id2 count data for 7th time window

data1 <- data1[data1$count_period > 0,]   # remove all dyads id2 was never seen in 1st time window
data2 <- data2[data2$count_period > 0,]   # remove all dyads id2 was never seen in 2nd time window
data3 <- data3[data3$count_period > 0,]   # remove all dyads id2 was never seen in 3rd time window
data4 <- data4[data4$count_period > 0,]   # remove all dyads id2 was never seen in 4th time window
data5 <- data5[data5$count_period > 0,]   # remove all dyads id2 was never seen in 5th time window
data6 <- data6[data6$count_period > 0,]   # remove all dyads id2 was never seen in 6th time window
#data7 <- data7[data7$count_period > 0,]   # remove all dyads id2 was never seen in 7th time window

colnames(data1)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data2)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data3)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data4)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data5)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
colnames(data6)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging
#colnames(data7)[13:14] <- c('count_all_2','count_period_2')                     # rename post merging

rm(counts,sightings) ; gc()

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

#table(data7$count_period_1)  # check id1 counts for 7th time window
#table(data7$count_period_2)  # check id2 counts for 7th time window

## count individuals
length(unique(data1$id_1))  # 650
length(unique(data2$id_1))  # 539
length(unique(data3$id_1))  # 172
length(unique(data4$id_1))  # 108
length(unique(data5$id_1))  # 62
length(unique(data6$id_1))  # 16
#length(unique(data7$id_1)) # 

#### add variable for sightings apart ####
data1$count_dyad <- (data1$count_period_1 + data1$count_period_2) - data1$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data1$apart <- data1$count_dyad - data1$event_count               # seen apart = total sightings - seen together
colnames(data1)[6] <- 'together'                                        # rename variable
table(data1$together)                                                   # check values
table(data1$apart)                                                      # check values

data2$count_dyad <- (data2$count_period_1 + data2$count_period_2) - data2$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data2$apart <- data2$count_dyad - data2$event_count               # seen apart = total sightings - seen together
colnames(data2)[6] <- 'together'                                        # rename variable
table(data2$together)                                                   # check values
table(data2$apart)                                                      # check values

data3$count_dyad <- (data3$count_period_1 + data3$count_period_2) - data3$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data3$apart <- data3$count_dyad - data3$event_count               # seen apart = total sightings - seen together
colnames(data3)[6] <- 'together'                                        # rename variable
table(data3$together)                                                   # check values
table(data3$apart)                                                      # check values

data4$count_dyad <- (data4$count_period_1 + data4$count_period_2) - data4$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data4$apart <- data4$count_dyad - data4$event_count               # seen apart = total sightings - seen together
colnames(data4)[6] <- 'together'                                        # rename variable
table(data4$together)                                                   # check values
table(data4$apart)                                                      # check values

data5$count_dyad <- (data5$count_period_1 + data5$count_period_2) - data5$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data5$apart <- data5$count_dyad - data5$event_count               # seen apart = total sightings - seen together
colnames(data5)[6] <- 'together'                                        # rename variable
table(data5$together)                                                   # check values
table(data5$apart)                                                      # check values

data6$count_dyad <- (data6$count_period_1 + data6$count_period_2) - data6$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data6$apart <- data6$count_dyad - data6$event_count               # seen apart = total sightings - seen together
colnames(data6)[6] <- 'together'                                        # rename variable
table(data6$together)                                                   # check values
table(data6$apart)                                                      # check values

#data7$count_dyad <- (data7$count_period_1 + data7$count_period_2) - data7$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
#data7$apart <- data7$count_dyad - data7$event_count               # seen apart = total sightings - seen together
#colnames(data7)[7] <- 'together'                                        # rename variable
#table(data7$together)                                                   # check values
#table(data7$apart)                                                      # check values

#### save data #####
data1 <- data1[,c(8,7,1:5,9,6,16,15,12,14)]  # select only necessary variables for 1st time window
data2 <- data2[,c(8,7,1:5,9,6,16,15,12,14)]  # select only necessary variables for 2nd time window
data3 <- data3[,c(8,7,1:5,9,6,16,15,12,14)]  # select only necessary variables for 3rd time window
data4 <- data4[,c(8,7,1:5,9,6,16,15,12,14)]  # select only necessary variables for 4th time window
data5 <- data5[,c(8,7,1:5,9,6,16,15,12,14)]  # select only necessary variables for 5th time window
data6 <- data6[,c(8,7,1:5,9,6,16,15,12,14)]  # select only necessary variables for 6th time window
#data7 <- data7[,c(8,7,1:5,9,6,16,15,12,14)]  # select only necessary variables for 7th time window

# write data files
write_csv(data1,'../data_processed/mpnp_period1_pairwiseevents.csv')
write_csv(data2,'../data_processed/mpnp_period2_pairwiseevents.csv')
write_csv(data3,'../data_processed/mpnp_period3_pairwiseevents.csv')
write_csv(data4,'../data_processed/mpnp_period4_pairwiseevents.csv')
write_csv(data5,'../data_processed/mpnp_period5_pairwiseevents.csv')

length(unique(data6$id_1))   # not usable data -- too few identified individuals
#write_csv(data6,'../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period6_pairwiseevents.csv')

#length(unique(data7$id_1))   # not usable data -- too few identified individuals
#write_csv(data7,'../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_period7_pairwiseevents.csv')

# clean up
rm(list=ls()[ls()!= 'windows']) ; gc()

########## Long time window -- first 5 short windows (can't just combine as have removed all where one or another was absent, so total count will be wrong) ##########
#### remove sightings outside long window ####
data <- readRDS('../data_processed/mpnp_bayesian_pairwiseevents.RDS')

nrow(data)/length(unique(data$period))
cumsum(1:length(unique(data$id_1)))[length(unique(data$id_1))] == length(unique(data$dyad_id)) # correct number of IDs and dyads

#### cut out last time window
table(data$period)
mpnp_long <- data[data$period < 6, c(1:4,6:8)] %>% distinct()
rm(data) ; gc()

## insert column for sighting period
#sightings <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/mpnp_eles_long.csv')
sightings <- read_csv('../data_processed/mpnp_eles_long_Jan23check.csv')
sightings <- sightings[!is.na(sightings$elephant),c(2,3)]  # select only elephant ID and date columns
sightings$date <- as.Date(sightings$date, '%d/%m/%Y')      # numeric date
sightings$day <- as.numeric(sightings$date)                # numeric date
sightings$day <- sightings$day - (min(sightings$day)-1)    # convert date to start at 1 on first day of study
summary(sightings$day)

periods <- windows$period_days

sightings$period <- NA                                     # create new variable for the time window
for(i in 1:nrow(sightings)){
  sightings$period[i] <- which(periods <= sightings$day[i])[length(which(periods <= sightings$day[i]))] # set time window value as last value in vector
}
periods ; #View(sightings[sample(1:nrow(sightings),20),])         # visual check that individual dates match time windows

# remove sightings after period 5
sightings <- sightings[sightings$period < 6,]

#### count elephant sightings within long window ####
# create data frame of counts per dyad per time window
counts <- data.frame(id = unique(sightings$elephant),
                     count = NA)
for(i in 1:nrow(counts)){
  individual <- sightings[sightings$elephant == counts$id[i],]    # create data frame of all sightings per individual
  counts$count[i] <- nrow(individual)                             # count all sightings = number of rows in produced data frame
}
rm(individual, i) ; gc()

summary(counts$count)      # check output

# join counts data frame with dyad data
colnames(counts)[1] <- 'id_1'                                              # rename ID column for merging
data <- left_join(x = mpnp_long, y = counts, by = 'id_1') ; rm(mpnp_long)  # add id1 count data
colnames(data)[8] <- 'count_1'                                             # rename post merging
colnames(counts)[1] <- 'id_2'                                              # rename ID column for merging
data <- left_join(x = data, y = counts, by = 'id_2')                       # add id2 count data for 1st time window
colnames(data)[9] <- 'count_2'                                             # rename post merging

rm(counts,sightings) ; gc()

table(data$count_1)  # check id1 counts for 1st time window
table(data$count_2)  # check id2 counts for 1st time window

#### add variable for sightings apart ####
data$count_dyad <- (data$count_1 + data$count_2) - data$event_count   # total sightings = total of both but remove duplicates where seen together (A + B - AB)
data$apart <- data$count_dyad - data$event_count                      # seen apart = total sightings - seen together
colnames(data)[5] <- 'together'                                       # rename variable
table(data$together)                                                  # check values
table(data$apart)                                                     # check values

#### save data #####
data <- data[,c(6,7,1:5,11,8:10)]  # rearrange

# write data files
write_csv(data,'../data_processed/mpnp_longtimewindow_pairwiseevents.csv')
