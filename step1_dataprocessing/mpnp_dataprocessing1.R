### Bayesian analysis of EFA data ####
# Script to process association data from Makgadikgadi Pans National Park, Botswana.
# Data collected: 1st November 2019 to 5th August 2021
# Collected by: Elephants for Africa (Dr Kate Evans)
# Data supplied by: EfA and Dr Kate Evans, 19th October 2021

### set up ####
# load packages
#library(tidyverse) ; library(dplyr) ; library(lubridate) ; library(janitor) ; library(hms) ; library(readxl)
library(tidyverse, '../packages/')
library(dplyr, '../packages/')
library(lubridate, '../packages/')
library(janitor, '../packages/')
library(hms, '../packages/')
library(readxl, '../packages/')

# set seed
set.seed(12345)

#### sightings data frame ####
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

# checking that all count columns within dataframe "s" are actually being used
check_na <- data.frame(col = 24:56,del = NA)
for(i in 1:nrow(check_na)) {
  check_na$del[i] <- nrow(s) - length(which(is.na(s[,i]) == 'FALSE'))
} # no columns to delete -- all have at least one recording
rm(check_na)

# individual data
#efa <- readxl::read_excel('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')
efa <- readxl::read_excel('../data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')
str(efa)
efa$time_cat <- lubridate::hour(efa$Time) # create variable which is just hour of the day
efa <- separate(efa, Time, into = c('wrong_date','time'), sep = ' ') # time column registers right time but wrong date
efa$time <- hms::as_hms(efa$time)                                    # extract time value only from time column
efa <- efa[,c(1:5,7,33,8:32)]                                        # remove incorrect date column
efa$Age_Range_ID <- as.factor(efa$Age_Range_ID)                      # make factor not character
efa$Activity_ID <- as.factor(efa$Activity_ID)                        # make factor not character
efa$Distance_To_Observer <- as.numeric(efa$Distance_To_Observer)     # make number not character
efa$Physical_Condition_ID <- as.factor(efa$Physical_Condition_ID)    # make factor not character
efa <- efa[,c(1:15,18:27,31:32)]                                     # rearrange
efa$Date <- lubridate::as_date(efa$Date)                             # make date
efa <- janitor::clean_names(efa)                                     # clean up

# correct IDs -- some have b instead B that gets treated as 2 separate elephants
efa <- efa %>% separate(col = elephant_id, into = c('letter','number'), sep = 1, remove = F)
unique(efa$letter)
efa$elephant_id <- ifelse(efa$letter == 'b', paste0('B',efa$number), efa$elephant_id)

## taking a look -- this is just data exploration to check that I know what all of the columns are doing and get an overview of the data ####
table(efa$elephant_id)
hist(table(efa$elephant_sighting_id), breaks = 30, las = 1)
length(which(table(efa$elephant_sighting_id) == 1)) # 2400 elephants seen singly
length(which(table(efa$elephant_sighting_id) == 2)) # 1017 elephants seen as pairs

table(efa$age_range_id) # 2 = <5, 3 = 5-9, 4 = 10-15, 5 = 16-20, 6 = 21-25, 7 = 26-35, 8 = 36+, 10 = Unknown (most are 4-7)

# reaction index -- check proportion of reactions
summary(efa$reaction_indices)
efa$reaction_indices <- factor(efa$reaction_indices, levels = c(1:3))
table(efa$reaction_indices) # none come towards the truck
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

# add additional encounter numbers for those that are missing
efa$encounter <- paste(efa$elephant_sighting_id, efa$date, efa$time, sep = '_')                     # single variable unique for all encounters
for(i in 1:nrow(efa)){                                                                              # set up for loop
  if(is.na(efa$elephant_sighting_id[i]) == TRUE){                                                   # only run loop if 
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

for(i in sort(unique(efa_long$encounter))) {
  encounter <- efa_long[efa_long$encounter == i,]
  efa_long$total_elephants[efa_long$encounter == i] <- nrow(encounter)                  # total eles
  efa_long$total_id[efa_long$encounter == i] <- length(which(encounter$letter == 'B'))  # total identified
}
rm(encounter, x, i) ; gc()
efa_long$perc_id <- 100*(efa_long$total_id/efa_long$total_elephants) # calculate proportion of elephants identified from total
summary(efa_long$perc_id)
mean(efa_long$total_elephants) ; sd(efa_long$total_elephants)

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

# save data
write_delim(efa_gps, '../data_processed/mpnp_eles_long_Jan23check.csv', delim = ',') # write to file
# efa <- read_csv('../data_processed/mpnp_eles_long_Jan23check.csv')
efa <- efa_gps
rm(efa_gps, efa_long, s, sighting_id) ; gc()

### nodes data frame ####
ele_nodes <- data.frame(id_no = sort(unique(efa$elephant)),   # create data frame for all information about an individual
                        age_class = NA,
                        age_category = NA,
                        count = NA,
                        dem_class = NA)

efa$age_range[which(efa$elephant == ele_nodes$id_no[1])]  # 6, 7 & 8 + 1UK (10) -- possible age classes for ID B1000
efa_age <- efa[,c(1,2,10)]                                # sighting, individual and age estimate
efa_age$age_range <- as.numeric(efa$age_range)            # make numeric
table(efa$age_range) ; table(efa_age$age_range)           # has converted 10 to 1
efa_age$age_range <- ifelse(efa_age$age_range == 1, 10, efa_age$age_range) # revert to value of 10 to avoid confusion
str(efa_age)

### mapping ####
plot(efa$longitude ~ efa$latitude, las = 1, xlab = 'latitude', ylab = 'longitude', pch = 19, col = rgb(0,0,1,0.1))
efa$longitude[which(efa$longitude < 24)] <- NA  # anything below 24 degrees is too far west
efa$latitude[which(efa$latitude < -20)]  <- NA  # anything closer to zero than -20 is too far north
plot(efa$longitude ~ efa$latitude, ylim = c(24.2,24.9), xlim = c(-20.8,-20.2), # replot without weird datapoints
     pch = 19, col = rgb(0,0,1,0.05), xlab = 'latitude', ylab = 'longitude', las = 1)

#### create dyadic data frame ####
## create gbi_matrix ####
### create group-by-individual matrix
eles_asnipe <- efa[,c(3,4,2,14)]                                # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                # convert to integer
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)  # start from 1, not 1st January 1970
eles_asnipe$Time <- hour(eles_asnipe$time)*60*60 + minute(eles_asnipe$time)*60 + second(eles_asnipe$time) # convert time values to seconds through day
eles_asnipe <- eles_asnipe[,c(5,6,3,4)]                         # select desired variables
colnames(eles_asnipe) <- c('Date','Time','ID','Location')       # rename for asnipe package
str(eles_asnipe)

# get_gbi generates a group by individual matrix. The function accepts a data.table with individual identifiers and a group column. The group by individual matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 4, pad = '0')    # ensure dates go in the right order even when as character
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value per sighting
eles_asnipe$observation <- as.integer(as.factor(eles_asnipe$encounter))  # unique count value per encounter
max(eles_asnipe$observation)                                             # 1437 (3765 if un-ID'd included)
eles_asnipe <- eles_asnipe[,c(3,7)]                                # ID and observation
eles_asnipe <- data.table::setDT(eles_asnipe)                      # convert to data table
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'observation', id = 'ID') # create group-by-individual matrix

## convert gbi_matrix to series of dyadic data frames ####
nrow(gbi_matrix) # 1437

# create observation number variable in full data set
eles_asnipe$observation <- as.integer(as.factor(eles_asnipe$observation))

# break down observations into blocks to run as parallel chunks
block_size <- 50 # set the number of observations per block to run at once
factor_levels <- 1:length(unique(eles_asnipe$observation)) # create a vector of numbers from 1 to the number of levels of the factor variable
factor_blocks <- split(factor_levels, ceiling(factor_levels / block_size)) # split the factor levels into blocks of the specified size

# create variable to identify which block number that observation falls into
block_table <- data.frame(block1 = factor_blocks[[1]])
for(i in 2:length(factor_blocks)){
  if(nrow(block_table) == length(factor_blocks[[i]])){
    block_table <- cbind(block_table, factor_blocks[[i]])
    colnames(block_table)[i] <- as.character(i)
  } else {
    length_na <- nrow(block_table) - length(factor_blocks[[i]])
    block_table <- cbind(block_table, c(factor_blocks[[i]], rep(NA, length_na)))
    colnames(block_table)[i] <- as.character(i)
  }
}
block_table <- pivot_longer(block_table, cols = everything(), names_to = 'block', values_to = 'observation')
block_table$block <- ifelse(block_table$block == 'block1', '1', block_table$block)
block_table$block <- as.numeric(block_table$block)
block_table

eles_asnipe <- left_join(eles_asnipe, block_table, by = "observation")

# create empty data frame
#nrows <- (length(unique(eles_asnipe$ID))-1) * nrow(eles_asnipe)
#gbi_df <- data.frame(node_1 = rep(NA, nrows),
#                     node_2 = rep(NA, nrows),
#                     social_event = rep(NA, nrows),
#                     obs_id = rep(NA, nrows),
#                     block = rep(NA, nrows))

# load library, create cluster
lapply(c("foreach", "doParallel"), require, character.only = TRUE)
library(parallel)
cl <- makeCluster(detectCores()) # can replace detectCores() with a number if known
registerDoParallel(cl)

# Function to be parallelised
process_blocks <- function(block_data) {
  num_rows <- (length(unique(eles_asnipe$ID))-1) * nrow(block_data)
  # Pre-allocate list (might need to make dataframe)
  rows <- data.frame( node_1 = rep(NA, num_rows),
                      node_2 = rep(NA, num_rows),
                      id_1 = rep(NA, num_rows),
                      id_2 = rep(NA, num_rows),
                      social_event = rep(NA, num_rows),
                      obs_id = rep(NA, num_rows),
                      block = rep(NA, num_rows))
  
  for (obs_id in min(block_data$observation):max(block_data$observation)) {
    #row_index <- 1
    for (i in which(gbi_matrix[obs_id, ] == 1)) {
      for (j in 1:ncol(gbi_matrix)) {
        if (i != j) {
          node_1 <- ifelse(i < j, i, j)
          node_2 <- ifelse(i < j, j, i)
          id_1 <- ifelse(i < j, colnames(gbi_matrix)[i], colnames(gbi_matrix)[j])
          id_2 <- ifelse(i < j, colnames(gbi_matrix)[j], colnames(gbi_matrix)[i])
          
          empty_row <- which(is.na(rows$node_1) == TRUE)[1]
          rows[empty_row, 1] <- node_1
          rows[empty_row, 2] <- node_2
          rows[empty_row, 3] <- id_1
          rows[empty_row, 4] <- id_2
          rows[empty_row, 5] <- ifelse(gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j], 1, 0)
          rows[empty_row, 6] <- obs_id
          rows[empty_row, 7] <- block_id
        }
      }
    }
    #row_index <- row_index + 1
    if(obs_id %% 10 == 0) {print(paste0('obs_id ', obs_id, ' in block ', block_id, ' finished at time ', Sys.time()))}
  }
  return(rows)
}

#### run loop ####
sorted_unique_blocks <- sort(unique(as.numeric(eles_asnipe$block)))
num_unique_blocks <- length(sorted_unique_blocks)

eles_asnipe$block <- as.character(eles_asnipe$block)

foreach(block_id = 1:num_unique_blocks) %dopar% {
  # break down data by block
  block_data <- eles_asnipe[eles_asnipe$block == block_id,]
  # set up data frame to write into
  #nrows <- (length(unique(block_data$ID))-1) * nrow(block_data)
  #gbi_df <- data.frame(node_1 = rep(NA, nrows),
  #                     node_2 = rep(NA, nrows),
  #                     id_1 = rep(NA, nrows),
  #                     id_2 = rep(NA, nrows),
  #                     social_event = rep(NA, nrows),
  #                     obs_id = rep(NA, nrows),
  #                     block = rep(NA, nrows))
  # run process_blocks
  gbi_df <- process_blocks(block_data)
  # trim down
  #gbi_df <- distinct(gbi_df)
  # save output
  saveRDS(gbi_df, paste0('../data_processed/mpnp_bayesian_allpairwiseevents_block',block_id,'_Feb23check.RDS'))
  #gbi_df
}

# run checks
check <- readRDS('../data_processed/mpnp_bayesian_allpairwiseevents_block2_Feb23check.RDS') %>% 
  distinct() %>% 
  filter(is.na(node_1) == FALSE)
check_matrix <- as.data.frame(gbi_matrix[51:100,])
check_matrix$group_size <- rowSums(check_matrix)
table(check_matrix$group_size)
sum(check_matrix$group_size)
sum(check$social_event)

# visual check who is present in each block against listed social events
for(i in 1:nrow(check_matrix)) {print(colnames(check_matrix)[which(check_matrix[i,] == 1)])}
View(check)
