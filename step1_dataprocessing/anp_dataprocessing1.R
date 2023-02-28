#### information ####
# Data collected by Amboseli Trust for Elephants (ATE) 1972-2021
# Data supplied by Vicki Fishlock, 24th February 2022

#### set up ####
#library(tidyverse, lib.loc = '../packages/')
#library(lubridate, lib.loc = '../packages/')
#library(janitor, lib.loc = '../packages/')
#library(hms, lib.loc = '../packages/')
#library(readxl, lib.loc = '../packages/')
#library(data.table, lib.loc = '../packages/')
#library(spatsoc, lib.loc = '../packages/')

library(tidyverse)
library(lubridate)
library(janitor)
library(hms)
library(readxl)
library(data.table)
library(spatsoc)

#### import sightings data ####
ate <- readxl::read_xlsx('../data_raw/Raw_ATE_MaleSightingsCleaned_Fishlock220808.xlsx') %>% 
  janitor::clean_names()# %>% 
  #select(-obs_casename, -row_num)
old <- readxl::read_xlsx('../data_raw/Raw_ATE_Sightings_Fishlock220224.xlsx') %>% janitor::clean_names()

str(ate)
str(old)
## hab_code_1 = habitat code for most of the group: 0 = Not described, 01 = short grass plain, 02 = Consimilis tall grasslands, 03 = Acacia tortilis woodland, 04 =	Acacia xanthophloea woodland, 05 = Salvadora/sueda, 06 = Palm woodland, 07 = Swamp edge woodland, 08 = Swamp edge, 09 =	Swamp, A =	Bushed grassland, B =	Open bushland north, C = Open bushland south, D	= Dense bushland north, E =	Acacia nubica, F = Dense bushland south, U = Unknown
## hab_code_2 = second habitat code if >1 category
## bull_q_r = quality of recognition for males (3 = knew all, 2 = knew > half, 1 = knew < half, U = unknown, blanks relate to data edits I haven't completed in the base data)
## obs_type = Group Type: B = male only, M = males and females, U = unknown
## grp_q_c = Group Quality of count: 3 = exact, 2 = good estimate, 1 = estimate, 0 = no estimate
## grp_size = N elephants encountered in group; -1 is missing data. Zeros are typos that I haven't had chance to correct at present
## bull_q_c = Male Quality of count: 3 = exact, 2 = good estimate, 1 = estimate, 0 = no estimate
## num_bulls = N independent males present. Zeros are typos I haven't had chance to correct at present
## bulls_1_2 = Males aged 10-24, present/absent (1 = present,0 = absent, U = unknown)
## bulls_3_5 = Males aged 25+, present/absent (1 = present,0 = absent, U = unknown)
## musth_male = Musth male present/ absent (1 = present,0 = absent, U = unknown)
## oestrus_fem = Oestrue female present/absent ((1 = present,0 = absent, U = unknown)
## act_code = Group activiity: 00 = not specified, 01 = walking while feeding, 02 = feeding, 03 = resting, 04 =	comfort (dusting, mudwallowing etc), 05 =	interacting, 06 =	drinking, 07 = walking, 08 = standing vigilant, 09 = more than one activity
ate$bull_q_c <- as.character(ate$bull_q_c)
ate$grp_q_c <- as.character(ate$grp_q_c)

locations <- old[,c('obs_id','casename','obs_num','utm_lat','utm_long','bull_q_r','obs_type')]
ate <- left_join(ate, locations, by = c('obs_id','casename'))

ate <- ate[,c('obs_id','casename','musth','obs_date','obs_time','obs_num','grid_code','utm_lat','utm_long','hab_code_1','hab_code_2','bull_q_r','obs_type','grp_q_c','grp_size','bull_q_c','num_bulls','bulls_1_2','bulls_3_5','musth_male','oestrus_fem','act_code')]

colnames(ate)[c(6,8,9,12,13)] <- c("obs_num_old","utm_lat_old","utm_long_old","bull_q_r_old","obs_type_old")

rm(locations,old)

## casename = MaleID number -- make character string obvious so clearly different from node_id later on
ate$id <- paste('M',ate$casename, sep = '')
ate$node_id <- as.integer(as.factor(ate$casename))

## obs_date = Date of observation -- convert to date value
ate$obs_date <- lubridate::as_date(ate$obs_date)

## obs_time = Time of observation -- convert to time value
ate <- separate(ate, obs_time, into = c('wrong_date','correct_time'), remove = F, sep = ' ')
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

ate$obs_num_old <- ifelse(ate$obs_num_old == '0','00', ate$obs_num_old)
ate$obs_num_old <- ifelse(ate$obs_num_old == '0a','0A', ate$obs_num_old)
ate$obs_num_old <- ifelse(ate$obs_num_old == '0b','0B', ate$obs_num_old)
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
no_gps <- unique(ate$obs_date[which(ate$utm_lat_old == 0)])
gps <- unique(ate$obs_date[which(ate$utm_lat_old != 0)])
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

## rearrange
ate <- ate[,c(1:3,25,26,4,27,28,8,29,9:24)]
colnames(ate)

## clean environment
rm(ate_nums, date_row, test, gps, i, no_gps, test_nums)

## save output
write_csv(ate, '../data_processed/anp_sightings_rawcombined.csv')

## progress report
print(paste0('sightings data imported at ', Sys.time()))

#### import nodes data ####
nodes_old <- readxl::read_excel('../data_raw/Raw_ATE_Males_Lee220121.xlsx') %>% janitor::clean_names()

nodes <- readxl::read_excel('../data_raw/Raw_ATE_LifeHistoryData_Fishlock220808.xlsx') %>% janitor::clean_names() %>% distinct()

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
### code to convert gbi matrix format to dyadic data frame, shared by Prof Dan Franks and available from @JHart96 GitHub repository (https://github.com/JHart96/bison_examples/blob/main/examples/convert_gbi.md)
# create observation number variable in full data set
ate$observation_number <- as.integer(as.factor(ate$obs_id))

## save workspace image
save.image('ate_dataprocessing1_preparallelisation.RData')
load('ate_dataprocessing1_preparallelisation.RData')

# set the size of the blocks
block_size <- 50

# create a vector of numbers from 1 to the number of levels of the factor variable
factor_levels <- 1:length(unique(ate$obs_id))

# split the factor levels into blocks of the specified size
factor_blocks <- split(factor_levels, ceiling(factor_levels / block_size))

# create variable in ate that identifies which block number that observation falls into
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
block_table <- pivot_longer(block_table, cols = everything(), names_to = 'block', values_to = 'observation_number')
block_table$block <- ifelse(block_table$block == 'block1', '1', block_table$block)
block_table$block <- as.numeric(block_table$block)
block_table

#ate_test <- left_join(ate_test, block_table, by = "observation_number")
ate <- left_join(ate, block_table, by = "observation_number")

# create empty data frame
nrows <- (length(unique(ate$id))-1) * nrow(ate)
gbi_df <- data.frame(node_1 = rep(NA, nrows),
                     node_2 = rep(NA, nrows),
                     social_event = rep(NA, nrows),
                     obs_id = rep(NA, nrows),
                     block = rep(NA, nrows))

# load library, create cluster
lapply(c("foreach", "doParallel"), require, character.only = TRUE)
library(parallel)
cl <- makeCluster(detectCores()) # can replace detectCores() with a number if known
registerDoParallel(cl)


# Function to be parallelised
process_blocks <- function(block_data) {
  num_rows <- (length(unique(ate$id))-1) * nrow(block_data)
  # Pre-allocate list (might need to make dataframe)
  rows <- data.frame( node_1 = rep(NA, num_rows),
                      node_2 = rep(NA, num_rows),
                      social_event = rep(NA, num_rows),
                      obs_id = rep(NA, num_rows),
                      block = rep(NA, num_rows))
  
  for (obs_id in min(block_data$observation_number):max(block_data$observation_number)) {
    #row_index <- 1
    for (i in which(gbi_matrix[obs_id, ] == 1)) {
      for (j in 1:ncol(gbi_matrix)) {
        if (i != j) {
          node_1 <- ifelse(i < j, i, j)
          node_2 <- ifelse(i < j, j, i)
          
          empty_row <- which(is.na(rows$node_1) == TRUE)[1]
          rows[empty_row, 1] <- node_1
          rows[empty_row, 2] <- node_2
          rows[empty_row, 3] <- ifelse(gbi_matrix[obs_id, i] == gbi_matrix[obs_id, j], 1, 0)
          rows[empty_row, 4] <- obs_id
          rows[empty_row, 5] <- block
        }
      }
    }
    #row_index <- row_index + 1
    #if(obs_id %% 10 == 0) {print(obs_id)}
    #if(obs_id %% 10 == 0) {print(Sys.time())}
    #print(paste0('obs_id ', obs_id, ' in block ', block, ' finished at time ', Sys.time()))
  }
  return(rows)
}

#### run loop ####
sorted_unique_blocks <- sort(unique(ate$block))
num_unique_blocks <- length(sorted_unique_blocks)

foreach(block_id = 1:num_unique_blocks) %dopar% {
  # break down data by block
  block <- sorted_unique_blocks[block_id]
  block_data <- ate[ate$block %in% block,]
  # set up data frame to write into
  nrows <- (length(unique(block_data$id))-1) * nrow(block_data)
  gbi_df <- data.frame(node_1 = rep(NA, nrows),
                       node_2 = rep(NA, nrows),
                       social_event = rep(NA, nrows),
                       obs_id = rep(NA, nrows),
                       block = rep(NA, nrows))
  # run process_blocks
  gbi_df <- process_blocks(block_data)
  process_blocks(block_data)
  # save output
  saveRDS(gbi_df, paste0('../data_processed/anp_bayesian_allpairwiseevents_block',block,'.RDS'))
  print(paste0('block ', block, ' finished at ', Sys.time()))
  gbi_df
}
