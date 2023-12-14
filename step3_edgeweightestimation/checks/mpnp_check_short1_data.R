library(tidyverse) ; library(dplyr) ; library(cmdstanr) ; library(igraph) ; library(janitor) ; library(lubridate) ; library(hms) ; library(readxl)

#### import data -- window 1 ####
### set up
time_window <- 1

### subset by time window
cdf1 <- read_csv(paste0('../data_processed/step1_dataprocessing/mpnp_period',time_window,'_pairwiseevents.csv'))

### set values for model
n1 <- nrow(cdf1)

### create data list
cls1 <- list(
  n_dyads    = n1,                  # total number of times one or other of the dyad was observed
  dyad_ids   = cdf1$dyad_id,        # identifier for each dyad
  together   = cdf1$together,       # count number of sightings seen together
  count_dyad = cdf1$count_dyad      # count total number of times seen
)

#### check for anomalies ####
## check for missing data
for(i in 1:length(cls1)){
  print(length(which(is.na(cls1[i]) == TRUE)))
}
print(length(which(is.nan(cdf1$dyad_id) == TRUE)))
print(length(which(is.nan(cdf1$together) == TRUE)))
print(length(which(is.nan(cdf1$count_dyad) == TRUE)))

## check for data where sightings together > sightings total
length(which(cls1$together > cls1$count_dyad))

## check for weird values (e.g. negatives)
barplot(table(cls1$together))
barplot(table(cls1$count_dyad))

## check high together values
11 / (cls1$count_dyad[cls1$together == 11])
 9 / (cls1$count_dyad[cls1$together ==  9])
 8 / (cls1$count_dyad[cls1$together ==  8])
 7 / (cls1$count_dyad[cls1$together ==  7])
 6 / (cls1$count_dyad[cls1$together ==  6])
 5 / (cls1$count_dyad[cls1$together ==  5])
 
## check for duplicate dyads or missing/extra individuals
nrow(cdf1) - length(unique(cdf1$dyad_id))
nrow(cdf1) - length(unique(paste0(cdf1$id_1, cdf1$id_2)))
id1 <- unique(cdf1$id_1) ; id2 <- unique(cdf1$id_2)
length(id1) ; length(id2)
id1 <- id1[2:length(id1)] ; id2 <- id2[1:(length(id2)-1)] ; length(which(id1 != id2))

#### import data -- window 5 for comparison ####
### set up
time_window <- 5

### subset by time window
cdf5 <- read_csv(paste0('../data_processed/step1_dataprocessing/mpnp_period',time_window,'_pairwiseevents.csv'))

### set values for model
n5 <- nrow(cdf5)
 
### create data list
cls5 <- list(
  n_dyads    = n5,                  # total number of times one or other of the dyad was observed
  dyad_ids   = cdf5$dyad_id,        # identifier for each dyad
  together   = cdf5$together,       # count number of sightings seen together
  count_dyad = cdf5$count_dyad      # count total number of times seen
)
 
#### check for anomalies ####
## check for missing data
for(i in 1:length(cls5)){
  print(length(which(is.na(cls5[i]) == TRUE)))
}

## check for data where sightings together > sightings total
length(which(cls5$together > cls5$count_dyad))

## check for weird values (e.g. negatives)
barplot(table(cls5$together))
barplot(table(cls5$count_dyad))

## check high together values
1 / (cls5$count_dyad[cls5$together == 1]) # most sightings for apair together at any time is 5

#### import data -- long time window for comparison ####
cdf_long <- read_csv('../data_processed/step1_dataprocessing/mpnp_longtimewindow_pairwiseevents.csv')

### set values for model
nlong <- nrow(cdf_long)

### create data list
cls_long <- list(
  n_dyads    = nlong,                  # total number of times one or other of the dyad was observed
  dyad_ids   = cdf_long$dyad_id,        # identifier for each dyad
  together   = cdf_long$together,       # count number of sightings seen together
  count_dyad = cdf_long$count_dyad      # count total number of times seen
)

#### check for anomalies ####
## check for missing data
for(i in 1:length(cls_long)){
  print(length(which(is.na(cls_long[i]) == TRUE)))
}

## check for data where sightings together > sightings total
length(which(cls_long$together > cls_long$count_dyad))

## check for weird values (e.g. negatives)
barplot(table(cls_long$together))
barplot(table(cls_long$count_dyad))

## check high together values
11 / (cls_long$count_dyad[cls_long$together == 11])
9 / (cls_long$count_dyad[cls_long$together ==  9])
8 / (cls_long$count_dyad[cls_long$together ==  8])
7 / (cls_long$count_dyad[cls_long$together ==  7])
6 / (cls_long$count_dyad[cls_long$together ==  6])
5 / (cls_long$count_dyad[cls_long$together ==  5])

