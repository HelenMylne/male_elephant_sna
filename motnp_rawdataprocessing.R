#### Information ####
# Script to process and plot association data from Mosi-Oa-Tunya National Park, Zambia.
# Data collected: 19th May 2016-21st December 2017
# Collected by: Mr David Youldon, Mr Dabwiso Sakala, Miss Helen Mylne and other volunteers/interns/facilitated research students working with ALERT during this time
# Data supplied by: ALERT and Mr David Youldon (11th August 2021) and via Volunteer Encounter (Bex Saunders, 19th October 2021)

#### Set up ####
library(tidyverse)  # data manipulation
library(lubridate)  # sort date columns out
library(zoo)        # sort date columns out
library(asnipe)     # generating networks
library(rethinking) # col.alpha
library(igraph)     # plotting networks
################ First Half of Script: process raw data into clean csv files ################
#### Observation Sessions ####
sessions <- readxl::read_excel("data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 1)
s <- sessions[,1:4]             # remove summary section
s <- janitor::clean_names(s)    # convert names to snake case
colnames(s)[4] <- 'mins'        # remove '4' from minutes name

s$start <- hms::as_hms(s$start) # select out times from date time cells. dates are wrong anyway.
s$end <- hms::as_hms(s$end)     # select out times from date time cells. dates are wrong anyway.

### categorise sessions: Morning = 5am-12am leave (up to 2pm return), All day = 6am leave return after 2pm, Afternoon = 12pm-6pm start, Evening = 6pm-10pm start, Night = 10pm-5am start
s$session <- ifelse(s$start >= 22*60*60, 'night', 
                    ifelse(s$start >= 18*60*60, 'evening', 
                           ifelse(s$start >= 12*60*60, 'afternoon',
                                  ifelse(s$end <= 5*60*60, 'night',
                                         ifelse(s$end <= 15*60*60, 'morning','day')))))

### create a new data frame giving the daily totals
days <- aggregate(s$mins, by = list(s$date), FUN = sum)
colnames(days) <- c('date','mins')

### write to csv
write_delim(days, 'data_processed/motnp_recording_days.csv',  delim = ';', col_names = T)
write_delim(s, 'data_processed/motnp_recording_sessions.csv', delim = ';', col_names = T)

### clear environment
rm(days,s,sessions)

#### IDs ####
### import data
id.raw <- readxl::read_excel("data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 2)
id <- id.raw[5:849,1:368]             # create a new data frame without top values or empty end columns
id <- id[!is.na(id[,1]),]             # remove empty rows
id <- id[!is.na(id[,5]),]             # remove unassigned numbers

### rename columns
colnames(id)[1:10] <- id.raw[4,1:10]  # rename first columns with the names from the original

months <- seq(from = dmy('1st May 2016'), to = dmy('1st Dec 2017'), by = 'month')  # create vector of days, 1st day of each month
#colnames(id)[11:30] <- months # label month columns with VALUE dates of 1st day of month
months <- paste('m',as.character(months), sep = '_')  # replace vector values with something readable
colnames(id)[11:30] <- months # label month columns with READABLE dates of 1st day

weeks <- seq(from = dmy('19th May 2016'), to = dmy('15th Nov 2017'), by = 'weeks') # create vector of days, every Thursday from start to end
#colnames(id)[31:108] <- weeks # label week columns with date VALUES of 1st day of each week, starting Thursdays
weeks <- paste('w',as.character(weeks), sep = '_') # replace weeks vector values with readable dates
colnames(id)[31:108] <- weeks       # label week columns with READABLE dates, starting each week on Thursdays

days <- id.raw[4,109:368]     # can't use sequence method because pattern of days is random
days.numbers <- data.frame(as.numeric(days))     # all fine except for last 3
colnames(days.numbers) <- 'nums'                 # rename for ease
days.numbers$dates <- as.Date(days.numbers$nums, origin = "1970-01-01") # mostly 70 years and 2 days ahead
days.numbers$dates.new <- as.Date(days.numbers$nums - (70*365.25+1),
                                  origin = "1970-01-01")                # adjust dates by 70y & 1d
days.numbers$dates.new[258:260] <- c(days.numbers$dates.new[257]+2,     # manually correct last 3
                                     days.numbers$dates.new[257]+5,     # manually correct last 3
                                     days.numbers$dates.new[257]+6)     # manually correct last 3
colnames(id)[109:368] <- days.numbers$dates.new-0.5                     # label day columns
days <- paste('d',as.character(days.numbers$dates.new-0.5), sep = '_')  # replace days vector values with readable dates
colnames(id)[109:368] <- days  # label day columns with readable dates

id <- janitor::clean_names(id) # convert names to snake case

### in columns of when elephant was seen or not, replace NA (not seen) with 0 and "Yes" (seen) with 1
for(i in 11:368){ 
  id[i] <- ifelse(is.na(id[i]), 0, 1)
}

### correct individual mistake -- elephant seen on 2 weeks but only 1 day
nkash <- id[434,]              # observed 20th July, reported in both weeks 13-19th July and 20th-26th
nkash$w_2017_07_13 <- 0        # remove sighting in week 13-19th July
nkash$weeks <- as.character(1) # correct count of number of weeks elephant observed
id[434,] <- nkash              # replace data row with corrected vector

### correct data types -- replace Excel counts with R counts
str(id)
id$months <- as.numeric(rowSums(id[11:30]))
id$weeks <- as.numeric(rowSums(id[31:108]))
id$days <- as.numeric(rowSums(id[109:368]))

### load age data
ages <- readxl::read_excel("data_raw/Raw_ALERT_ElephantAges_Saunders211019.xlsx", sheet = 1)  # read in age data
names <- colnames(ages) ; names[8] <- 'comments' ; colnames(ages) <- names                    # rename column 8
ages <- ages[!is.na(ages$id_no),]      # remove empty rows -- 4 less ages than ids
table(ages$age_category)               # Some values incorrectly saved as dates not categories
ages$age_category <- ifelse(ages$age_category == '44228', '1-2',
                            ifelse(ages$age_category == '44257', '2-3',
                                   ifelse(ages$age_category == '44289', '3-4',
                                          ifelse(ages$age_category == '44320', '4-5',
                                                 ifelse(ages$age_category == '44352', '5-6',
                                                        ifelse(ages$age_category == '44383', '6-7',
                                                               ifelse(ages$age_category == '44415', '7-8',
                                                                      ifelse(ages$age_category == '44447', '8-9',
                                                                             ifelse(ages$age_category == '44478', '9-10',
                                                                                    ages$age_category)))))))))
ages$age_category <- ifelse(ages$age_category == '8-9 months', '0-3',   # group together calves <3 years
                            ifelse(ages$age_category == '<1', '0-3',
                                   ifelse(ages$age_category == '1', '0-3',
                                          ifelse(ages$age_category == '1-2', '0-3',
                                                 ifelse(ages$age_category == '2', '0-3',
                                                        ifelse(ages$age_category == '2-3', '0-3',
                                                               ages$age_category))))))
ages$age_category <- ifelse(ages$age_category == '20-39', '20-35',
                            ifelse(ages$age_category == '15-20', '15-19', ages$age_category))
table(ages$age_category)               # Now correct apart from missing values

### Add values for missing elephants -- Helen estimated ages from images sent by Bex, 11th November 2021
ages$id_no[which(ages$age_category == 'missing')]       # "M0026" "M0154" "M0196" "M0240"
ages$age_class[which(ages$age_category == 'missing')]   # calf, NA, adult, pubescent

ages$age_category[which(ages$id_no == 'M0026')] <- '1-2'
ages$age_category[which(ages$id_no == 'M0196')] <- '20-25'
ages$age_category[which(ages$id_no == 'M0240')] <- '10-15'

ages$id_no[which(is.na(ages$age_category))]      # "F0052"  "F0070" "F0137" "U0008" "M0138" "M0223"  "M0227"
ages$age_class[which(is.na(ages$age_category))]  #  "Adult" "Pubescent" NA  "Calf"  "Pubescent" "Adult"  "Adult"

ages$age_category[which(ages$id_no == 'F0052')] <- '35-50'
ages$age_category[which(ages$id_no == 'F0070')] <- '6-7'
ages$age_category[which(ages$id_no == 'M0223')] <- '20-25'
ages$age_category[which(ages$id_no == 'M0227')] <- '20-25'

### 6 elephants don't match ages and ID, or are only a number, no other data -- remove from ages spreadsheet.
ages <- ages[ages$id_no != 'F0137',]
ages <- ages[ages$id_no != 'M0013',]
ages <- ages[ages$id_no != 'M0116',]
ages <- ages[ages$id_no != 'M0125',]
ages <- ages[ages$id_no != 'M0154',]
ages <- ages[ages$id_no != 'M0207',]

### Correct mismatches for merging of spreadsheets
id$name[88]   <- 'Kaitlyn'                                 # name misspelled as 'Kaitlun' in ID sheet
id$family[11] <- 'U23 F69' ; ages$family[11] <- 'U23 F69'  # Margaret different families recorded in 2 spreadsheets

### Remove elephants with wrong name in ID/ages sheet -- confusion over names suggests potentially incorrect sightings data
id <- id[id$id_no != 'F0157',] ; ages <- ages[ages$id_no != 'F0157',]
id <- id[id$id_no != 'M0138',] ; ages <- ages[ages$id_no != 'M0138',]

### Add age data into id sheet
id <- merge(x = id, y = ages)             # 468 elephants
id <- id[,c(1:5,369,6:10,370,11:368)]

### Write out processed data
write_delim(id, "data_processed/motnp_id.csv", delim = ",", col_names = T)

### Clear Environment
rm(ages, days.numbers, id, id.raw, names, nkash, days, i, months, weeks)
#### Encounters ####
encounters <- readxl::read_excel("data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 3)
e <- encounters[c(1:3,8:27)]
e <- janitor::clean_names(e)
colnames(e)[c(4,5,7,8,23)] <- c('gps_s','gps_e','total_id','perc_id','present')
group.size <- paste('e',1:76, sep = '')

d <- e %>% separate(present,
                    as.character(c(1:max(e$total_id, na.rm = T))))
colnames(d)[23:98] <- group.size

for(i in 23:98) {
  for(j in 1:nrow(d)){
    d[j,i] <- ifelse(d[j,i] == 'B1', NA,
                     ifelse(d[j,i] == 'B2', NA,
                            ifelse(d[j,i] == 'B3', NA,
                                   ifelse(d[j,i] == 'B4', NA,
                                          ifelse(d[j,i] == 'B6', NA,
                                                 ifelse(d[j,i] == 'B7', NA,
                                                        ifelse(d[j,i] == '', NA, d[j,i])))))))
  }
}
d <- unite(data = d, col = present, 23:98, remove = F, sep = '_', na.rm=T)
e$present <- d$present
d <- e %>% separate(col = present, sep = '_',
                    into = as.character(group.size, na.rm = T))
colnames(d)[23:98] <- group.size
d$e1 <- ifelse(d$e1 == '', NA, d$e1)

# recount total_id and perc_id
d$count <- rowSums(!is.na(d[23:98]))
d$count_mismatch <- d$count-d$total_id
summary(d$count_mismatch) # range -7 to +8
d$total_elephants <- as.numeric(d$total_elephants)
d$perc_id_new <- 100*d$count / d$total_elephants
d <- d[c(1:7,99:101,8:98)]
colnames(d[7:11]) <- c('total_id_dy','total_id_hkm','total_id_diff','perc_id_dy','perc_id_hkm')

d <- d %>% 
  dplyr::rename(total_id_dy = total_id,
                total_id_hkm = count,
                total_id_diff = count_mismatch,
                perc_id_dy = perc_id,
                perc_id_hkm = perc_id_new)

# correct mis-recorded locations -- 3 which are an order of magnitude from the rest in either east or south direction, one for which the east value started 17 instead of 25 but otherwise did not match south so assume that changing the 25 to 17 was error in recording.
d$gps_e[232] <- d$gps_e[232]/10 # value was 25470333 -- 2547033 puts it within cluster of all other points
d$gps_e[480] <- d$gps_e[480]*10 # value was 254665 -- 2546650 puts it within cluster of all other points
d$gps_s[672] <- d$gps_s[672]*10 # value was 175296 -- 1752960 puts it within cluster of all other points
d$gps_e[624] <- d$gps_e[624]+800000 # value started 17 rather than 25 (1745269), but then matched format of other east coordinates -- 2545269 puts it within cluster of all other points

### write csv
write_delim(d, 'data_processed/motnp_encounters.csv', na = 'NA', col_names = T, delim = ',')

### clear environment, leave d for mapping
rm(e, encounters, group.size, i, j)
#### Map sightings ####
par(mfrow = c(1,1))

### determine positions of grid lines
min(d$gps_e, na.rm=T) ; max(d$gps_e, na.rm=T) # 2542980 and 2554583 -- grid lines from 2543000 to 2555000
min(d$gps_s, na.rm=T) ; max(d$gps_s, na.rm=T) # 1745598 and 1755287 -- grid lines from 1745000 to 1755000

### highlight corrected points to ensure they look sensible
plot(d$gps_e, d$gps_s, col = col.alpha(rangi2, 0.3), pch = 19, xaxt = 'n', yaxt = 'n',
     xlab = 'east', ylab = 'south',
     xlim = c(min(d$gps_e, na.rm=T),max(d$gps_e, na.rm=T)),
     ylim = c(max(d$gps_s, na.rm=T),min(d$gps_s, na.rm=T))) # upside down so more south is bottom
for(i in c(232,480,624,672)) points(d$gps_e[i], d$gps_s[i], col = 'red', pch = 19) # fall within common sighting locations now.
abline(v = seq(2543000,2555000,1000), col = col.alpha('black',0.3))
abline(h = seq(1745000,1755000,1000), col = col.alpha('black',0.3))
axis(1, at = seq(2543000,2555000,1000), las = 2,
     labels = as.character(seq(25.43,25.55,0.01)))
axis(2, at = seq(1745000,1755000,1000), las = 1,
     labels = as.character(seq(-17.45,-17.55,-0.01)))

### points based on herd type
plot(d$gps_e, d$gps_s, col = 'white', pch = 19, xaxt = 'n', yaxt = 'n',
     xlab = 'east', ylab = 'south',
     xlim = c(min(d$gps_e, na.rm=T),max(d$gps_e, na.rm=T)),
     ylim = c(max(d$gps_s, na.rm=T),min(d$gps_s, na.rm=T)))
for(i in 1:length(d$date)) points(d$gps_e[i], d$gps_s[i],
                                  pch = ifelse(d$type[i] == 'MO', 1,
                                               ifelse(d$type[i] == 'MX', 2,
                                                      ifelse(d$type[i] == 'BH', 4, 0))),
                                  col = ifelse(d$type[i] == 'MO','blue',
                                               ifelse(d$type[i] == 'MX','darkgreen',
                                                      ifelse(d$type[i] == 'BH','magenta','black'))))
abline(v = seq(2543000,2555000,1000), col = col.alpha('black',0.3))
abline(h = seq(1745000,1755000,1000), col = col.alpha('black',0.3))
axis(1, at = seq(2543000,2555000,1000), las = 2,
     labels = as.character(seq(25.43,25.55,0.01)))
axis(2, at = seq(1745000,1755000,1000), las = 1,
     labels = as.character(seq(-17.45,-17.55,-0.01)))

### males only
males <- rbind(d[d$type == 'MO',], d[d$type == 'MX',])
plot(gps_s ~ gps_e, data = males,
     col = col.alpha(rangi2,0.4), pch = 19, xaxt = 'n', yaxt = 'n',
     xlab = 'east', ylab = 'south',
     xlim = c(min(d$gps_e, na.rm=T),max(d$gps_e, na.rm=T)),
     ylim = c(max(d$gps_s, na.rm=T),min(d$gps_s, na.rm=T)))
abline(v = seq(2543000,2555000,1000), col = col.alpha('black',0.3))
abline(h = seq(1745000,1755000,1000), col = col.alpha('black',0.3))
axis(1, at = seq(2543000,2555000,1000), las = 2,
     labels = as.character(seq(25.43,25.55,0.01)))
axis(2, at = seq(1745000,1755000,1000), las = 1,
     labels = as.character(seq(-17.45,-17.55,-0.01)))

### clear environment
rm(d, males, i)

#### Create dyadic network ####
### first load data d (encounters)
d <- read_csv('data_processed/motnp_encounters.csv')
for(j in 2:ncol(d)){
  for(i in 1:nrow(d)){
    d[i,j] <- ifelse(d[i,j] == 'NA',NA,d[i,j])
  }
}
str(d) # d$total_elephants left as character due to presence of some estimate values (e.g. 50+). If need to convert this to numeric, will first need to remove additional symbols.

### consider only the first sighting of each group, not the repeated measures every 15 minutes
typ1 <- d[d$typ == 1,] ; typNA <- d[is.na(d$typ),]
first <- rbind(typNA, typ1) ; first <- first[!is.na(first$date),]
first <- first[,c(1,3:9,23:98)]   # remove columns counting number of each demographic group
first$encounter <- 1:nrow(first) ; first <- first[,c(85,1:84)] # index variable of sightings
first$date <- lubridate::as_date(first$date)

# make long, so every row is a sighting of an individual elephant
eles_long <- gather(data = first, key = number, value = elephant, e1:e76, factor_key = TRUE)
eles_long <- eles_long[!is.na(eles_long$elephant),]

# 2 incorrectly labelled elephants: no such code as N or D so N40 and D128 must be wrong. On keyboard, N next to M and D next to F --> assume that D128=F128 and N40=M40.
eles_long$elephant[which(eles_long$elephant == 'D128')] <- 'F128'  # row 2297
eles_long$elephant[which(eles_long$elephant == 'N40')] <- 'M40'    # row 1859

# write to csv for use in plotting
write_delim(eles_long, 'data_processed/motnp_eles_long.csv',col_names = T, delim = ',')

### make data frame that matches required structure for asnipe::get_associations_points_tw()
eles_asnipe <- eles_long[,c(1,3,4,5,11)]
eles_asnipe$location <- as.factor(paste(eles_asnipe$gps_s, '_', eles_asnipe$gps_e))
eles_asnipe <- eles_asnipe[,c(1,2,5,6)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$Date <- as.numeric(eles_asnipe$Date)
eles_asnipe$ID <- as.factor(eles_asnipe$ID)

eles_asnipe$Time <- ifelse(eles_asnipe$Time > 1, NA, eles_asnipe$Time) # time = proportion of day so anything >1 has to be wrong. 1 visit is recorded at time 1599 -- cannot be corrected, so convert to missing value.
eles_asnipe$Time <- eles_asnipe$Time*86400   # convert time values to seconds so actually readable
hist(eles_asnipe$Time, breaks = 30, las = 1) # can see the times we were out from the lunchtime reduction

# generate group-by-individual matrix
a <- get_associations_points_tw(eles_asnipe,            # data to input
                                time_window = 60,       # >1 minute apart = different groups
                                which_days = NULL,      # all days
                                which_locations = NULL) # all locations
gbi <- a[[1]] # group-by-individual = first part of a

n <- get_network(gbi, data_format = 'GBI', association_index = 'SRI') # 472x472 matrix containing SRI values for all dyads
n <- as.data.frame(n)
n$id1 <- row.names(n) ; n <- n[,c(473, 1:472)] # add column at start to indicate identity of first elephant in dyad

sna <- gather(data = n, key = id2, value = sri, 2:473, factor_key = TRUE) # long format data.frame containing all dyad pairs (all included twice because F1-F2 == F2-F1)
sna$sri <- ifelse(sna$id1 == sna$id2, NA, sna$sri)
sna$sex1 <- sna$id1 ; sna$sex2 <- sna$id2
sna <- separate(sna, sex1, into = c('sex1','number1'), sep = 1)
sna <- separate(sna, sex2, into = c('sex2','number2'), sep = 1)
sna <- sna[,c(1:4,6)]
sna <- sna[!is.na(sna$sri),]
sna$id2 <- as.character(sna$id2)

### remove duplicates
sna_no_duplicates <- sna %>% 
  rowwise() %>% 
  mutate(dyad = paste(sort(c(id1,id2)), collapse = '_')) %>% 
  ungroup() %>% 
  distinct_at(vars(dyad))

sna$dyad <- NA
for(i in 1:nrow(sna)){
  sna$dyad[i] <- paste(sort(c(sna$id1[i],sna$id2[i])), collapse = '_')
}

sna_no_duplicates <- separate(sna_no_duplicates, dyad, into = c('id1','id2'), sep = '_', remove = F)
sna_no_duplicates <- separate(sna_no_duplicates, id1, into = c('sex1','number1'), sep = 1, remove = F)
sna_no_duplicates <- separate(sna_no_duplicates, id2, into = c('sex2','number2'), sep = 1, remove = F)
sna_no_duplicates <- sna_no_duplicates[,c(1,4,3,2,5,6)]

sna_no_duplicates$sri <- NA
for(i in 1:nrow(sna_no_duplicates)) {
  sna_no_duplicates$sri[i] <- sna$sri[match(x = sna_no_duplicates$dyad[i], table = sna$dyad)]
}

write_delim(sna_no_duplicates, "data_processed/motnp_dyads_noduplicates.csv",
            delim = ",", na = "", append = FALSE, col_names = T, quote_escape = "double")

associations <- sna[sna$sri > 0,] # considering only dyads for which there was at least some association
associations <- associations[!is.na(associations$sri),]

par(mfrow = c(2,1))
hist(sna$sri,          xlab = 'SRI', las = 1, ylab = '', main = 'Including all dyads')
hist(associations$sri, xlab = 'SRI', las = 1, ylab = '', main = 'After removing unassociated dyads')

### clear environment, reset plot window
rm(a,associations,d,eles_asnipe,eles_long,first,gbi,n,sna,sna_no_duplicates,typ1,typNA,i,j)
par(mfrow = c(1,1))

#### Create plotting data -- very slow so do not repeat once files are created ####
### notes data frame
ele_nodes <- read_delim('data_processed/motnp_id_22.01.13.csv', delim = ',')
ele_nodes <- ele_nodes[,c(1:7,11:12)]

### edges data frame
ele_links <- read.csv("data_processed/motnp_6thJan2022/motnp_dyads_22.01.06.csv", header=T, as.is=T)
ele_links$type <- paste(ele_links$sex1,ele_links$sex2, sep = '')
ele_links$type <- ifelse(ele_links$type == 'FM', 'MF',
                         ifelse(ele_links$type == 'FU','UF',
                                ifelse(ele_links$type == 'MU', 'UM', ele_links$type)))
ele_links$num1 <- ele_links$id1 ; ele_links$num2 <- ele_links$id2
ele_links <- separate(ele_links, num1, into = c('s1','num1'), sep = 1)
ele_links <- separate(ele_links, num2, into = c('s2','num2'), sep = 1)
ele_links$num1 <- sprintf("%04d", as.integer(ele_links$num1))
ele_links$num2 <- sprintf("%04d", as.integer(ele_links$num2))
ele_links$from <- paste(ele_links$s1,ele_links$num1,sep='')
ele_links$to <- paste(ele_links$s2,ele_links$num2,sep='')
ele_links <- ele_links[,c(12,13,7,3:5,1,2,6)]
colnames(ele_links)[4] <- 'weight'

### sightings per elephant
eles_long <- read_delim('data_processed/motnp_eles_long_22.01.13.csv', delim = ',')
eles_long$time <- as.numeric(eles_long$time)
eles_long$id_no <- eles_long$elephant
eles_long <- separate(eles_long, id_no, into = c('sex','num'), sep = 1)
eles_long$num <- sprintf("%04d", as.integer(eles_long$num))
eles_long$id_no <- paste(eles_long$sex, eles_long$num, sep='')
eles_long <- eles_long[,c(1:9,11,12,14)]

### make column for shortened id number (not padded with zeroes) for use as labels
ele_nodes$id  <- ele_nodes$id_no
ele_nodes     <- separate(ele_nodes, id, into = c('sex','num'), sep = 1)
ele_nodes$num <- as.numeric(ele_nodes$num)
ele_nodes$id  <- paste(ele_nodes$sex,ele_nodes$num,sep='')

### cleaning
ele_links <- ele_links[which(ele_links$from != 'F0157'),] # remove F0157 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'F0157'),] # remove F0157 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'F0158'),] # remove F0158 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'F0158'),] # remove F0158 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'F0176'),] # remove F0176 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'F0176'),] # remove F0176 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'M0013'),] # remove M0013 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'M0013'),] # remove M0013 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'M0125'),] # remove M0125 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'M0125'),] # remove M0125 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'M0138'),] # remove M0138 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'M0138'),] # remove M0138 from links -- doesn't exist in IDs
ele_nodes <- ele_nodes[which(ele_nodes$id_no!= 'M0175'),] # remove M0175 from nodes -- doesn't exist in sightings

length(unique(ele_nodes$id_no)) ; length(unique(ele_links$from)) ; length(unique(ele_links$to)) # 466, 465, 465 (One less for from and to because F1 is only in "from" and U9 is only in "to". No elephant is from and to itself.)

### Insert column for count of sightings
ele_nodes$count <- NA
for(i in 1:nrow(ele_nodes)){
  ele_nodes$count[i] <- sum(eles_long$id_no == ele_nodes$id_no[i])
}

### Add information about individuals in each dyad
ele_links$ageclass1 <- NA ; ele_links$ageclass2 <- NA
ele_links$age1   <- NA ; ele_links$age2   <- NA
ele_links$days1  <- NA ; ele_links$days2  <- NA
ele_links$count1 <- NA ; ele_links$count2 <- NA
# WARNING -- THIS STEP TAKES ABOUT 30 MINUTES TO RUN
for(i in 1:nrow(ele_links)){
  for(j in 1:nrow(ele_nodes)){
    if(ele_nodes$id_no[j] == ele_links$from[i]) {ele_links$ageclass1[i] <- ele_nodes$age_class[j]}
    if(ele_nodes$id_no[j] == ele_links$to[i])   {ele_links$ageclass2[i] <- ele_nodes$age_class[j]}
    if(ele_nodes$id_no[j] == ele_links$from[i]) {ele_links$age1[i]      <- ele_nodes$age_category[j]}
    if(ele_nodes$id_no[j] == ele_links$to[i])   {ele_links$age2[i]      <- ele_nodes$age_category[j]}
    if(ele_nodes$id_no[j] == ele_links$from[i]) {ele_links$days1[i]     <- ele_nodes$days[j]}
    if(ele_nodes$id_no[j] == ele_links$to[i])   {ele_links$days2[i]     <- ele_nodes$days[j]}
    if(ele_nodes$id_no[j] == ele_links$from[i]) {ele_links$count1[i]    <- ele_nodes$count[j]}
    if(ele_nodes$id_no[j] == ele_links$to[i])   {ele_links$count2[i]    <- ele_nodes$count[j]}
  }
}

ele_links$age.dyad <- ifelse(ele_links$ageclass1 == 'Adult',
                             ifelse(ele_links$ageclass2 == 'Adult', 'AA',
                                    ifelse(ele_links$ageclass2 == 'Pubescent', 'AP',
                                           ifelse(ele_links$ageclass2 == 'Juvenile', 'AJ', 'AC'))),
                             ifelse(ele_links$ageclass1 == 'Pubescent',
                                    ifelse(ele_links$ageclass2 == 'Adult', 'AP',
                                           ifelse(ele_links$ageclass2 == 'Pubescent', 'PP',
                                                  ifelse(ele_links$ageclass2 == 'Juvenile', 'PJ', 'PC'))),
                                    ifelse(ele_links$ageclass1 == 'Juvenile',
                                           ifelse(ele_links$ageclass2 == 'Adult', 'AJ',
                                                  ifelse(ele_links$ageclass2 == 'Pubescent', 'PJ',
                                                         ifelse(ele_links$ageclass2 == 'Juvenile', 'JJ', 'CJ'))), 'CC')))

### create column for combined age and sex in both data frames
ele_nodes$dem_class <- paste(sep = '',
                             ifelse(ele_nodes$age_class == 'Adult','A',
                                    ifelse(ele_nodes$age_class == 'Pubescent','P',
                                           ifelse(ele_nodes$age_class == 'Juvenile', 'J','C'))),
                             ele_nodes$sex)
ele_links$dem_class1 <- paste(sep = '',
                              ifelse(ele_links$ageclass1 == 'Adult','A',
                                     ifelse(ele_links$ageclass1 == 'Pubescent','P',
                                            ifelse(ele_links$ageclass1 == 'Juvenile','J','C'))),
                              ele_links$sex1)
ele_links$dem_class2 <- paste(sep = '',
                              ifelse(ele_links$ageclass2 == 'Adult','A',
                                     ifelse(ele_links$ageclass2 == 'Pubescent','P',
                                            ifelse(ele_links$ageclass2 == 'Juvenile','J','C'))),
                              ele_links$sex2)

### write CSV files of both data frames so this section doesn't need to be run again
write_delim(ele_links, path = 'data_processed/motnp_elelinks.csv', delim = ',', col_names = T)
write_delim(ele_nodes, path = 'data_processed/motnp_elenodes.csv', delim = ',', col_names = T)

### clear environment
rm(ele_links, ele_nodes, eles_long, i, j)
dev.off()

#### Plotting Networks ####
### load data
ele_links <- read_csv('data_processed/motnp_elelinks.csv')
ele_links$age.dyad   <- as.factor(ele_links$age.dyad)
ele_links$sex1       <- as.factor(ele_links$sex1)
ele_links$sex2       <- as.factor(ele_links$sex2)
ele_links$dem_class1 <- as.factor(ele_links$dem_class1)
ele_links$dem_class2 <- as.factor(ele_links$dem_class2)
str(ele_links)

ele_nodes <- read_csv('data_processed/motnp_elenodes.csv')
ele_nodes$sex       <- as.factor(ele_nodes$sex)
ele_nodes$age_class <- as.factor(ele_nodes$age_class)
ele_nodes$dem_class <- as.factor(ele_nodes$dem_class)
str(ele_nodes)

### all elephants
e_edges <- ele_links[,c(1:4)] ; e_edges <- e_edges[!is.na(e_edges$from),]
e_nodes <- ele_nodes ; e_nodes <- e_nodes[!is.na(e_nodes$id_no),]
e.net <- igraph::graph_from_data_frame(d = e_edges,         # data frame of network edges
                                       vertices = e_nodes,  # node ID and attributes of nodes
                                       directed = F)        # direction of edge strengths
graph_attr(e.net, 'layout') <- layout_with_lgl
# plot all edges
set.seed(11)
plot(e.net, edge.width = e_edges$weight,
     #vertex.shape = ifelse(e_nodes$sex == 'Female', 'square',
     #                     ifelse(e_nodes$sex == 'Male','circle','raster')),
     vertex.color= ifelse(e_nodes$age_class == 'Adult','seagreen1',
                          ifelse(e_nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(e_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e_nodes$sex,
     edge.curved = 0,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
     vertex.label.cex = 0.5,
     vertex.label.family = 'Helvetica')
legend(x = -1.2, y = 1, legend = c('Adult','Pubescent','Juvenile','Calf'), # takes several seconds to appear
       pch = 19, col = c('seagreen1','skyblue','yellow','magenta'))
# plot only strongest edges
summary(e_edges$weight) # 3rd quartile = 0
quantile(e_edges$weight, c(0.9, 0.95))  # 90th percentile = 0.119, 95th percentile = 0.179
set.seed(11)
plot(e.net,
     edge.width = ifelse(e_edges$weight < 0.12, 0, e_edges$weight),
     edge.color = ifelse(e_edges$weight < 0.12, 0, rgb(0.1,0.1,0.1,alpha=0.5)),
     edge.curved = 0,
     vertex.color= ifelse(e_nodes$age_class == 'Adult','seagreen1',
                          ifelse(e_nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(e_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e_nodes$sex,
     vertex.label.cex = 0.5,
     vertex.label.family = 'Helvetica')
legend(x = -1.2, y = 1, legend = c('Adult','Pubescent','Juvenile','Calf'), # takes several seconds to appear
       pch = 19, col = c('seagreen1','skyblue','yellow','magenta'))

## all elephants seen 10 or more times
e10_n <- ele_nodes[ele_nodes$count >9,]
e10_l <- ele_links[ele_links$count1 >9,]
e10_l <- e10_l[e10_l$count2 >9,]
e10.net <- igraph::graph_from_data_frame(d = e10_l, vertices = e10_n, directed = F)
graph_attr(e10.net, 'layout') <- layout_in_circle
# plot all edges
plot(e10.net, edge.width = e10_l$weight, edge.curved = 0,
     vertex.color= ifelse(e10_n$age_class == 'Adult','seagreen1',
                          ifelse(e10_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e10_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 4,
     vertex.label = e10_n$sex,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.2))
# plot only strongest edges
summary(e10_l$weight) # median = 0.0617, 3rd quartile = 0.122
quantile(e10_l$weight, c(0.9, 0.95))  # 90th percentile = 0.192, 95th percentile = 0.240
plot(e10.net, edge.curved = 0,
     edge.width = ifelse(e10_l$weight < 0.2, 0, e10_l$weight),
     edge.color = ifelse(e10_l$weight < 0.2, 0, rgb(0.1,0.1,0.1,alpha=0.4)),
     vertex.color= ifelse(e10_n$age_class == 'Adult','seagreen1',
                          ifelse(e10_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e10_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 4,
     vertex.label = e10_n$sex)

# compare edge strengths between sexes
AP <- e10_l[e10_l$age.dyad == 'AA' | e10_l$age.dyad == 'AP' | e10_l$age.dyad == 'PP',]
MM <- AP[AP$type == 'MM',]
summary(MM$weight)
FF <- AP[AP$type == 'FF',]
summary(FF$weight)
MF <- AP[AP$type == 'MF',]
summary(MF$weight)
# plot different types of interaction to see how males and females compare
par(mfrow = c(3,1))
hist(MM$weight, las = 1, breaks = 30, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Male-Male Interactions')
hist(MF$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Male-Female Interactions')
hist(FF$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Female-Female Interactions')
# plot again, but cut off top of 0 bar just so as to see the rest of it better
hist(MM$weight, las = 1, breaks = 30, xlab = 'weight', xlim = c(0,1), ylim = c(0,300),
     main = 'Adult and Pubescent Male-Male Interactions')
hist(MF$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1), ylim = c(0,400),
     main = 'Adult and Pubescent Male-Female Interactions')
hist(FF$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1), ylim = c(0,150),
     main = 'Adult and Pubescent Female-Female Interactions')

# compare edge strengths between ages
AA <- e10_l[e10_l$age.dyad == 'AA',]
summary(AA$weight)
AP <- e10_l[e10_l$age.dyad == 'AP',]
summary(AP$weight)
PP <- e10_l[e10_l$age.dyad == 'PP',]
summary(PP$weight)
# plot different types of interaction to see how males and females compare
hist(AA$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Male and Female Adult-Adult Interactions')
hist(AP$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Male and Female Adult-Pubescent Interactions')
hist(PP$weight, las = 1, breaks = 30, xlab = 'weight', xlim = c(0,1),
     main = 'Male and Female Pubescent-Pubescent Interactions')
# plot again, but cut off top of 0 bar just so as to see the rest of it better
hist(AA$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1), ylim = c(0,500),
     main = 'Male and Female Adult-Adult Interactions')
hist(AP$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1), ylim = c(0,300),
     main = 'Male and Female Adult-Pubescent Interactions')
hist(PP$weight, las = 1, breaks = 30, xlab = 'weight', xlim = c(0,1), ylim = c(0,50),
     main = 'Male and Female Pubescent-Pubescent Interactions')

## all elephants seen 20 or more times
e20_n <- ele_nodes[ele_nodes$count >19,]
e20_l <- ele_links[ele_links$count1 >19,]
e20_l <- e20_l[e20_l$count2 >19,]
e20.net <- igraph::graph_from_data_frame(d = e20_l, vertices = e20_n, directed = F)
graph_attr(e20.net, 'layout') <- layout_in_circle
# plot all edges
plot(e20.net, edge.width = e20_l$weight, edge.curved = 0,
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$dem_class,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.4))
plot(e20.net, edge.width = e20_l$weight, edge.curved = 0,
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$dem_class,
     edge.color = ifelse(e20_l$dem_class1 == 'AM' & e20_l$dem_class2 == 'AM' |
                           e20_l$dem_class2 == 'AM' & e20_l$dem_class1 == 'PM' |
                           e20_l$dem_class1 == 'PM' & e20_l$dem_class2 == 'AM' |
                           e20_l$dem_class2 == 'PM' & e20_l$dem_class1 == 'PM',
                         rgb(1.0,0.0,0.0,alpha=0.4),
                         rgb(0.1,0.1,0.1,alpha=0.4)))
# plot only strongest edges
summary(e20_l$weight) # median = 0.073, 3rd quartile = 0.127
quantile(e20_l$weight, c(0.9,0.95)) # 90th = 0.200, 95th = 0.248
plot(e20.net, edge.curved = 0,
     edge.width = ifelse(e20_l$weight < 0.2, 0, e20_l$weight),
     edge.color = ifelse(e20_l$weight < 0.2, 0, rgb(0.1,0.1,0.1,alpha=0.4)),
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$dem_class)

# compare edge strengths between sexes
AP <- e20_l[e20_l$age.dyad == 'AA' | e20_l$age.dyad == 'AP' | e20_l$age.dyad == 'PP',]
MM <- AP[AP$type == 'MM',]
summary(MM$weight)
FF <- AP[AP$type == 'FF',]
summary(FF$weight)
MF <- AP[AP$type == 'MF',]
summary(MF$weight)
# plot different types of interaction to see how males and females compare
hist(MM$weight, las = 1, breaks = 30, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Male-Male Interactions')
hist(MF$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Male-Female Interactions')
hist(FF$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Female-Female Interactions')
# plot again, but cut off top of 0 bar just so as to see the rest of it better
hist(MM$weight, las = 1, breaks = 30, xlab = 'weight', xlim = c(0,1), ylim = c(0,80),
     main = 'Adult and Pubescent Male-Male Interactions')
hist(MF$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1), ylim = c(0,120),
     main = 'Adult and Pubescent Male-Female Interactions')
hist(FF$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1), ylim = c(0,40),
     main = 'Adult and Pubescent Female-Female Interactions')

# compare edge strengths between ages
AA <- e20_l[e20_l$age.dyad == 'AA',]
summary(AA$weight)
AP <- e20_l[e20_l$age.dyad == 'AP',]
summary(AP$weight)
PP <- e20_l[e20_l$age.dyad == 'PP',]
summary(PP$weight)
# plot different types of interaction to see how males and females compare
hist(AA$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Male and Female Adult-Adult Interactions')
hist(AP$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Male and Female Adult-Pubescent Interactions')
hist(PP$weight, las = 1, breaks = 30, xlab = 'weight', xlim = c(0,1),
     main = 'Male and Female Pubescent-Pubescent Interactions')
# plot again, but cut off top of 0 bar just so as to see the rest of it better
hist(AA$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1), ylim = c(0,150),
     main = 'Male and Female Adult-Adult Interactions')
hist(AP$weight, las = 1, breaks = 50, xlab = 'weight', xlim = c(0,1), ylim = c(0,80),
     main = 'Male and Female Adult-Pubescent Interactions')
hist(PP$weight, las = 1, breaks = 30, xlab = 'weight', xlim = c(0,1), ylim = c(0,10),
     main = 'Male and Female Pubescent-Pubescent Interactions')

### males only -- warning, slow to plot
males <- ele_nodes[ele_nodes$sex=='M'&ele_nodes$age_class!='Calf'&ele_nodes$age_class!='Juvenile',]
males <- males[!is.na(males$id_no),]

m_links <- ele_links[ele_links$type == 'MM',]
m_links <- m_links[m_links$age1 == 'Adult' | m_links$age1 == 'Pubescent',]
m_links <- m_links[m_links$age2 == 'Adult' | m_links$age2 == 'Pubescent',]

m_edges <- m_links[,c(1:4)] ; m_edges <- m_edges[!is.na(m_edges$from),]
m_nodes <- males ; m_nodes <- m_nodes[!is.na(m_nodes$id_no),]

m.net <- igraph::graph_from_data_frame(d = m_edges,         # data frame of network edges
                                       vertices = m_nodes,  # node ID and attributes of nodes
                                       directed = F)        # direction of edge strengths
graph_attr(m.net, 'layout') <- layout_with_lgl
# plot all edges
set.seed(11)
plot(m.net, edge.width = m_edges$weight,
     vertex.color = ifelse(m_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.shape = ifelse(m_nodes$count < 5, 'square', 'circle'),
     vertex.label = m_nodes$id,
     vertex.size = ifelse(m_nodes$count < 5, 3, 8),
     edge.curved = 0,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.3))
# plot only strongest edges
summary(m_edges$weight) # 3rd quartile = 0
quantile(m_edges$weight, c(0.9,0.95)) # 90th = 0.117, 95th = 0.178
set.seed(11)
plot(m.net,
     edge.width = ifelse(m_edges$weight < 0.18, 0, m_edges$weight),
     edge.color = ifelse(m_edges$weight < 0.19, 0, rgb(0.1,0.1,0.1,alpha=0.4)),
     edge.curved = 0,
     vertex.color = ifelse(m_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.shape = ifelse(m_nodes$count < 5, 'square', 'circle'),
     vertex.label = m_nodes$id,
     vertex.size = ifelse(m_nodes$count < 5, 3, 8))

### males with 10 or more sightings only -- quick to plot
males10 <- males[males$count > 9,]
males10 <- males10[!is.na(males10$id_no),]

m10_links <- m_links[m_links$count1 > 9,]
m10_links <- m10_links[m10_links$count2 > 9,]

m10_edges <- m10_links[,c(1:4)] ; m10_edges <- m10_edges[!is.na(m10_edges$from),]
m10_nodes <- males10 ; m10_nodes <- m10_nodes[!is.na(m10_nodes$id_no),]

summary(m10_edges$weight) # median = 0.688, 3rd quartile = 0.127
quantile(m10_edges$weight, c(0.9,0.95)) # 90th = 0.193, 95th = 0.236

AA <- m10_links[m10_links$age.dyad == 'AA',]
AP <- m10_links[m10_links$age.dyad == 'AP',]
PP <- m10_links[m10_links$age.dyad == 'PP',]
summary(AA$weight) # median = 0.076, 3rd quartile = 0.131
summary(AP$weight) # median = 0.061, 3rd quartile = 0.121
summary(PP$weight) # median = 0.030, 3rd quartile = 0.099

m10.net <- igraph::graph_from_data_frame(d = m10_edges,         # data frame of network edges
                                         vertices = m10_nodes,  # node ID and attributes of nodes
                                         directed = F)          # direction of edge strengths
graph_attr(m10.net, 'layout') <- layout_in_circle
# plot all edges
plot(m10.net, edge.width = m10_edges$weight,
     vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.label = m10_nodes$id,
     vertex.size = m10_nodes$count/2,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
     edge.curved = 0)
graph_attr(m10.net, 'layout') <- layout_with_lgl
set.seed(11)
plot(m10.net, edge.width = m10_edges$weight,
     vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.label = m10_nodes$id,
     vertex.size = m10_nodes$count/2,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
     edge.curved = 0)
# plot only strongest edges
set.seed(11)
plot(m10.net, 
     edge.width = ifelse(m10_edges$weight < 0.19, 0, m10_edges$weight),
     edge.color = ifelse(m10_edges$weight < 0.19, 0, rgb(0.1,0.1,0.1,alpha=0.6)),
     vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.label = m10_nodes$id,
     vertex.size = m10_nodes$count/2,
     edge.curved = 0)

### females -- compare strength of ties to males, slow to plot
females <- anti_join(ele_nodes, males)
females <- females[!is.na(females$id_no),]

f_links <- ele_links[ele_links$dem_class1 != 'AM' & ele_links$dem_class2 != 'AM',]
f_links <- f_links[f_links$dem_class1 != 'PM' & f_links$dem_class2 != 'PM',]
f_links <- f_links[f_links$from != 'F0157' | f_links$to != 'F0157',]

f_edges <- f_links[,c(1:4)] ; f_edges <- f_edges[!is.na(f_edges$from),]
f_nodes <- females[females$id_no != 'F0157',] ; f_nodes <- f_nodes[!is.na(f_nodes$id_no),]

f.net <- igraph::graph_from_data_frame(d = f_edges,         # data frame of network edges
                                       vertices = f_nodes,  # node ID and attributes of nodes
                                       directed = F)        # direction of edge strengths
graph_attr(f.net, 'layout') <- layout_with_lgl
# plot all edges
set.seed(11)
plot(f.net, edge.width = f_edges$weight,
     vertex.color = ifelse(f_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f_nodes$id,
     vertex.size = f_nodes$count^0.8, # note the power here -- very large variation in count
     edge.curved = 0,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5))
# plot only strongest edges
summary(f_edges$weight) # 3rd quartile = 0
quantile(f_edges$weight, c(0.9, 0.95)) # 90th = 0.129, 95th = 0.199
set.seed(11)
plot(f.net,
     edge.width = ifelse(f_edges$weight < 0.13, 0, f_edges$weight),
     edge.color = ifelse(f_edges$weight < 0.13, 0, rgb(0.1,0.1,0.1,alpha=0.4)),
     vertex.color = ifelse(f_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f_nodes$id,
     vertex.size = f_nodes$count^0.8, # note the power here -- very large variation in count
     edge.curved = 0)

### females with only 10 sightings -- quick to plot
females10 <- females[females$count > 9,]
females10 <- females10[!is.na(females10$id_no),]

f10_links <- f_links[f_links$count1 > 9,]
f10_links <- f10_links[f10_links$count2 > 9,]

f10_edges <- f10_links[,c(1:4)] ; f10_edges <- f10_edges[!is.na(f10_edges$from),]
f10_nodes <- females10 ; f10_nodes <- f10_nodes[!is.na(f10_nodes$id_no),]

f10.net <- igraph::graph_from_data_frame(d = f10_edges,         # data frame of network edges
                                         vertices = f10_nodes,  # node ID and attributes of nodes
                                         directed = F)          # direction of edge strengths
graph_attr(f10.net, 'layout') <- layout_with_lgl
# plot all edges
set.seed(11)
plot(f10.net, edge.width = f10_edges$weight,
     vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f10_nodes$id_no,
     vertex.size = f10_nodes$count/1.5,
     edge.curved = 0,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5))
graph_attr(f10.net, 'layout') <- layout_in_circle
set.seed(11)
plot(f10.net, edge.width = f10_edges$weight,
     vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f10_nodes$id,
     vertex.size = f10_nodes$count/2,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
     edge.curved = 0)

# plot only strongest edges
summary(f10_edges$weight) # median = 0.059, 3rd quartile = 0.13
quantile(f10_edges$weight, c(0.9, 0.95)) # 90th = 0.222, 95th = 0.320
set.seed(11)
plot(f10.net,
     edge.width = ifelse(f10_edges$weight < 0.2, 0, f10_edges$weight),
     edge.color = ifelse(f10_edges$weight < 0.2, 0, rgb(0.1,0.1,0.1,alpha=0.5)),
     vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f10_nodes$id,
     vertex.size = f10_nodes$count/2,
     edge.curved = 0)

### comparative plots male vs female
# plot all edges
set.seed(11) ; plot(m10.net,
                    edge.width = m10_edges$weight,
                    edge.color = rgb(0.1,0.1,0.1,alpha=20/nrow(m10_nodes)),
                    edge.curved = 0,
                    vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
                    vertex.label = m10_nodes$id,
                    vertex.size = m10_nodes$count/2)
set.seed(11) ; plot(f10.net,
                    edge.width = f10_edges$weight,
                    edge.color = rgb(0.1,0.1,0.1,alpha=20/nrow(f10_nodes)),
                    edge.curved = 0,
                    vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                                          ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                                 ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
                    vertex.label = f10_nodes$id,
                    vertex.size = f10_nodes$count/2)
# plot only strongest edges
quantile(m10_edges$weight, c(0.75, 0.9)) # 0.127, 0.192
quantile(f10_edges$weight, c(0.75, 0.9)) # 0.131, 0.222
set.seed(11) ; plot(m10.net,
                    edge.width = ifelse(m10_edges$weight < 0.2, 0,
                                        m10_edges$weight),
                    edge.color = ifelse(m10_edges$weight < 0.2, 0,
                                        rgb(0.1,0.1,0.1,alpha=50/nrow(m10_nodes))),
                    edge.curved = 0,
                    vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
                    vertex.label = m10_nodes$id,
                    vertex.size = m10_nodes$count/2)
set.seed(11) ; plot(f10.net,
                    edge.width = ifelse(f10_edges$weight < 0.2, 0,
                                        f10_edges$weight),
                    edge.color = ifelse(f10_edges$weight < 0.2, 0,
                                        rgb(0.1,0.1,0.1,alpha=50/nrow(f10_nodes))),
                    edge.curved = 0,
                    vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                                          ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                                 ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
                    vertex.label = f10_nodes$id,
                    vertex.size = f10_nodes$count/2)

### B7
b7.n <- ele_nodes[ele_nodes$herd == "B7",] ; b7.n <- b7.n[!is.na(b7.n$id_no),]
b7.e <- ele_links[ele_links$from == 'F0052' | ele_links$from == 'F0060' | ele_links$from == 'F0098' | 
                    ele_links$from == 'U0017' | ele_links$from == 'U0021',]
b7.e <- b7.e[b7.e$to == 'F0052' | b7.e$to == 'F0060' | b7.e$to == 'F0098' | 
               b7.e$to == 'U0017' | b7.e$to == 'U0021',]
b7 <- graph_from_data_frame(d = b7.e, vertices = b7.n, directed = F)
set.seed(5)
plot(b7, edge.curved = 0, vertex.label = b7.n$id,
     edge.width = b7.e$weight*8)

#### Comparing SRI ####
### load data
ele_links <- read_csv('data_processed/motnp_elelinks.csv')
ele_links$age.dyad   <- as.factor(ele_links$age.dyad)
ele_links$sex1       <- as.factor(ele_links$sex1)
ele_links$sex2       <- as.factor(ele_links$sex2)
ele_links$dem_class1 <- as.factor(ele_links$dem_class1)
ele_links$dem_class2 <- as.factor(ele_links$dem_class2)
str(ele_links)

ele_links$dem_dyad <- paste(ele_links$dem_class1, ele_links$dem_class2, sep = '_')
sort(unique(ele_links$dem_dyad))
# put adults before pubescents, juveniles and calves
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'PF_AF', 'AF_PF',
                             ifelse(ele_links$dem_dyad == 'PF_AM', 'AM_PF',ele_links$dem_dyad))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'PM_AM', 'AM_PM',
                             ifelse(ele_links$dem_dyad == 'PM_AF', 'AF_PM',ele_links$dem_dyad))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'PU_AM', 'AM_PU',
                             ifelse(ele_links$dem_dyad == 'PU_AF', 'AF_PU',ele_links$dem_dyad))

ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'JF_AF', 'AF_JF',
                             ifelse(ele_links$dem_dyad == 'JF_AM', 'AM_JF',ele_links$dem_dyad))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'JM_AM', 'AM_JM',
                             ifelse(ele_links$dem_dyad == 'JM_AF', 'AF_JM',ele_links$dem_dyad))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'JU_AM', 'AM_JU',
                             ifelse(ele_links$dem_dyad == 'JU_AF', 'AF_JU',ele_links$dem_dyad))

ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CF_AF', 'AF_CF',
                             ifelse(ele_links$dem_dyad == 'CF_AM', 'AM_CF',ele_links$dem_dyad))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CM_AM', 'AM_CM',
                             ifelse(ele_links$dem_dyad == 'CM_AF', 'AF_CM',ele_links$dem_dyad))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CU_AM', 'AM_CU',
                             ifelse(ele_links$dem_dyad == 'CU_AF', 'AF_CU',ele_links$dem_dyad))

# put pubescents before juveniles and calves
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'JF_PM', 'PM_JF',
                             ifelse(ele_links$dem_dyad == 'JM_PM', 'PM_JM',
                                    ifelse(ele_links$dem_dyad == 'JU_PM', 'PM_JU',ele_links$dem_dyad)))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'JF_PF', 'PF_JF',
                             ifelse(ele_links$dem_dyad == 'JM_PF', 'PF_JM',
                                    ifelse(ele_links$dem_dyad == 'JU_PF', 'PF_JU',ele_links$dem_dyad)))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'JF_PU', 'PU_JF',
                             ifelse(ele_links$dem_dyad == 'JM_PU', 'PU_JM',
                                    ifelse(ele_links$dem_dyad == 'JU_PU', 'PU_JU',ele_links$dem_dyad)))

ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CF_PM', 'PM_CF',
                             ifelse(ele_links$dem_dyad == 'CM_PM', 'PM_CM',
                                    ifelse(ele_links$dem_dyad == 'CU_PM', 'PM_CU',ele_links$dem_dyad)))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CF_PF', 'PF_CF',
                             ifelse(ele_links$dem_dyad == 'CM_PF', 'PF_CM',
                                    ifelse(ele_links$dem_dyad == 'CU_PF', 'PF_CU',ele_links$dem_dyad)))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CF_PU', 'PU_CF',
                             ifelse(ele_links$dem_dyad == 'CM_PU', 'PU_CM',
                                    ifelse(ele_links$dem_dyad == 'CU_PU', 'PU_CU',ele_links$dem_dyad)))

# put juveniles before calves
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CF_JF', 'JF_CF',
                             ifelse(ele_links$dem_dyad == 'CF_JM', 'JM_CF',
                                    ifelse(ele_links$dem_dyad == 'CF_JU', 'JU_CF',ele_links$dem_dyad)))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CM_JF', 'JF_CM',
                             ifelse(ele_links$dem_dyad == 'CM_JM', 'JM_CM',
                                    ifelse(ele_links$dem_dyad == 'CM_JU', 'JU_CM',ele_links$dem_dyad)))
ele_links$dem_dyad <- ifelse(ele_links$dem_dyad == 'CU_JF', 'JF_CU',
                             ifelse(ele_links$dem_dyad == 'CU_JM', 'JM_CU',
                                    ifelse(ele_links$dem_dyad == 'CU_JU', 'JU_CU',ele_links$dem_dyad)))

#ele_nodes <- read_csv('data_processed/motnp_elenodes.csv')
#ele_nodes$sex       <- as.factor(ele_nodes$sex)
#ele_nodes$age_class <- as.factor(ele_nodes$age_class)
#ele_nodes$dem_class <- as.factor(ele_nodes$dem_class)
#str(ele_nodes)

associations <- ele_links[ele_links$weight > 0,]

links10 <- ele_links[ele_links$count1 > 9 & ele_links$count2 > 9,]
assoc10 <- associations[associations$count1 > 9 & associations$count2 > 9,]

### Compare types
boxplot(weight ~ type, data = ele_links, xlab = 'dyad sex', las = 1,
        main = 'association strength including all dyads')
boxplot(weight ~ type, data = associations, xlab = 'dyad sex', las = 1, notch = T,
        main = 'association strength including all dyads \nobserved together at least once')
boxplot(weight ~ type, data = links10, xlab = 'dyad sex', las = 1, notch = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times')
boxplot(weight ~ type, data = assoc10, xlab = 'dyad sex', las = 1, notch = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times and \nobserved together at least once')






### Compare age dyads
boxplot(weight ~ age.dyad, data = ele_links, xlab = 'dyad age', las = 1, notch = T,
        main = 'association strength including all dyads')
boxplot(weight ~ age.dyad, data = associations, xlab = 'dyad age', las = 1, notch = T,
        main = 'association strength including all dyads\nobserved together at least once')
boxplot(weight ~ age.dyad, data = links10, xlab = 'dyad sex', las = 1, notch = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times')
boxplot(weight ~ age.dyad, data = assoc10, xlab = 'dyad sex', las = 1, notch = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times and \nobserved together at least once')






### Compare full demography
boxplot(weight ~ dem_dyad, data = ele_links, xlab = 'dyad type', las = 1, horizontal = T,
        main = 'association strength including all dyads')
boxplot(weight ~ dem_dyad, data = associations, xlab = 'dyad age', las = 1, horizontal = T,
        main = 'association strength including all dyads\nobserved together at least once')
boxplot(weight ~ dem_dyad, data = links10, xlab = 'dyad sex', las = 1, horizontal = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times')
boxplot(weight ~ dem_dyad, data = assoc10, xlab = 'dyad sex', las = 1, horizontal = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times and \nobserved together at least once')




### Calculate averages per individual and compare demographic groups








################ Second Half of Script: process sightings into dyadic associations -- very slow (~ 48 hours to run) ################
### import data
# elephant encounters
eles <- read_delim(file = 'data_processed/motnp_eles_long_22.01.06.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')
colnames(eles)[14] <- 'herd_type'
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]
str(eles)

# nodes
nodes <- read_delim(file = 'data_processed/motnp_elenodes_22.01.06.csv', delim = ',')
colnames(nodes)
str(nodes)

### create group-by-individual matrix
eles_asnipe <- eles[,c(3,4,2,5)]
eles_asnipe$Date <- as.integer(eles_asnipe$date)
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
which(is.na(eles_asnipe$Time))                                         # 161 698 1122 1469 1770
eles_asnipe[c(161,698,1122,1469,1770),]                                # all 1 sighting of B7+M44

eles_asnipe <- eles_asnipe[,c(5,6,3,4)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
str(eles_asnipe)

# get_gbi generates a group by individual matrix. The function accepts a data.table with individual identifiers and a group column. The group by individual matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_')
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter))
max(eles_asnipe$group) # 574 -- agrees with number of different sightings for which elephants were identified
eles_asnipe <- eles_asnipe[,c(3,7)]
eles_asnipe <- data.table::setDT(eles_asnipe)
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')

### code to convert gbi matrix format to dyadic data frame, shared by Prof Dan Franks and available from @JHart96 GitHub repository (https://github.com/JHart96/bison_examples/blob/main/examples/convert_gbi.md) -- NOTE: this step takes ~ 1.5 days to run
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1:nrow(gbi_matrix)) {
  for (i in which(gbi[obs_id, ] == 1)) {
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
write_delim(gbi_df, 'data_processed/motnp_bayesian_allpairwiseevents_22.01.06.csv', delim = ',')
# test <- read_delim('data_processed/motnp_bayesian_allpairwiseevents_22.01.06.csv', delim = ',')  # check that this has worked because you don't want to have to run it again!

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
readr::write_delim(dyads, 'data_processed/motnp_bayesian_trimmedpairwiseevents_22.01.10.csv', delim = ',')

