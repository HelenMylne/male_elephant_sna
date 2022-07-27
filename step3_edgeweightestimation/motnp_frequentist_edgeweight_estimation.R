#### Information ####
# Script to process association data from Mosi-Oa-Tunya National Park, Zambia.
# Combines together two scripts (21.12.07_ALERT_frequentist.R and 22.01.04_ALERT_bayesian.R) to go from raw data to analysable in a single script.
# 21.12.07_ALERT_frequentist.R : Remove anything based on frequentist SRI calculated values. Did have concerns that only 574 of the potential 1313 sightings were contained in the elephant data set ("eles_long") but I have confirmed that this is correct -- only 574 sightings had confirmed IDs.
# 22.01.04_ALERT_bayesian.R : Remove workings to leave only clean script and use full length eles_long input data instead of half-length.

# Data collected: 19th May 2016-21st December 2017
# Collected by: Mr David Youldon, Mr Dabwiso Sakala, Miss Helen Mylne and other volunteers/interns/facilitated research students working with ALERT during this time
# Data supplied by: ALERT and Mr David Youldon (11th August 2021) and via Volunteer Encounter (Bex Saunders, 19th October 2021)

# NOTE: Complete script will take up about 2 days to run on a standard 6-core i7 processor. Avoid running sections unnecessarily by checking for pre-saved data frames created in previous runs.

#### Set up ####
library(tidyverse)  # data manipulation
library(lubridate)  # sort date columns out
library(zoo)        # sort date columns out
library(asnipe)     # generating networks
library(rethinking) # col.alpha
library(igraph)     # plotting networks
#### Observation Sessions ####
sessions <- readxl::read_excel("../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 1)
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
write_delim(days, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_recording_days.csv',  delim = ';', col_names = T)
write_delim(s, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_recording_sessions.csv', delim = ';', col_names = T)

### clear environment
rm(days,s,sessions)

#### IDs ####
### import data
id.raw <- readxl::read_excel("../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 2)
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
ages <- readxl::read_excel("../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ALERT_ElephantAges_Saunders211019.xlsx", sheet = 1)  # read in age data
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
write_delim(id, "../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_id.csv", delim = ",", col_names = T)

### Clear Environment
rm(ages, days.numbers, id, id.raw, names, nkash, days, i, months, weeks)
#### Encounters ####
### import data
encounters <- readxl::read_excel("../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 3)
e <- encounters[c(1:3,8:27)]

### rename columns
e <- janitor::clean_names(e)
colnames(e)[c(4,5,7,8,23)] <- c('gps_s','gps_e','total_id','perc_id','present')
group.size <- paste('e',1:76, sep = '')

### separate column of elephants present into a whole set of columns
d <- e %>% separate(present,
                    as.character(c(1:max(e$total_id, na.rm = T))))
colnames(d)[23:98] <- group.size

# remove herd sightings (elephant IDs still present but removes herd ID -- e.g. a sighting of the matriarch from breeding herd 7 would show as "B7(F52)", and this will change to just say "F52")
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

### convert "total_elephants" to numeric -- can't just do "d$total_elephants <- as.numeric(d$total_elephants)" because makes all "50+" (or similar) values NA.
counts <- data.frame(total = d$total_elephants)
unique(counts$total) # convert manually everything with a + in it: "50+" "22+" "10+" "4+"  "27+" "30+" "55+" "15+" "3+"  "78+" "2+"  "61+" "37+" "52+" "40+" "60+" "5+"  "9+" "43+" "20+" "7+"  "11+"  "6+"  "23+" "12+" "42+" "18+" "25+" "70+" "28+" "14+" "17+" "77+" "97+" "13+" "16+" "19+" "51+" "24+" "1+"  "8+"  "36+" "34+" "41+" "59+" "88+" "29+" "66+" "32+" "21+" "38+"  "48+"  "39+"
plus <- data.frame(total = c("50+","22+","10+","4+","27+","30+","55+","15+","3+","78+","2+","61+","37+","52+","40+","60+","5+","9+","43+","20+","7+","11+","6+","23+","12+","42+","18+","25+","70+","28+","14+","17+","77+","97+","13+","16+","19+","51+","24+","1+","8+","36+","34+","41+","59+","88+","29+","66+","32+","21+","38+","48+","39+"))
plus$type <- 'minimum'
plus$value <- c(50,22,10,4,27,30,55,15,3,78,2,61,37,52,40,60,5,9,43,20,7,11,6,23,12,42,18,25,70,28,14,17,77,97,13,16,19,51,24,1,8,36,34,41,59,88,29,66,32,21,38,48,39)
counts <- left_join(x = counts, y = plus, by = 'total')
counts$type  <- ifelse(is.na(counts$value) == TRUE, 'estimate', 'minimum')
counts$value <- ifelse(is.na(counts$value) == TRUE, counts$total, counts$value)
d$total_elephants_numeric <- as.numeric(counts$value)
d$total_elephants_uncert  <- as.character(counts$type)

# recount total_id and perc_id -- some of these counts don't match to the number of elephants identified, and others the count given does not make up the correct percentage of those observed
d$count <- rowSums(!is.na(d[23:98]))
d$count_mismatch <- d$count-d$total_id
summary(d$count_mismatch) # range -7 to +8
# d$total_elephants <- as.numeric(d$total_elephants) # can't do this -- makes all where there is a "+" NA

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
write_delim(d, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_encounters_correctid.csv', na = 'NA', col_names = T, delim = ',')

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
d <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_encounters_correctid.csv')
for(j in 2:ncol(d)){
  for(i in 1:nrow(d)){
    d[i,j] <- ifelse(d[i,j] == 'NA',NA,d[i,j])
  }
}
str(d) # d$total_elephants left as character due to presence of some estimate values (e.g. 50+). If need to convert this to numeric, will first need to remove additional symbols.

### consider only the first sighting of each group, not the repeated measures every 15 minutes
typ1 <- d[d$typ == 1,] ; typNA <- d[is.na(d$typ),]
first <- rbind(typNA, typ1) ; first <- first[!is.na(first$date),]
first <- first[,c(1,3:12,26:101)]   # remove columns counting number of each demographic group
first$encounter <- 1:nrow(first) ; first <- first[,c(88,1:87)]   # index variable of sightings
first$date <- lubridate::as_date(first$date)

# make long, so every row is a sighting of an individual elephant
eles_long <- gather(data = first, key = number, value = elephant, e1:e76, factor_key = TRUE)
eles_long <- eles_long[!is.na(eles_long$elephant),]

# 2 incorrectly labelled elephants: no such code as N or D so N40 and D128 must be wrong. On keyboard, N next to M and D next to F --> assume that D128=F128 and N40=M40.
eles_long$elephant[which(eles_long$elephant == 'D128')] <- 'F128'  # row 2297
eles_long$elephant[which(eles_long$elephant == 'N40')] <- 'M40'    # row 1859

# write to csv for use in plotting
write_delim(eles_long, '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_eles_long_correctid.csv',col_names = T, delim = ',')

### make data frame that matches required structure for asnipe::get_associations_points_tw()
eles_asnipe <- eles_long[,c(1,3:5,14)]
eles_asnipe$location <- as.factor(paste(eles_asnipe$gps_s, eles_asnipe$gps_e, sep = '_'))
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

write_delim(sna_no_duplicates, "../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_dyads_noduplicates.csv",
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
ele_nodes <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_id.csv', delim = ',')
ele_nodes <- ele_nodes[,c(1:7,11:12)]

### edges data frame
ele_links <- read.csv("../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_dyads_noduplicates.csv", header=T, as.is=T)
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
eles_long <- read_delim('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_eles_long_correctid.csv', delim = ',')
eles_long$time <- as.numeric(eles_long$time)
eles_long$id_no <- eles_long$elephant
eles_long <- separate(eles_long, id_no, into = c('sex','num'), sep = 1)
eles_long$num <- sprintf("%04d", as.integer(eles_long$num))
eles_long$id_no <- paste(eles_long$sex, eles_long$num, sep='')
eles_long <- eles_long[,c(1:12,14:15,17)]

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
write_delim(ele_links, path = '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elelinks.csv', delim = ',', col_names = T)
write_delim(ele_nodes, path = '../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv', delim = ',', col_names = T)

### clear environment
rm(ele_links, ele_nodes, eles_long, i, j)
dev.off()

#### Plotting Networks ####
### load data
ele_links <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elelinks.csv')
ele_links$age.dyad   <- as.factor(ele_links$age.dyad)
ele_links$sex1       <- as.factor(ele_links$sex1)
ele_links$sex2       <- as.factor(ele_links$sex2)
ele_links$dem_class1 <- as.factor(ele_links$dem_class1)
ele_links$dem_class2 <- as.factor(ele_links$dem_class2)
str(ele_links)

ele_nodes <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_elenodes.csv')
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
t <- table(ele_links$dem_dyad)
barplot(t[order(t, decreasing=T)], horiz = T, las = 1, cex.names = 0.5,
        col = c(rep('blue',6), rep('grey',3), rep('blue',3),'grey',rep('blue',3),rep('grey',4),'blue',
                'grey','blue','grey','blue',rep('grey',7),'blue',rep('grey',8),'blue',rep('grey',9),'blue',
                rep('grey',13))) ; abline(v = 0)

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

kruskal.test(x = links10$weight, g = links10$type) # chi-squared = 38.864, df = 5, p-value = 2.529e-07
median(links10[links10$type == 'FF',]$weight) # 0.05416667
median(links10[links10$type == 'MF',]$weight) # 0.05582137
median(links10[links10$type == 'MM',]$weight) # 0.0659859
median(links10[links10$type == 'UF',]$weight) # 0.06660666
median(links10[links10$type == 'UM',]$weight) # 0.06950355
median(links10[links10$type == 'UU',]$weight) # 0.06666667
links10$type2 <- factor(links10$type, levels = c('FF','MF','MM','UF','UU','UM')) # this seems like the opposite order of what I'd expect...

### Dunn test
library(FSA)
DT <- dunnTest(weight ~ type2,
              data=links10,
              method="bh")      # Adjusts p-values for multiple comparisons;
DT
#Dunn (1964) Kruskal-Wallis multiple comparison

#Comparison           Z      P.unadj        P.adj
#1     FF - MF  0.53239138 5.944550e-01 8.106204e-01
##2     FF - MM -3.08662524 2.024426e-03 7.591599e-03 -- significantly different... from what...?
##3     MF - MM -4.93838488 7.877226e-07 1.181584e-05 -- significantly different... from what...?
##4     FF - UF -2.90885788 3.627517e-03 1.088255e-02 -- significantly different... from what...?
##5     MF - UF -3.93226180 8.415035e-05 6.311276e-04 -- significantly different... from what...?
##6     MM - UF -0.55020177 5.821810e-01 8.732715e-01 -- significantly different... from what...?
##7     FF - UM -2.55689173 1.056121e-02 2.640302e-02 -- significantly different... from what...?
##8     MF - UM -3.75522581 1.731853e-04 8.659263e-04 -- significantly different... from what...?
#9     MM - UM  0.17460071 8.613934e-01 9.229215e-01
#10    UF - UM  0.62733967 5.304366e-01 8.840610e-01
#11    FF - UU -1.54921314 1.213305e-01 2.274947e-01
#12    MF - UU -1.83578347 6.638969e-02 1.422636e-01
#13    MM - UU -0.28275307 7.773661e-01 8.969609e-01
#14    UF - UU -0.02132596 9.829856e-01 9.829856e-01
#15    UM - UU -0.34465882 7.303509e-01 9.129386e-01

### Pairwise Mann–Whitney
PT <- pairwise.wilcox.test(links10$weight, links10$type2, p.adjust.method="fdr")
PT
#    FF      MF      MM      UF      UM     
#MF 0.83372 -       -       -       -      
#MM 0.00838 7e-06   -       -       -      
#UF 0.01630 0.00071 0.73816 -       -      
#UM 0.03643 0.00071 0.98263 0.73816 -      
#UU 0.27019 0.15126 0.83372 0.98263 0.83372

PT <- pairwise.wilcox.test(links10$weight, links10$type2, p.adjust.method="holm")
PT
#    FF      MF      MM      UF      UM     
#MF 1.0000 -      -      -      -     
#MM 0.0268 7e-06  -      -      -     
#UF 0.0598 0.0018 1.0000 -      -     
#UM 0.1457 0.0018 1.0000 1.0000 -     
#UU 1.0000 0.6353 1.0000 1.0000 1.0000

### Compare age dyads
assoc10$age.dyad <- factor(assoc10$age.dyad, levels = c('AA','AP','AJ','AC','PP','PJ','PC','JJ','CJ','CC'))
boxplot(weight ~ age.dyad, data = ele_links, xlab = 'dyad age', las = 1, notch = T,
        main = 'association strength including all dyads')
boxplot(weight ~ age.dyad, data = associations, xlab = 'dyad age', las = 1, notch = T,
        main = 'association strength including all dyads\nobserved together at least once')
boxplot(weight ~ age.dyad, data = links10, xlab = 'dyad sex', las = 1, notch = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times')
boxplot(weight ~ age.dyad, data = assoc10, xlab = 'dyad sex', las = 1, notch = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times and \nobserved together at least once')

kruskal.test(x = links10$weight, g = links10$age.dyad) # chi-squared = 91.329, df = 9, p-value = 8.807e-16
median(links10[links10$age.dyad == 'CC',]$weight) # 0.03697312
median(links10[links10$age.dyad == 'PP',]$weight) # 0.0407855
median(links10[links10$age.dyad == 'PC',]$weight) # 0.05359877
median(links10[links10$age.dyad == 'AP',]$weight) # 0.05674847
median(links10[links10$age.dyad == 'JJ',]$weight) # 0.05714286
median(links10[links10$age.dyad == 'PJ',]$weight) # 0.06395825
median(links10[links10$age.dyad == 'AC',]$weight) # 0.06438896
median(links10[links10$age.dyad == 'AA',]$weight) # 0.06738869
median(links10[links10$age.dyad == 'AJ',]$weight) # 0.06833717
median(links10[links10$age.dyad == 'CJ',]$weight) # 0.07648137
links10$age_dyad <- factor(links10$age.dyad, levels = c('CC','PP','PC','AP','JJ','PJ','AC','AA','AJ','CJ'))

males <- links10[links10$type == 'MM',]
males <- males[males$age.dyad == 'AA' | males$age.dyad == 'AP' | males$age.dyad == 'PP',]
males <- males[males$age1 != '9-10' & males$age2 != '9-10',]
males$age1.fct <- as.integer(as.factor(males$age1))
males$age2.fct <- as.integer(as.factor(males$age2))
males$age_pair <- as.factor(ifelse(males$age1.fct < males$age2.fct,
                                   paste(males$age1.fct, males$age2.fct, sep = '_'),
                                   paste(males$age2.fct, males$age1.fct, sep = '_')))
males
table(males$age_pair)
par(cex.axis = 1)
boxplot(males$weight ~ males$age_pair, las = 1, xlab = 'weight', ylab = 'age pairing (1 = 10-15y, 2 = 16-19y, 3 = 20-25y, 4 = 26-39y, 5 = 40+', horizontal = T)

males$age_pair_words <- ifelse(males$age_pair == '1_1','yp-yp',
                               ifelse(males$age_pair == '1_2','yp-op',
                                      ifelse(males$age_pair == '1_3','yp-ya',
                                             ifelse(males$age_pair == '1_4','yp-ma',
                                                    ifelse(males$age_pair == '1_5','yp-oa','other')))))
males$age_pair_words <- ifelse(males$age_pair_words != 'other', males$age_pair_words, 
                               ifelse(males$age_pair == '2_2','op-op',
                                      ifelse(males$age_pair == '2_3','op-ya',
                                             ifelse(males$age_pair == '2_4','op-ma',
                                                    ifelse(males$age_pair == '2_5','op-oa', 'other')))))
males$age_pair_words <- ifelse(males$age_pair_words != 'other', males$age_pair_words, 
                               ifelse(males$age_pair == '3_3','ya-ya',
                                      ifelse(males$age_pair == '3_4','ya-ma',
                                             ifelse(males$age_pair == '3_5','ya-oa', 'other'))))
males$age_pair_words <- ifelse(males$age_pair_words != 'other', males$age_pair_words,
                               ifelse(males$age_pair == '4_4','ma-ma',
                                      ifelse(males$age_pair == '4_5','ma-oa', 'other')))
males$age_pair_words <- ifelse(males$age_pair_words == 'other', ifelse(males$age_pair == '5_5','oa-oa','check'), males$age_pair_words)
unique(males$age_pair_words)
which(males$age_pair_words == 'other')
which(males$age_pair_words == 'check')

par(mar = c(5,6,4,2), mgp = c(3,1,0))
males$age_pair_words <- factor(males$age_pair_words,
                               levels = c('oa-oa',
                                          'ma-oa','ma-ma',
                                          'ya-oa','ya-ma','ya-ya',
                                          'op-oa','op-ma','op-ya','op-op',
                                          'yp-oa','yp-ma','yp-ya','yp-op','yp-yp'))
boxplot(males$weight ~ males$age_pair_words, ylim = c(0,0.5), #notch = T, 
        ylab = 'age pairing (yp = young pubescent, ma = mid-adult, oa = old adult) \n', horizontal = T,
        las = 1, xlab = 'weight', outpch = 16, outcol = 'grey',
        col = c('skyblue','grey','grey','grey','grey','skyblue','grey','grey','grey',
                'skyblue','grey','grey','skyblue','grey','skyblue'))

ggplot(males, aes(x = weight, y = age_pair_words))+
  geom_boxplot(fill = c('skyblue','grey','grey','grey','grey','skyblue','grey','grey','grey',
                        'skyblue','grey','grey','skyblue','grey','skyblue'))+
  theme_classic()+
  labs(y = 'age pairing (yp = young pubescent, ma = mid-adult, oa = old adult)')+
  scale_y_discrete(limits = rev(levels(males$age_pair_words)))+
  scale_x_continuous(limits = c(-0.03,0.5))+
  geom_curve(aes(x = -0.015, y = 'oa-oa', xend = -0.015, yend = 'yp-yp'),
             arrow = arrow(length = unit(2,'mm'), type = 'closed', ends = 'first'),
             curvature = 0)+
  geom_text(aes(x = -0.03, y = 7.5, label = "increasing age"), stat = "unique", angle = 90)

### Pairwise Mann–Whitney
PT <- pairwise.wilcox.test(links10$weight, links10$age_dyad, p.adjust.method="fdr")
PT
#   CC      PP      PC      AP      JJ      PJ      AC      AA      AJ     
#PP 0.55580 -       -       -       -       -       -       -       -      
#PC 0.00935 0.12663 -       -       -       -       -       -       -      
#AP 0.00015 0.02819 0.74672 -       -       -       -       -       -      
#JJ 0.06972 0.11050 0.35274 0.38307 -       -       -       -       -      
#PJ 0.00034 0.00737 0.16123 0.16726 0.74431 -       -       -       -      
#AC 9.1e-07 0.00145 0.12663 0.06513 0.64603 0.69738 -       -       -      
#AA 3.1e-11 1.7e-05 0.00566 1.9e-05 0.77569 0.83265 0.12213 -       -      
#AJ 9.5e-09 3.5e-05 0.00835 0.00145 0.84658 0.70764 0.12663 0.66911 -      
#CJ 1.9e-05 0.00036 0.01175 0.00865 0.78466 0.24237 0.07173 0.16726 0.32600

PT <- pairwise.wilcox.test(links10$weight, links10$age_dyad, p.adjust.method="holm")
PT
#   CC      PP      PC      AP      JJ      PJ      AC      AA      AJ     
#PP 1.00000 -       -       -       -       -       -       -       -      
#PC 0.10240 1.00000 -       -       -       -       -       -       -      
#AP 0.00099 0.32133 1.00000 -       -       -       -       -       -          AP ≠ CC
#JJ 0.81344 1.00000 1.00000 1.00000 -       -       -       -       -      
#PJ 0.00254 0.07339 1.00000 1.00000 1.00000 -       -       -       -          PJ ≠ CC
#AC 2.6e-06 0.01340 1.00000 0.75264 1.00000 1.00000 -       -       -          AC ≠ CC or PP
#AA 3.1e-11 6.3e-05 0.05397 0.00010 1.00000 1.00000 1.00000 -       -          AA ≠ CC or PP or AP
#AJ 1.9e-08 0.00022 0.08628 0.01340 1.00000 1.00000 1.00000 1.00000 -          AJ ≠ CC or PP or AP
#CJ 0.00010 0.00285 0.13160 0.09226 1.00000 1.00000 0.84163 1.00000 1.00000    CJ ≠ CC or PP


### Compare full demography
links10$dem_dyad <- as.factor(links10$dem_dyad)
boxplot(weight ~ dem_dyad, data = ele_links, ylab = 'dyad type', las = 1, horizontal = T,
        main = 'association strength including all dyads')
boxplot(weight ~ dem_dyad, data = associations, ylab = 'dyad age', las = 1, horizontal = T,
        main = 'association strength including all dyads\nobserved together at least once')
boxplot(weight ~ dem_dyad, data = links10, ylab = 'dyad sex', las = 1, horizontal = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times')
boxplot(weight ~ dem_dyad, data = assoc10, ylab = 'dyad sex', las = 1, horizontal = T,
        main = 'association strength including all dyads where\nboth members were observed at least 10 times and \nobserved together at least once')

kruskal.test(x = links10$weight, g = links10$dem_dyad) # chi-squared = 337.54, df = 52, p-value < 2.2e-16
tapply(X = links10$weight, INDEX = links10$dem_dyad, FUN = median)
links10$dem_dyad2 <- factor(links10$dem_dyad,
                            levels = c('AF_CF','AF_JF','AM_CF','CF_CM','CF_CU','CM_CM','JF_JM','JM_CF','JU_CF','PF_CF','PM_CF',
                                       'JF_CU',
                                       'JF_CM',
                                       'PF_CM','PM_PM',
                                       'AM_JF','AF_CM',
                                       'PF_JF','PF_PF','PF_PM','AM_PF',
                                       'AM_CM','PM_JF','PM_CM','AF_PF','CM_CU','AF_PM','PF_JU',
                                       'AF_JU','JM_CM','JU_CU','JM_JU','AM_PM','AF_AF','AF_AM','PF_CU','PM_JM','PM_JU','CU_CU','JU_CM','PM_CU','AM_JU','AM_CU',
                                       'AM_JM','AF_CU','PF_JM','AM_AM','JU_JU',
                                       'AF_JM',
                                       'JF_JU','JM_CU','JF_CF','JM_JM'))
boxplot(weight ~ dem_dyad2, data = links10, ylab = 'dyad sex', las = 1, horizontal = T, par(cex.axis = 0.5),
        main = 'association strength including all dyads where\nboth members were observed at least 10 times')


### Pairwise Mann–Whitney
PT <- pairwise.wilcox.test(links10$weight, links10$dem_dyad2, p.adjust.method="holm", paired = F)
PT
## sig diffs:
# AF_CF: PF_PM, AM_PF, AM_CM, AF_PF, CM_CU, AF_PM, PF_JU, AF_JU, JU_CU, AM_PM, AF_AF, AF_AM, PF_CU, PM_JM, PM_JU, CU_CU, PM_CU, AM_JU, AM_CU, AM_JM, AF_CU, PF_JM, AM_AM, AF_JM, JM_CU
# AF_JF:
# AM_CF: PF_CM, PM_PM, AF_CM, PF_PF, PF_PM, AM_PF, AM_CM, JM_CU, PM_CM, AF_PF, CM_CU, AF_PM, PF_JU, AF_JU, JM_CM, JU_CU, JM_JU, AM_PM, AF_AF, AF_AM, PF_CU, PM_JM, PM_JU, CU_CU, JU_CM, PM_CU, AM_JU, AM_CU, AM_JM, AF_CU, PF_JM, AM_AM
# CF_CM: JM_CU
# CF_CU: JM_CU
# CM_CM: JM_CU
# JF_JM:
# JM_CF:
# JU_CF:
# PF_CF: AM_PM, PM_JM, PM_JU, PM_CU, AM_JU, AM_CU, AM_JM, AF_CU, AM_AM, AF_JM, JM_CU
# PM_CF: AM_JM, AM_AM, AF_JM, JM_CU
# JF_CU:
# JF_CM:
# PF_CM: AM_AM, JM_CU
# PM_PM: AF_JM, JM_CU
# AM_JF:
# AF_CM: AF_CU, AM_AM, AF_JM, JM_CU
# PF_JF: 
# PF_PF: AF_JM, JM_CU
# PF_PM: AF_JM, JM_CU
# AM_PF: AF_CU, AM_AM, AF_JM, JM_CU
# AM_CM: AF_CU, AM_AM, AF_JM, JM_CU
# PM_JF:
# PM_CM: AF_JM, JM_CU
# AF_PF: AM_AM, AF_JM, JM_CU
# CM_CU: 
# AF_PM: AF_JM, JM_CU
# PF_JU:
# AF_JU:
# JM_CM:
# JU_CU:
# JM_JU:
# AM_PM: AM_AM, AF_JM, JM_CU
# AF_AF:
# AF_AM: AM_AM, AF_JM, JM_CU
# PF_CU:
# PM_JM: 
# PM_JU:
# CU_CU:
# JU_CM:
# PM_CU: 
# AM_JU: JM_CU
# AM_CU: JM_CU
# AM_JM:
# AF_CU:
# PF_JM:
# AM_AM:
# JU_JU: 
# AF_JM:
# JF_JU: 
# JM_CU: 
# JF_CF: 





### Calculate averages per individual and compare demographic groups







