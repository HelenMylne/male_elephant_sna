## Bayesian analysis of ALERT data
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

################ First Half of Script: recreating and checking data produced in 21.12.07_ALERT_frequentist.R ################
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

### Add age data into id sheet
id <- merge(x = id, y = ages)             # 468 elephants
id <- id[,c(1:5,369,6:10,370,11:368)]

### Remove elephants with wrong name in ID/ages sheet -- confusion over names suggests potentially incorrect sightings data
id <- id[id$id_no != 'F0157',] ; ages <- ages[ages$id_no != 'F0157',]
id <- id[is.na(id$comments),]

### Write out processed data
write_delim(id, "data_processed/motnp_id_22.01.06.csv", delim = ",", col_names = T)

### Clear Environment
rm(ages, days.numbers, id, id.raw, names, nkash, days, i, months, weeks)
#### Encounters ####
### import data
encounters <- readxl::read_excel("data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 3)
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

### recount total_id and perc_id -- some of these counts don't match to the number of elephants identified, and others the count given does not make up the correct percentage of those observed
d$count <- rowSums(!is.na(d[23:98]))
d$count_mismatch <- d$count-d$total_id
summary(d$count_mismatch) # range -7 to +8
d$perc_id_new <- 100*d$count / d$total_elephants_numeric
d <- d[c(1:6,99:100,7,101:102,8,103,9:98)]
colnames(d)[9:13] <- c('total_id_dy','total_id_hkm','total_id_diff','perc_id_dy','perc_id_hkm')
str(d)

### correct mis-recorded locations -- 3 which are an order of magnitude from the rest in either east or south direction, one for which the east value started 17 instead of 25 but otherwise did not match south so assume that changing the 25 to 17 was error in recording.
d$gps_e[232] <- d$gps_e[232]/10 # value was 25470333 -- 2547033 puts it within cluster of all other points (assume typo)
d$gps_e[480] <- d$gps_e[480]*10 # value was 254665 -- 2546650 puts it within cluster of all other points (assume Excel treated as shorter number)
d$gps_s[672] <- d$gps_s[672]*10 # value was 175296 -- 1752960 puts it within cluster of all other points (assume Excel treated as shorter number)
d$gps_e[624] <- d$gps_e[624]+800000 # value started 17 rather than 25 (1745269), but then matched format of other east coordinates -- 2545269 puts it within cluster of all other points

### write csv
write_delim(d, 'data_processed/motnp_encounters_correctid_22.01.06.csv', na = 'NA', col_names = T, delim = ',')

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

#### Convert sightings data to long format for analysis -- checked, no error in creating eles_long, 574 is correct ####
### import encounter data
d <- read_csv('data_processed/motnp_6thJan2022/motnp_encounters_correctid_22.01.06.csv')
for(j in 2:ncol(d)){
  for(i in 1:nrow(d)){
    d[i,j] <- ifelse(d[i,j] == 'NA',NA,d[i,j])
  }
}
str(d) # d$total_elephants left as character due to presence of some estimate values (e.g. 50+). If need to convert this to numeric, will first need to remove additional symbols.

### consider only the first sighting of each group, not the repeated measures every 15 minutes
typ1 <- d[d$typ == 1,]                # anything measured as typ ≥2 is a repeat sighting of the same individuals
typNA <- d[is.na(d$typ),]             # Did not start to repeat sightings data every 15 minutes immediately, so typ=NA before we started this
first <- rbind(typNA, typ1)           # combine those that are first sightings into a single data frame
first <- first[!is.na(first$date),]   # remove data points where no date recorded
first <- first[,c(1,3:14,28:103)]     # remove columns counting number of each demographic group -- unnecessary and inaccurate
first$encounter <- 1:nrow(first) ; first <- first[,c(90,1:89)]   # index variable of sightings
first$date <- lubridate::as_date(first$date)                     # put date column into correct format
length(unique(first$encounter))       # 1078 separate encounters

### make long, so every row is a sighting of an individual elephant
eles_long <- gather(data = first, key = number, value = elephant, e1:e76, factor_key = TRUE)
length(unique(eles_long$encounter))   # 1078
eles_long <- eles_long[!is.na(eles_long$elephant),] # filter out NA values: previously created 76 columns of individuals per sighting, but most columns were blank (only one sighting had 76 identified elephants) so now have to remove all blank cells
length(unique(eles_long$encounter))   # 574 -- this is where half the elephants disappear. 574 encounters is correct despite there being 1313 total sightings recorded in raw data of which 1078 independent, because only 574 of these 1078 independent sightings contained identifiable individuals.
length(first[first$total_id_hkm > 0,]$encounter) # 574 observations -- confirm that 574 is the correct number of sightings for which an ID was possible.

### 2 incorrectly labelled elephants: no such code as N or D so N40 and D128 must be wrong. On keyboard, N next to M and D next to F --> assume that D128=F128 and N40=M40 (M40 one of the most commonly observed elephants).
eles_long$elephant[which(eles_long$elephant == 'D128')] <- 'F128'  # row 2297
eles_long$elephant[which(eles_long$elephant == 'N40')] <- 'M40'    # row 1859

### write to csv for use in plotting
write_delim(eles_long, 'data_processed/motnp_eles_long_22.01.06.csv',col_names = T, delim = ',')

### clear environment
rm(d, typ1, typNA, first, eles_long)
#### Create frequentist SRI dyads using group-by-individual matrix (asnipe) -- quite slow to run (~ 10 minutes) #####
### import data
eles_long <- read_delim('data_processed/motnp_eles_long_22.01.06.csv', delim = ',')
str(eles_long)

### make data frame that matches required structure for asnipe::get_associations_points_tw()
eles_long <- eles_long[,c(2:5,16)]
eles_long$location <- as.factor(paste(eles_long$gps_s, eles_long$gps_e, sep = '_'))
eles_long <- eles_long[,c(1,2,5,6)]
colnames(eles_long) <- c('Date','Time','ID','Location')
eles_long$Date <- as.numeric(eles_long$Date)
eles_long$ID <- as.factor(eles_long$ID)

eles_long$Time <- ifelse(eles_long$Time > 1, NA, eles_long$Time) # time = proportion of day so anything >1 has to be wrong. 1 visit is recorded at time 1599 -- cannot be corrected, so convert to missing value.
eles_long$Time <- eles_long$Time*86400     # convert time values to seconds so actually readable
hist(eles_long$Time, breaks = 30, las = 1) # can see the times we were out from the lunchtime reduction

# generate group-by-individual matrix
a <- get_associations_points_tw(eles_long,            # data to input
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
colnames(sna_no_duplicates)
sna_no_duplicates <- sna_no_duplicates[,c(1,2,5,3,6)]

sna_no_duplicates$sri_frequentist <- NA
for(i in 1:nrow(sna_no_duplicates)) {
  sna_no_duplicates$sri_frequentist[i] <- sna$sri[match(x = sna_no_duplicates$dyad[i], table = sna$dyad)]
}

### plot histogram of average cost
associations <- sna_no_duplicates[sna_no_duplicates$sri_frequentist > 0,] # considering only dyads for which there was at least some association
associations <- associations[!is.na(associations$sri_frequentist),]

par(mfrow = c(2,1))
hist(sna_no_duplicates$sri_frequentist, xlab = 'SRI', las = 1, ylab = '', main = 'Including all dyads')
hist(associations$sri_frequentist, xlab = 'SRI', las = 1, ylab = '', main = 'After removing unassociated dyads')

### write to csv file
write_delim(sna_no_duplicates, "data_processed/motnp_dyads_22.01.06.csv",
            delim = ",", na = "", append = FALSE, col_names = T, quote_escape = "double")

### clear environment, reset plot window
rm(a,associations,d,eles_asnipe,eles_long,first,gbi,n,sna,sna_no_duplicates,typ1,typNA,i,j)
par(mfrow = c(1,1))

#### Create frequentist plotting data -- slow to run (~ 45 minutes), avoid repeating unnecessarily ####
### notes data frame
ele_nodes <- read_delim('data_processed/motnp_id_22.01.06.csv', delim = ',')
ele_nodes <- ele_nodes[,c(1:7,11:12)]

### edges data frame
ele_links <- read.csv("data_processed/motnp_dyads_22.01.06.csv", header=T, as.is=T)
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

### make long (zero-padded) label for each elephant in eles_long to match ele_nodes
eles_long <- read_delim('data_processed/motnp_eles_long_correctid.csv', delim = ',')
eles_long$time <- as.numeric(eles_long$time)
eles_long$id_no <- eles_long$elephant
eles_long <- separate(eles_long, id_no, into = c('sex','num'), sep = 1)
eles_long$num <- sprintf("%04d", as.integer(eles_long$num))
eles_long$id_no <- paste(eles_long$sex, eles_long$num, sep='')
eles_long <- eles_long[,c(1:12,14:15,17)]

### make short (unpadded) label for each elephant in ele_nodes to match eles_long
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
ele_links <- ele_links[which(ele_links$from != 'M0021'),] # remove M0013 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'M0021'),] # remove M0013 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'M0125'),] # remove M0125 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'M0125'),] # remove M0125 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'M0138'),] # remove M0138 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'M0138'),] # remove M0138 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'M0223'),] # remove M0138 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'M0223'),] # remove M0138 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$from != 'M0227'),] # remove M0138 from links -- doesn't exist in IDs
ele_links <- ele_links[which(ele_links$to   != 'M0227'),] # remove M0138 from links -- doesn't exist in IDs
ele_nodes <- ele_nodes[which(ele_nodes$id_no!= 'M0175'),] # remove M0175 from nodes -- doesn't exist in sightings

length(unique(ele_nodes$id_no)) ; length(unique(ele_links$from)) ; length(unique(ele_links$to)) # 463, 462, 462 (One less for links than nodes because F1 is only in "from" and U9 is only in "to". No elephant is from and to itself.)

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
# WARNING -- THIS STEP TAKES ABOUT 30-40 MINUTES TO RUN
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
write_delim(ele_links, path = 'data_processed/motnp_elelinks_22.01.06.csv', delim = ',', col_names = T)
write_delim(ele_nodes, path = 'data_processed/motnp_elenodes_22.01.06.csv', delim = ',', col_names = T)

### clear environment
rm(ele_links, ele_nodes, eles_long, i, j)
