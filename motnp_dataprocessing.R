## Bayesian analysis of ALERT data
#### Information ####
# Script to process association data from Mosi-Oa-Tunya National Park, Zambia.

# Data collected: 19th May 2016-21st December 2017
# Collected by: Mr David Youldon, Mr Dabwiso Sakala, Miss Helen Mylne and other volunteers/interns/facilitated research students working with ALERT during this time
# Data supplied by: ALERT and Mr David Youldon (11th August 2021) and via Volunteer Encounter (Bex Saunders, 19th October 2021)

# NOTE: Complete script will take up about 2 days to run on a standard 6-core i7 processor.

#### Set up ####
library(tidyverse)  # data manipulation
library(lubridate)  # sort date columns out
library(zoo)        # sort date columns out
library(asnipe)     # generating networks

################ 1) Create basic data frames from raw data ################
#### Observation Sessions ####
sessions <- readxl::read_excel("data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 1) # read in raw data
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
write_delim(days, 'data_processed/motnp_recording_days_22.01.13.csv',  delim = ';', col_names = T)
write_delim(s, 'data_processed/motnp_recording_sessions_22.01.13.csv', delim = ';', col_names = T)

### clear environment
rm(days,s,sessions)

#### IDs ####
### import data
id.raw <- readxl::read_excel("data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 2) # read in raw data
id <- id.raw[5:849,1:368]             # create a new data frame without top values or empty end columns
id <- id[!is.na(id[,1]),]             # remove empty rows
id <- id[!is.na(id[,5]),]             # remove unassigned numbers

### rename columns
colnames(id)[1:10] <- id.raw[4,1:10]  # rename first columns with the names from the original

months <- seq(from = dmy('1st May 2016'), to = dmy('1st Dec 2017'), by = 'month')  # create vector of days, 1st day of each month
months <- paste('m',as.character(months), sep = '_')  # replace vector values with something readable
colnames(id)[11:30] <- months # label month columns with user-friendly dates of 1st day

weeks <- seq(from = dmy('19th May 2016'), to = dmy('15th Nov 2017'), by = 'weeks') # create vector of days, every Thursday from start to end
weeks <- paste('w',as.character(weeks), sep = '_') # replace weeks vector values with readable dates
colnames(id)[31:108] <- weeks       # label week columns with user-friendly dates, starting each week on Thursdays

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

### Add values for missing elephants -- ages estimated from additional images sent 11th November 2021
ages$id_no[which(ages$age_category == 'missing')]       # "M0026" "M0154" "M0196" "M0240"
ages$age_class[which(ages$age_category == 'missing')]   # calf, NA, adult, pubescent

ages$age_category[which(ages$id_no == 'M0026')] <- '0-3'
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
id <- id[,c(1:5,369,6:10,370,11:368)]     # rearrange columns

### Remove elephants with wrong name in ID/ages sheet -- confusion over names suggests potentially incorrect sightings data
id <- id[id$id_no != 'F0157',] ; ages <- ages[ages$id_no != 'F0157',]
id <- id[is.na(id$comments),]

### Write out processed data
write_delim(id, "data_processed/motnp_id_22.01.13.csv", delim = ",", col_names = T)

### Clear Environment
rm(ages, days.numbers, id, id.raw, names, nkash, days, i, months, weeks)
#### Encounters ####
### import data
encounters <- readxl::read_excel("data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 3) # read in raw data
e <- encounters[c(1:3,8:27)]    # only necessary columns

### rename columns
e <- janitor::clean_names(e)
colnames(e)[c(4,5,7,8,23)] <- c('gps_s','gps_e','total_id','perc_id','present')
group.size <- paste('e',1:76, sep = '')

### separate column of elephants present into a whole set of columns
d <- e %>% separate(present, as.character(c(1:max(e$total_id, na.rm = T))))
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

### convert "total_elephants" to numeric (first need to remove "+" from any estimated values)
counts <- data.frame(total = d$total_elephants)
unique(counts$total) # convert manually everything with a + in it: "50+" "22+" "10+" "4+"  "27+" "30+" "55+" "15+" "3+"  "78+" "2+"  "61+" "37+" "52+" "40+" "60+" "5+"  "9+" "43+" "20+" "7+"  "11+"  "6+"  "23+" "12+" "42+" "18+" "25+" "70+" "28+" "14+" "17+" "77+" "97+" "13+" "16+" "19+" "51+" "24+" "1+"  "8+"  "36+" "34+" "41+" "59+" "88+" "29+" "66+" "32+" "21+" "38+"  "48+"  "39+"
plus <- data.frame(total = c("50+","22+","10+","4+","27+","30+","55+","15+","3+","78+","2+","61+","37+","52+","40+","60+","5+","9+","43+","20+","7+","11+","6+","23+","12+","42+","18+","25+","70+","28+","14+","17+","77+","97+","13+","16+","19+","51+","24+","1+","8+","36+","34+","41+","59+","88+","29+","66+","32+","21+","38+","48+","39+"))
plus$type <- 'minimum' # anything with a plus is the minimum number that could have been there
plus$value <- c(50,22,10,4,27,30,55,15,3,78,2,61,37,52,40,60,5,9,43,20,7,11,6,23,12,42,18,25,70,28,14,17,77,97,13,16,19,51,24,1,8,36,34,41,59,88,29,66,32,21,38,48,39) # set values without +
counts <- left_join(x = counts, y = plus, by = 'total') # add values back into counts data frame
counts$type  <- ifelse(is.na(counts$value) == TRUE, 'estimate', 'minimum')      # if not a minimum, then actual estimate
counts$value <- ifelse(is.na(counts$value) == TRUE, counts$total, counts$value) # set values to total numbers where NA
d$total_elephants_numeric <- as.numeric(counts$value)  # make numeric
d$total_elephants_uncert  <- as.character(counts$type) # indicate whether count is confident or not

### recount total_id and perc_id -- some of these counts don't match to the number of elephants identified, and others the count given does not make up the correct percentage of those observed
d$count <- rowSums(!is.na(d[23:98]))     # recount number of elephants identified as present (skip over NA values)
d$count_mismatch <- d$count-d$total_id   # identify difference between ALERT reported count and total identified
summary(d$count_mismatch)                # range -7 to +8
d$perc_id_new <- 100*d$count / d$total_elephants_numeric # New percentage of elephants identified
d <- d[c(1:6,99:100,7,101:102,8,103,9:98)]               # rearrange columns
colnames(d)[9:13] <- c('total_id_dy','total_id_hkm','total_id_diff','perc_id_dy','perc_id_hkm') # rename columns: dy = David Youldon (original values), hkm = Helen Mylne (new values)
str(d)

### correct mis-recorded locations -- 3 which are an order of magnitude from the rest in either east or south direction, one for which the east value started 17 instead of 25 but otherwise did not match south so assume that changing the 25 to 17 was error in recording.
d$gps_e[232] <- d$gps_e[232]/10 # value was 25470333 -- 2547033 puts it within cluster of all other points (assume typo)
d$gps_e[480] <- d$gps_e[480]*10 # value was 254665 -- 2546650 puts it within cluster of all other points (assume Excel treated as shorter number)
d$gps_s[672] <- d$gps_s[672]*10 # value was 175296 -- 1752960 puts it within cluster of all other points (assume Excel treated as shorter number)
d$gps_e[624] <- d$gps_e[624]+800000 # value started 17 rather than 25 (1745269), but then matched format of other east coordinates -- 2545269 puts it within cluster of all other points

### correct mis-recorded herd types -- some MO that should be BH or MX
# In MO sightings: "F102" "F104" "F107" "F19"  "F21"  "F57"  "F63"  "F64"  "F68"  "F99"
d$unique <- 1:nrow(d)

mo <- d[d$type == 'MO',] # any females recorded in MO sighting should be changed to MX.
sort(unique(mo$e1))    # F102, F19, F63, F68, U66
sort(unique(mo$e2))    # F107, F57, F64
sort(unique(mo$e3))    # F21, F99
sort(unique(mo$e4))    # F104
sort(unique(mo$e5))    # No more females -- Females always recorded in spreadsheet first so will come through in first columns. If none by e5, will be none with higher e values.
which(mo$e1 == 'F102') # 238
which(mo$e2 == 'F57')  # 238
which(mo$e1 == 'F19')  # 286
which(mo$e3 == 'F21')  # 286
which(mo$e4 == 'F104') # 286
which(mo$e1 == 'F63')  # 26
which(mo$e2 == 'F64')  # 26
which(mo$e3 == 'F99')  # 26
which(mo$e1 == 'F68')  # 41
which(mo$e1 == 'U66')  # 403
which(mo$e2 == 'F107') # 118, 119

mo$unique[which(mo$e1 == 'F63')]  # 50
mo$unique[which(mo$e2 == 'F64')]  # 50
mo$unique[which(mo$e3 == 'F99')]  # 50
mo$unique[which(mo$e1 == 'F68')]  # 88
mo$unique[which(mo$e2 == 'F107')] # 276 277
mo$unique[which(mo$e1 == 'F102')] # 628
mo$unique[which(mo$e2 == 'F57')]  # 628
mo$unique[which(mo$e1 == 'F19')]  # 736
mo$unique[which(mo$e3 == 'F21')]  # 736
mo$unique[which(mo$e4 == 'F104')] # 736
mo$unique[which(mo$e1 == 'U66')]  # 1006

d$type[which(d$unique == 50)]   <- 'MX'  # 3 females + 5 males
d$type[which(d$unique == 88)]   <- 'MX'  # adult female + adult male
d$type[which(d$unique == 276)]  <- 'MX'  # 6 males + female
d$type[which(d$unique == 277)]  <- 'MX'  # 6 males + female
d$type[which(d$unique == 628)]  <- 'BH'  # adult female + calf
d$type[which(d$unique == 736)]  <- 'BH'  # breeding herd 3 + additional female
d$type[which(d$unique == 1006)] <- 'UK'  # a calf with a whole group of males -- this is too unclear to reassign

View(d[c(50,88,276,277,628,736,1006),]) # still an issue with 50 -- 1 elephant, but 8 identified, suggests that these may have got scrambled together somehow, but nothing to be done about it

bh <- d[d$type == 'BH',]
eles_long <- read_csv('data_processed/motnp_eles_long_22.01.06.csv')
ele_nodes <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
eles_bh <- eles_long[eles_long$type == 'BH',]
sort(unique(eles_bh$elephant)) # "M103" "M12"  "M120" "M126" "M137" "M139" "M148" "M160" "M164" "M174" "M176" "M179" "M180" "M182" "M187" "M189" "M19"  "M190" "M191" "M2"   "M203" "M204" "M206" "M215" "M217" "M225" "M23"  "M238" "M24"  "M240" "M245" "M248" "M249" "M26" "M6"   "M60"  "M69"  "M7"   "M71"  "M72"  "M84"  "M86"  "M87"  "M88"  "M9"   "M96"  "M97" 
# check which of these are adult males and therefore shouldn't be in BH:
ele_nodes[which(ele_nodes$id_no == 'M0103'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0012'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0120'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0126'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0137'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0139'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0148'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0160'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0164'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0174'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0176'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0179'),c(4,5)] # Juvenile -- fine
ele_nodes[which(ele_nodes$id_no == 'M0180'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0182'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0187'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0189'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0019'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0190'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0191'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0002'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0203'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0204'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0206'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0215'),c(4,5)] # Juvenile -- fine
ele_nodes[which(ele_nodes$id_no == 'M0217'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0225'),c(4,5)] # Juvenile -- fine
ele_nodes[which(ele_nodes$id_no == 'M0023'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0238'),c(4,5)] # Juvenile -- fine
ele_nodes[which(ele_nodes$id_no == 'M0024'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0240'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0245'),c(4,5)] # Juvenile -- fine
ele_nodes[which(ele_nodes$id_no == 'M0248'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0249'),c(4,5)] # Juvenile -- fine
ele_nodes[which(ele_nodes$id_no == 'M0026'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0006'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0069'),c(4,5)] # Juvenile -- fine
ele_nodes[which(ele_nodes$id_no == 'M0007'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0071'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0072'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0084'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0086'),c(4,5)] # Calf -- fine
ele_nodes[which(ele_nodes$id_no == 'M0087'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0088'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0009'),c(4,5)] # Adult -- change to MX
ele_nodes[which(ele_nodes$id_no == 'M0096'),c(4,5)] # Pubescent -- fine
ele_nodes[which(ele_nodes$id_no == 'M0097'),c(4,5)] # Pubescent -- fine
# adult males listed in sightings of breeding herds = M12, M120, M139, M174, M190, M191, M2, M203, M217, M88, M9

which(eles_bh$elephant == 'M12')  # 461
which(eles_bh$elephant == 'M120') # 501
which(eles_bh$elephant == 'M139') # 477
which(eles_bh$elephant == 'M174') # 484
which(eles_bh$elephant == 'M190') # 493
which(eles_bh$elephant == 'M191') # 496
which(eles_bh$elephant == 'M2')   # 433
which(eles_bh$elephant == 'M203') # 480
which(eles_bh$elephant == 'M217') # 502
which(eles_bh$elephant == 'M88')  # 435
which(eles_bh$elephant == 'M9')   # 382

eles_bh$number[which(eles_bh$elephant == 'M9')]   # 382, e5
eles_bh$number[which(eles_bh$elephant == 'M2')]   # 433, e7
eles_bh$number[which(eles_bh$elephant == 'M88')]  # 435, e7
eles_bh$number[which(eles_bh$elephant == 'M12')]  # 461, e8
eles_bh$number[which(eles_bh$elephant == 'M139')] # 477, e10
eles_bh$number[which(eles_bh$elephant == 'M203')] # 480, e10
eles_bh$number[which(eles_bh$elephant == 'M174')] # 484, e11
eles_bh$number[which(eles_bh$elephant == 'M190')] # 493, e14
eles_bh$number[which(eles_bh$elephant == 'M191')] # 496, e15
eles_bh$number[which(eles_bh$elephant == 'M120')] # 501, e17
eles_bh$number[which(eles_bh$elephant == 'M217')] # 502, e18

bh$unique[which(bh$e7  == 'M2')]    # 73, 382
bh$unique[which(bh$e10 == 'M139')]  # 73, 382
bh$unique[which(bh$e11 == 'M174')]  # 73, 382
bh$unique[which(bh$e14 == 'M190')]  # 73, 382
bh$unique[which(bh$e15 == 'M191')]  # 73, 382
bh$unique[which(bh$e7  == 'M88')]   # 80, 396
bh$unique[which(bh$e17 == 'M120')]  # 147 148, 721 722
bh$unique[which(bh$e18 == 'M217')]  # 147 148, 721 722
bh$unique[which(bh$e5  == 'M9')]    # 156, 759
bh$unique[which(bh$e10 == 'M203')]  # 157 158 159, 769 770 771
bh$unique[which(bh$e8  == 'M12')]   # 175, 947

View(bh[73,])               

View(d[which(d$unique == 382),])  # 3F, 3U + 9M --> MX
View(d[which(d$unique == 396),])  # B7 + M88 --> MX
View(d[which(d$unique == 721),])  # 2 AM, 9 F, 4 PM, 4 U --> MX
View(d[which(d$unique == 722),])  # same sighting as above, 15 minutes apart
View(d[which(d$unique == 759),])  # B3 + M9 --> MX
View(d[which(d$unique == 769),])  # 1 AM, 4 F, 3 PM, 2 U --> MX
View(d[which(d$unique == 770),])  # same sighting as above, 15 minutes apart
View(d[which(d$unique == 771),])  # same sighting as above, 15 minutes apart
View(d[which(d$unique == 947),])  # 1AM, 1PM, 5F, 1U --> MX

d$type[which(d$unique == 382)] <- 'MX'
d$type[which(d$unique == 396)] <- 'MX'
d$type[which(d$unique == 721)] <- 'MX'
d$type[which(d$unique == 722)] <- 'MX'
d$type[which(d$unique == 759)] <- 'MX'
d$type[which(d$unique == 769)] <- 'MX'
d$type[which(d$unique == 770)] <- 'MX'
d$type[which(d$unique == 771)] <- 'MX'
d$type[which(d$unique == 947)] <- 'MX'


bh <- d[d$type == 'BH',]
which(bh$e5 == 'M9')   # none

mo <- d[d$type == 'MO',]
which(mo$e1 == 'F102') # none

unique(d$type) # "MO" "BH" "MX" "UK"

colnames(d)[14] <- 'herd_type'

### write csv
write_delim(d, 'data_processed/motnp_encounters_22.01.13.csv', na = 'NA', col_names = T, delim = ',')

s <- d

### clear environment, leave d for mapping
rm(bh, counts, d, e, ele_nodes, eles_bh, eles_long, encounters, mo, plus, group.size, i, j)
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

#### Convert sightings data to long format for analysis ####
### import encounter data
d <- read_csv('data_processed/motnp_encounters_22.01.13.csv')
for(j in 2:ncol(d)){
  for(i in 1:nrow(d)){
    d[i,j] <- ifelse(d[i,j] == 'NA',NA,d[i,j])
  }
}
str(d)

### consider only the first sighting of each group, not the repeated measures every 15 minutes
typ1 <- d[d$typ == 1,]                # anything measured as typ ≥2 is a repeat sighting of the same individuals
typNA <- d[is.na(d$typ),]             # Did not start to repeat sightings data every 15 minutes immediately, so typ=NA before we started this
first <- rbind(typNA, typ1)           # combine those that are first sightings into a single data frame
first <- first[!is.na(first$date),]   # remove data points where no date recorded
first <- first[,c(1,3:14,28:103)]     # remove columns counting number of each demographic group -- unnecessary and inaccurate
first$encounter <- 1:nrow(first) ; first <- first[,c(90,1:89)]   # index variable of sightings
first$date <- lubridate::as_date(first$date)                     # put date column into correct format
length(unique(first$encounter))       # 1078 separate encounters
str(first)

### make long, so every row is a sighting of an individual elephant
eles_long <- gather(data = first, key = number, value = elephant, e1:e76, factor_key = TRUE)
length(unique(eles_long$encounter))   # 1078
eles_long <- eles_long[!is.na(eles_long$elephant),] # filter out NA values
length(unique(eles_long$encounter))   # 574 -- number of independent sightings containing identified individuals

### 2 incorrectly labelled elephants: no such code as N or D so N40 and D128 must be wrong. On keyboard, N next to M and D next to F --> assume that D128=F128 and N40=M40 (M40 one of the most commonly observed elephants).
eles_long$elephant[which(eles_long$elephant == 'D128')] <- 'F128'  # row 2297
eles_long$elephant[which(eles_long$elephant == 'N40')] <- 'M40'    # row 1859

### write to csv for use in plotting
write_delim(eles_long, 'data_processed/motnp_eles_long_22.01.13.csv',col_names = T, delim = ',')

### clear environment
rm(d, typ1, typNA, first, eles_long)
#### Create nodes data frame ####
### nodes data frame
ele_nodes <- read_delim('data_processed/motnp_id_22.01.13.csv', delim = ',')  # read in ID data
ele_nodes <- ele_nodes[,c(1:7,11:12)]                                         # select desired columns

### make long (zero-padded) label for each elephant in eles_long to match ele_nodes
eles_long <- read_delim('data_processed/motnp_eles_long_22.01.13.csv', delim = ',')  # read in individual sightings data
eles_long$time <- as.numeric(eles_long$time)   # fix time variable
eles_long$id_no <- eles_long$elephant          # add id_no variable to pad out
eles_long <- separate(eles_long, id_no, into = c('sex','num'), sep = 1)  # separate id_no variable for padding
eles_long$num <- sprintf("%04d", as.integer(eles_long$num))              # pad out
eles_long$id_no <- paste(eles_long$sex, eles_long$num, sep='')           # recombine
eles_long <- eles_long[,c(1:8,10,13,14,16,17,19)] # remove unneccsary columns

### make short (unpadded) label for each elephant in ele_nodes to match eles_long
ele_nodes$id  <- ele_nodes$id_no                                          # create variable
ele_nodes     <- separate(ele_nodes, id, into = c('sex','num'), sep = 1)  # separate to remove 0s
ele_nodes$num <- as.numeric(ele_nodes$num)                                # remove 0s
ele_nodes$id  <- paste(ele_nodes$sex,ele_nodes$num,sep='')                # recombine

### cleaning
ele_nodes <- ele_nodes[which(ele_nodes$id_no != 'M0175'),] # remove M0175 from nodes -- doesn't exist in sightings

### Insert column for count of sightings
ele_nodes$count <- NA
for(i in 1:nrow(ele_nodes)){
  ele_nodes$count[i] <- sum(eles_long$id_no == ele_nodes$id_no[i]) # count all instances of individual presence
}

### create column for combined age and sex in both data frames
ele_nodes$dem_class <- paste0(ifelse(ele_nodes$age_class == 'Adult','A',
                                     ifelse(ele_nodes$age_class == 'Pubescent','P',
                                            ifelse(ele_nodes$age_class == 'Juvenile', 'J','C'))),
                              ele_nodes$sex)

### write CSV file
write_delim(ele_nodes, path = 'data_processed/motnp_elenodes_22.01.13.csv', delim = ',', col_names = T)

### clear environment
rm(ele_nodes, eles_long, i, j)
################ 2) Create dyadic data frame of sightings -- very slow (~ 48 hours to run) ################
#### Create dyadic data for all sightings ####
### import data
# elephant encounters
eles <- read_delim(file = 'data_processed/motnp_eles_long_22.01.13.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_') # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]            # rearrange variables
str(eles)

# nodes
nodes <- read_delim(file = 'data_processed/motnp_elenodes_22.01.13.csv', delim = ',') # read in node data
colnames(nodes)
str(nodes)

### create group-by-individual matrix
eles_asnipe <- eles[,c(3,4,2,5)]                                       # date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)                       # make date numeric
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
which(is.na(eles_asnipe$Time))                                         # 161 698 1122 1469 1770
eles_asnipe[c(161,698,1122,1469,1770),]                                # all 1 sighting of B7+M44

eles_asnipe <- eles_asnipe[,c(5,6,3,4)]                                # create data frame to produce gbi matrix from
colnames(eles_asnipe) <- c('Date','Time','ID','Location')              # rename variables for get_gbi
eles_asnipe$ID <- as.character(eles_asnipe$ID)                         # correct data type
str(eles_asnipe)

# get_gbi generates a group-by-individual matrix. The function accepts a data.table with individual identifiers and a group column. The gbi matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')  # 0-pad dates
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value for each sighting
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter)) # unique factor for every sighting
max(eles_asnipe$group)                         # 574 -- agrees with number of different sightings for which elephants were identified
eles_asnipe <- eles_asnipe[,c(3,7)]            # create data table for gbi matrix
eles_asnipe <- data.table::setDT(eles_asnipe)  # create data table for gbi matrix
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')  # create gbi matrix

### code to convert gbi matrix format to dyadic data frame, shared by Prof Dan Franks and available from @JHart96 GitHub repository (https://github.com/JHart96/bison_examples/blob/main/examples/convert_gbi.md) -- NOTE: this step takes ~ 1.5 days to run
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric()) # create empty data frame
for (obs_id in 1:nrow(gbi_matrix)) {              # run through every sighting in gbi matrix
  for (i in which(gbi_matrix[obs_id, ] == 1)) {   # for that sighting, take every individual that was present
    for (j in 1:ncol(gbi_matrix)) {               # run through every other elephant and record if they were also present
      if (i != j) {
        # make sure node_1 < node_2: not necessary but makes things a bit easier
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
                                           obs_id = obs_id) # add a row to gbi_df data frame containing 1/0 for both individuals present or only reference one
      }
    }
  }
}
gbi_df # check structure of gbi_df
write_delim(gbi_df, 'data_processed/motnp_bayesian_allpairwiseevents_22.01.13.csv', delim = ',')
# test <- read_delim('data_processed/motnp_bayesian_allpairwiseevents_22.01.13.csv', delim = ',')  # check that this has worked because you don't want to have to run it again!

### add elephant ID numbers to assigned index factors
str(gbi_df)
gbi_id <- data.frame(id_1 = colnames(gbi_matrix), node_1 = as.numeric(1:472))
gbi_check <- left_join(x = gbi_df, y = gbi_id, by = 'node_1')
gbi_id <- data.frame(id_2 = colnames(gbi_matrix), node_2 = as.numeric(1:472))
gbi_check <- left_join(x = gbi_check, y = gbi_id, by = 'node_2')

### correct obs_id to match encounter numbers in other spreadsheets (some encounters have missing data: encounter numbers 1,2,3,5,6,8... where obs_id 1,2,3,4,5,6...). Warning: this step can take quite a lot of memory - may need to clear some things from environment
eles$obs_id <- as.integer(as.factor(eles$encounter)) ; eles <- eles[,c(1,13,2:12)] # create factor variable per encounter
gbi_encounter <- left_join(x = gbi_check, y = eles, by = 'obs_id')  # add sighting information to gbi_df
length(unique(gbi_encounter$encounter)) # 574 -- correct

### remove duplicate rows where an dyad is recording as being observed together twice during the same sighting
gbi_encounter$unique <- paste(gbi_encounter$node_1, gbi_encounter$node_2, gbi_encounter$social_event, gbi_encounter$obs_id, sep = '_') # create unique variable for every dyad and encounter
gbi_distinct <- dplyr::distinct(gbi_encounter) # 40708781 obs
colnames(gbi_distinct) # "node_1","node_2","social_event","obs_id","id_1","id_2","encounter","elephant","date","time","location","gps_s","gps_e",'herd_type","total_elephants_numeric","total_elephants_uncert","total_id_hkm","perc_id_hkm","unique"
head(gbi_distinct, 10)
gbi_distinct <- gbi_distinct[,c(1:7,9:18)]                                   # rearrange columns
colnames(gbi_distinct)[c(7,16:17)] <- c('encounter_id','total_id','perc_id') # rename columns

### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power
gbi_distinct$dyad <- paste(gbi_distinct$id_1, gbi_distinct$id_2, sep = '_') # create character variable unique to each dyad
gbi_distinct$dyad_id <- as.integer(as.factor(gbi_distinct$dyad))            # create integer variable unique to each dyad
gbi_distinct$location_id <- as.integer(as.factor(gbi_distinct$location))    # create integer variable unique to each location
gbi_distinct <- dplyr::distinct(gbi_distinct)                               # remove any duplicate sightings of same individual at same place and time
head(gbi_distinct)

### convert to Binomial model data format -- aggregate all sightings of each dyad together into a count
df_agg <- gbi_distinct %>%
  group_by(id_1, id_2) %>%  # group by dyad pair
  summarise(event_count=sum(social_event), dyad_id=cur_group_id()) %>% # count total number of times dyad were seen together
  mutate(node_1_id=as.integer(as.factor(id_1)), node_2_id=as.integer(as.factor(id_2))) # record node IDs
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
df_agg <- df_agg[,c(1:4)]                              # select desired columns
df_agg$node_1 <- as.integer(as.factor(df_agg$id_1))    # F1 = 1, F10 = 2, F100 = 3, F101 = 4, F102 = 5... U99 = N
df_agg$node_2 <- as.integer(as.factor(df_agg$id_2))+1  # add 1 so starts at 2 (otherwise F10 is 2 in id_1, and 1 in id_2)
head(df_agg,10) ; tail(df_agg,10)

### add data about nodes
colnames(nodes)
nodes <- nodes[,c(1,3:5,9,11:13)]                      # select desired columns
nodes$id_1 <- nodes$id ; nodes$id_2 <- nodes$id        # add columns for combining datasets
colnames(nodes) ; colnames(df_agg)
dyads <- left_join(x = df_agg, y = nodes, by = 'id_1') # add all information about elephant 1
colnames(dyads)[c(2,7:15)] <- c('id_2','id_pad_1','name_1','age_class_1','age_category_1','sex_1','id1_deletecolumn','count_1','dem_class_1','deletecolumn1')
dyads <- left_join(x = dyads, y = nodes, by = 'id_2')  # add all information about elephant 2
colnames(dyads)[c(1,16:24)] <- c('id_1','id_pad_2','name_2','age_class_2','age_category_2','sex_2','id2_deletecolumn','count_2','dem_class_2','deletecolumn2')
dyads <- dyads[,c(4,1,2,5,6,3,7,16,8,17,9,18,10,19,11,20,13,22,14,23)] # rearrange variables
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
                 dyads$id_1 != "M21" & dyads$id_1 != "M223" & dyads$id_1 != "M227", ] # remove individuals with incorrect data
dyads <- dyads[dyads$id_2 != "F157" & dyads$id_2 != "F158" & dyads$id_2 != "F176" & 
                 dyads$id_2 != "M125" & dyads$id_2 != "M13" & dyads$id_2 !=  "M138" & 
                 dyads$id_2 != "M21" & dyads$id_2 != "M223" & dyads$id_2 != "M227", ] # remove individuals with incorrect data
length(which(is.na(dyads$id_pad_1)))                 # 0 entries where elephants have no information
length(which(is.na(dyads$id_pad_2)))                 # 0 entries where elephants have no information

### write csv
readr::write_delim(dyads, 'data_processed/motnp_bayesian_trimmedpairwiseevents_22.01.13.csv', delim = ',') # write to file

#### Create dyadic data with MO sightings -- unused in analysis ####
### import data
# elephant encounters
eles <- read_delim(file = 'data_processed/motnp_eles_long_22.01.13.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')
unique(eles$herd_type)       # "MO" "BH" "MX" "UK"
mo <- eles[eles$herd_type == 'MO',]
mo <- mo[!is.na(mo$encounter),]
mo <- mo[,c(1,16,2,3,17,4,5,14,7,8,10,13)]
str(mo)
length(unique(mo$encounter)) # 247
sort(unique(mo$elephant))    # no females or unknowns

# nodes
nodes <- read_delim(file = 'data_processed/motnp_elenodes_22.01.13.csv', delim = ',')
nodes <- nodes[nodes$dem_class == 'AM' | nodes$dem_class == 'PM',]
str(nodes)

### create group-by-individual matrix
eles_asnipe <- mo[,c(3,4,2,5)]
eles_asnipe$Date <- as.integer(eles_asnipe$date)
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
which(is.na(eles_asnipe$Time))                                         # all filled

eles_asnipe <- eles_asnipe[,c(5,6,3,4)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
str(eles_asnipe)

# get_gbi generates a group by individual matrix. The function accepts a data.table with individual identifiers and a group column. The group by individual matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_')
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter))
max(eles_asnipe$group) # 247 -- number of different MO sightings for which elephants were identified
eles_asnipe <- eles_asnipe[,c(3,7)]
eles_asnipe <- data.table::setDT(eles_asnipe)
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')
sum(gbi_matrix) # 681

### convert gbi matrix format to dyadic data frame -- NOTE: this step takes several minutes
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1:nrow(gbi_matrix)) {
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
}
gbi_df
test <- gbi_df[gbi_df$social_event == 1,]
length(unique(test$obs_id)) # 135

write_delim(gbi_df, 'data_processed/motnp_bayesian_allpairwiseevents_maleonly_22.01.13.csv', delim = ',')
# gbi_df <- read_csv('data_processed/motnp_bayesian_allpairwiseevents_maleonly_22.01.13.csv')

### add elephant ID numbers to assigned index factors
str(gbi_df)
gbi_id <- data.frame(id_1 = colnames(gbi_matrix), node_1 = as.numeric(1:123))
gbi_check <- left_join(x = gbi_df, y = gbi_id, by = 'node_1')
gbi_id <- data.frame(id_2 = colnames(gbi_matrix), node_2 = as.numeric(1:123))
gbi_check <- left_join(x = gbi_check, y = gbi_id, by = 'node_2')

# no M10, M108, M109...
nodes[nodes$id_no == 'M0010',] # 3 observations
nodes[nodes$id_no == 'M0108',] # 1 observation
nodes[nodes$id_no == 'M0109',] # 2 observations
d <- read_csv('data_processed/motnp_eles_long_22.01.13.csv')
m10  <- d[d$elephant == 'M10' ,] # 3 observations -- not MO
m108 <- d[d$elephant == 'M108',] # 1 observation  -- not MO
m109 <- d[d$elephant == 'M109',] # 2 observations -- not MO
# missing mlaes were not seen in MO groups

### correct obs_id to match encounter numbers in other spreadsheets (some encounters have missing data: encounter numbers 1,2,3,5,6,8... where obs_id 1,2,3,4,5,6...). Warning: sometimes this step can exceed vector memory and not run.
eles$obs_id <- as.integer(as.factor(eles$encounter))
eles <- eles[,c(1,18,2:17)]
gbi_encounter <- left_join(x = gbi_check, y = eles, by = 'obs_id')
length(unique(gbi_encounter$encounter)) # 247 -- correct

### remove duplicate rows where an dyad is recording as being observed together twice during the same sighting
gbi_encounter$unique <- paste(gbi_encounter$node_1, gbi_encounter$node_2, gbi_encounter$social_event, gbi_encounter$obs_id, sep = '_')
gbi_distinct <- dplyr::distinct(gbi_encounter)
colnames(gbi_distinct) # "node_1","node_2","social_event","obs_id","id_1","id_2","encounter","elephant","date","time","location","gps_s","gps_e","total_elephants","total_id_dy","total_id_hkm","total_id_diff","perc_id_dy","perc_id_hkm","herd_type","number","elephant","location","unique"
gbi_distinct <- gbi_distinct[,c(1:14,16,19,20,23)]
colnames(gbi_distinct)[c(7,15:16)] <- c('encounter_id','total_id','perc_id')

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
  mutate(node_1_id = as.integer(as.factor(id_1)), node_2_id = as.integer(as.factor(id_2)))
length(df_agg$id_1) == cumsum(1:122)[122] # check have correct number of dyads -- number will be the (n-1)th value of the triangular number sequence in which n = total number of elephants in analysis (123). If TRUE, correct number of pairs.
head(df_agg) ; tail(df_agg)

### node_1_id and node_2_id wrong but unnecessary in these data sets -- only needed in one for all group types
df_agg <- df_agg[,c(1:4)]

### write csv
readr::write_delim(df_agg, 'data_processed/motnp_bayesian_trimmedpairwiseevents_maleonly_22.01.13.csv', delim = ',')

#### Create dyadic data with MX sightings -- unused in analysis ####
### import data
# elephant encounters
eles <- read_delim(file = 'data_processed/motnp_eles_long_22.01.13.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')
unique(eles$herd_type)       # "MO" "BH" "MX" "UK"
mx <- eles[eles$herd_type == 'MX',]
mx <- mx[!is.na(mx$encounter),]
mx <- mx[,c(1,16,2,3,17,4,5,14,7,8,10,13)]
str(mx)
length(unique(mx$encounter)) # 223
length(unique(mx$elephant))  # 451 -- very few elephants never seen in mixed group
mx_eles <- sort(unique(mx$elephant))

# nodes
nodes <- read_delim(file = 'data_processed/motnp_elenodes_22.01.13.csv', delim = ',')
mx_nodes <- sort(unique(nodes$id))
mx_nodes
mx_eles

# remove all individuals never seen in mixed groups
which(mx_nodes != mx_eles) # 60 onwards
mx <- mx[mx$elephant != 'F158',]
which(mx_nodes != mx_eles) # 62 onwards
mx <- mx[mx$elephant != 'F176',]
mx <- mx[mx$elephant != 'M125',]
mx <- mx[mx$elephant != 'M13',]
mx <- mx[mx$elephant != 'M138',]
mx <- mx[mx$elephant != 'M21',]
mx <- mx[mx$elephant != 'M223',]
mx <- mx[mx$elephant != 'M227',]
nodes <- nodes[nodes$id != 'M101',]
nodes <- nodes[nodes$id != 'M105',]
nodes <- nodes[nodes$id != 'M107',]
nodes <- nodes[nodes$id != 'M125',]
nodes <- nodes[nodes$id != 'M129',]
nodes <- nodes[nodes$id != 'M132',]
nodes <- nodes[nodes$id != 'M144',]
nodes <- nodes[nodes$id != 'M150',]
nodes <- nodes[nodes$id != 'M16',]
nodes <- nodes[nodes$id != 'M195',]
nodes <- nodes[nodes$id != 'M199',]
nodes <- nodes[nodes$id != 'M201',]
nodes <- nodes[nodes$id != 'M210',]
nodes <- nodes[nodes$id != 'M214',]
nodes <- nodes[nodes$id != 'M225',]
nodes <- nodes[nodes$id != 'M238',]
nodes <- nodes[nodes$id != 'M252',]
nodes <- nodes[nodes$id != 'M55',]
nodes <- nodes[nodes$id != 'M7',]
nodes <- nodes[nodes$id != 'M93',]
nodes <- nodes[nodes$id != 'U63',]

mx_eles <- sort(unique(mx$elephant))
mx_nodes <- sort(unique(nodes$id))
which(mx_nodes != mx_eles) # all the same individuals

str(nodes)

### create group-by-individual matrix
eles_asnipe <- mx[,c(3,4,2,5)]
eles_asnipe$Date <- as.integer(eles_asnipe$date)
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
which(is.na(eles_asnipe$Time))                                         # all filled

eles_asnipe <- eles_asnipe[,c(5,6,3,4)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
str(eles_asnipe)

# get_gbi generates a group by individual matrix. The function accepts a data.table with individual identifiers and a group column. The group by individual matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_')
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter))
max(eles_asnipe$group) # 247 -- number of different MX sightings for which elephants were identified
eles_asnipe <- eles_asnipe[,c(3,7)]
eles_asnipe <- data.table::setDT(eles_asnipe)
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')
sum(gbi_matrix) # 3097

### convert gbi matrix format to dyadic data frame -- NOTE: this step takes several hours, run overnight
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1:nrow(gbi_matrix)) {
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
}
gbi_df
test <- gbi_df[gbi_df$social_event == 1,]
length(unique(test$obs_id)) # 213...? Should this be 223?

write_delim(gbi_df, 'data_processed/motnp_bayesian_allpairwiseevents_mixed_22.01.13.csv', delim = ',')

### add elephant ID numbers to assigned index factors
str(gbi_df)
gbi_id <- data.frame(id_1 = colnames(gbi_matrix), node_1 = as.numeric(1:443))
gbi_check <- left_join(x = gbi_df, y = gbi_id, by = 'node_1')
gbi_id <- data.frame(id_2 = colnames(gbi_matrix), node_2 = as.numeric(1:443))
gbi_check <- left_join(x = gbi_check, y = gbi_id, by = 'node_2')
gbi_check[gbi_check$id_1 == 'F1',]
gbi_check[gbi_check$id_1 == 'M1',]

### correct obs_id to match encounter numbers in other spreadsheets (some encounters have missing data: encounter numbers 1,2,3,5,6,8... where obs_id 1,2,3,4,5,6...). Warning: sometimes this step can exceed vector memory and not run.
eles$obs_id <- as.integer(as.factor(eles$encounter))
eles <- eles[,c(1,18,2:17)]
gbi_encounter <- left_join(x = gbi_check, y = eles, by = 'obs_id')
length(unique(gbi_encounter$encounter)) # 223 -- correct

### remove duplicate rows where an dyad is recording as being observed together twice during the same sighting
gbi_encounter$unique <- paste(gbi_encounter$node_1, gbi_encounter$node_2, gbi_encounter$social_event, gbi_encounter$obs_id, sep = '_')
gbi_distinct <- dplyr::distinct(gbi_encounter)
colnames(gbi_distinct) # "node_1","node_2","social_event","obs_id","id_1","id_2","encounter","date","time","gps_s","gps_e","total_elephants","total_elephants_numeric","total_elephants_uncert","total_id_dy","total_id_hkm","total_id_diff","perc_id_dy","perc_id_hkm","herd_type","number","elephant","location","unique"
gbi_distinct <- gbi_distinct[,c(1:14,16,19,20,23)]
colnames(gbi_distinct)[c(7,15:16)] <- c('encounter_id','total_id','perc_id')

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
  mutate(node_1_id = as.integer(as.factor(id_1)), node_2_id = as.integer(as.factor(id_2)))
length(df_agg$id_1) == cumsum(1:442)[442] # check have correct number of dyads -- number will be the (n-1)th value of the triangular number sequence in which n = total number of elephants in analysis (443). If TRUE, correct number of pairs.
head(df_agg) ; tail(df_agg)

### check values
sum(df_agg$event_count) == sum(gbi_matrix) # not even close.... but actually should they be? gbi_matrix shows a 1 per elephant, but that will be multiplied by GroupSize-1 for every event (dyad)
summary(df_agg$event_count) # 0-27
df_agg[which(df_agg$event_count == 27),] # F52 + U17 = mother and calf -- very reasonable to be seen 27 times together
nodes[nodes$id_no == 'F0052',] # count = 39, so 27 if anything is low for that dyad
nodes[nodes$id_no == 'U0017',] # count = 39, so 27 if anything is low for that dyad
# could be 39 sightings but including typ = 2+? --> no because they had already been removed before creating count column...

df_agg[which(df_agg$id_1 == 'F52' & df_agg$id_2 == 'F60'),] # event_count = 25
nodes[nodes$id_no == 'F0060',]                              # count = 36, 25 again seems low

B7 <- eles[eles$elephant == 'F52' | eles$elephant == 'F60' | eles$elephant == 'U17' | eles$elephant == 'U21' | eles$elephant == 'F98',]
B7 <- B7[,c(1:9,11,14,15,17)]
head(B7)
length(which(B7$elephant == 'F52')) # 39
length(which(B7$elephant == 'U17')) # 39
length(which(B7$elephant == 'F60')) # 36
length(which(B7$elephant == 'U21')) # 35
length(which(B7$elephant == 'F98')) # 31

f52 <- sort(B7$encounter[which(B7$elephant == 'F52')])
u17 <- sort(B7$encounter[which(B7$elephant == 'U17')])
f52 ; u17 # U17 listed twice for encounter 487 and absent from 679, but otherwise all 39 are the same -- 27 is definitely wrong. Or not -- could be 

length(which(B7$elephant == 'F52' & B7$herd_type == 'MX')) # 28
length(which(B7$elephant == 'U17' & B7$herd_type == 'MX')) # 27 -- there it is! OK so this is likely still wrong because should have been there for the 28th sighting too, but this is why there are only 27 sightings.

length(which(B7$elephant == 'F60' & B7$herd_type == 'MX')) # 26
length(which(B7$elephant == 'U21' & B7$herd_type == 'MX')) # 26
df_agg[which(df_agg$event_count == 26),] # F60 and U21 -- perfect

f52 <- sort(B7$encounter[which(B7$elephant == 'F52' & B7$herd_type == 'MX')])
f60 <- sort(B7$encounter[which(B7$elephant == 'F60' & B7$herd_type == 'MX')])
f52 # 129 141 146 157 165 180 232 396 422 468 470 480 482 485 515 534 539 629 644 671 679 738 835 857 866 + 72 101 103 = 25 + 3 non-matching
f60 # 129 141 146 157 165 180 232 396 422 468 470 480 482 485 515 534 539 629 644 671 679 738 835 857 866 + 128        = 25 + 1 non-matching
# counts are correct

### write csv
readr::write_delim(df_agg, 'data_processed/motnp_bayesian_trimmedpairwiseevents_mixed_22.01.13.csv', delim = ',')

#### Create dyadic data with BH sightings -- unused in analysis ####
### import data
# elephant encounters
eles <- read_delim(file = 'data_processed/motnp_eles_long_22.01.13.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')
bh <- eles[eles$herd_type == 'BH',]
bh <- bh[!is.na(bh$encounter),]
bh <- bh[,c(1,16,2,3,17,4,5,14,7,8,10,13)]
str(bh)
length(unique(bh$encounter)) # 97
length(unique(bh$elephant))  # 117
bh_eles <- sort(unique(bh$elephant))

s <- read_csv('data_processed/motnp_encounters_22.01.13.csv')
typNA <- s[is.na(s$typ),] ; typ1 <- s[s$typ == 1,]
s <- rbind(typNA, typ1) ; s <- s[!is.na(s$date),] ; s <- s[!is.na(s$e1),]
length(which(s$herd_type == 'BH')) # 97

# nodes
nodes <- read_delim(file = 'data_processed/motnp_elenodes_22.01.13.csv', delim = ',')
nodes <- nodes[nodes$dem_class != 'AM',]
bh_nodes <- sort(unique(nodes$id))
bh_nodes # all elephants
bh_eles  # not many elephants only seen in breeding herds

### remove all individuals never seen in breeding herds
# remove from nodes: "F100" "F101"  "F103" "F105" "F106" "F108" "F111" "F112" "F113" "F114" "F115" "F116" "F117" "F118" "F119" "F121" "F123" "F125" "F126" "F127" "F128"  "F130" "F131" "F134" "F135" "F138" "F139" "F14"  "F140" "F143" "F144" "F145" "F146" "F148" "F149" "F15"  "F150" "F151" "F152" "F154" "F156"  "F23"  "F24"  "F27"  "F31"  "F32"  "F35"  "F37" "F39"  "F4"  "F42" "F45" "F48"  "F49"  "F5"  "F53"  "F54"  "F55"  "F56"  "F58"  "F59"  "F61"  "F63"  "F64"  "F65"  "F66"  "F68"  "F70"  "F71"  "F72"  "F73"  "F78"  "F79" "F80"  "F81"  "F82"  "F83"  "F87"  "F88"  "F89"  "F90"  "F91"  "F92" "F93"  "F95" "F97" "F99"  "M102" "M104" "M109" "M110" "M114" "M117" "M121" "M122" "M128" "M131" "M135" "M136" "M147" "M15"  "M153" "M155" "M157" "M158" "M159" "M161" "M162" "M164" "M165" "M166" "M167" "M169" "M17"  "M173" "M178" "M181" "M183" "M185" "M186" "M187" "M189" "M195" "M198" "M199" "M202" "M208" "M210" "M212" "M213" "M215" "M216" "M218" "M22"  "M220" "M221" "M226" "M228"  "M235"  "M241"  "M246" "M247"  "M25"  "M255"  "M38"  "M64" "M70"  "M77"  "M89"  "M91"  "M96"  "M97"   "U10"  "U11" "U13"  "U14"  "U16"  "U18"  "U19"  "U2"   "U20"  "U22"  "U24"  "U27"  "U28"  "U29"  "U3"   "U30"  "U32"  "U33"  "U34"  "U35"  "U36"  "U37"  "U38"  "U39"  "U4"   "U40"  "U41"   "U44"  "U45"  "U46"  "U47"  "U49"  "U50"  "U51"  "U52"  "U53"  "U56"  "U57"  "U59"  "U6"   "U64"  "U65"  "U66"  "U67"  "U7"  
# remove from bh: "F157" 

bh <- bh[bh$elephant != 'F157',]
nodes <- nodes[nodes$id != "F100" & nodes$id != "F101" & nodes$id != "F103" & nodes$id != "F105" & nodes$id != "F106" &
           nodes$id != "F108" & nodes$id != "F111" & nodes$id != "F112" & nodes$id != "F113" & nodes$id != "F114" &
           nodes$id != "F115" & nodes$id != "F116" & nodes$id != "F117" & nodes$id != "F118" & nodes$id != "F119" &
           nodes$id != "F121" & nodes$id != "F123" & nodes$id != "F125" & nodes$id != "F126" & nodes$id != "F127" &
           nodes$id != "F128" & nodes$id != "F130" & nodes$id != "F131" & nodes$id != "F134" & nodes$id != "F135" &
           nodes$id != "F138" & nodes$id != "F139" & nodes$id != "F14"  & nodes$id != "F140" & nodes$id != "F143" &
           nodes$id != "F144" & nodes$id != "F145" & nodes$id != "F146" & nodes$id != "F148" & nodes$id != "F149" &
           nodes$id != "F15"  & nodes$id != "F150" & nodes$id != "F151" & nodes$id != "F152" & nodes$id != "F154" &
           nodes$id != "F156" & nodes$id != "F23"  & nodes$id !=  "F24" & nodes$id != "F27"  & nodes$id != "F31"  &
           nodes$id != "F32"  & nodes$id != "F35"  & nodes$id !=  "F37" & nodes$id != "F39"  & nodes$id !=  "F4"  &
           nodes$id != "F42"  & nodes$id != "F45"  & nodes$id != "F48"  & nodes$id != "F49"  & nodes$id !=  "F5"  &
           nodes$id != "F53"  & nodes$id != "F54"  & nodes$id != "F55"  & nodes$id != "F56"  & nodes$id != "F58"  &
           nodes$id != "F59"  & nodes$id != "F61"  & nodes$id != "F63"  & nodes$id != "F64"  & nodes$id != "F65"  &
           nodes$id != "F66"  & nodes$id != "F68"  & nodes$id != "F70"  & nodes$id != "F71"  & nodes$id != "F72"  &
           nodes$id != "F73"  & nodes$id != "F78"  & nodes$id != "F79"  & nodes$id != "F80"  & nodes$id != "F81"  &
           nodes$id != "F82"  & nodes$id != "F83"  & nodes$id != "F87"  & nodes$id != "F88"  & nodes$id != "F89"  &
           nodes$id != "F90"  & nodes$id != "F91"  & nodes$id != "F92"  & nodes$id != "F93"  & nodes$id != "F95"  &
           nodes$id != "F97"  & nodes$id != "F99"  & nodes$id != "M102" & nodes$id != "M104" & nodes$id != "M109" &
           nodes$id != "M110" & nodes$id != "M114" & nodes$id != "M117" & nodes$id != "M121" & nodes$id != "M122" &
           nodes$id != "M128" & nodes$id != "M131" & nodes$id != "M135" & nodes$id != "M136" & nodes$id != "M147" &
           nodes$id != "M15"  & nodes$id != "M153" & nodes$id != "M155" & nodes$id != "M157" & nodes$id != "M158" &
           nodes$id != "M159" & nodes$id != "M161" & nodes$id != "M162" & nodes$id != "M164" & nodes$id != "M165" &
           nodes$id != "M166" & nodes$id != "M167" & nodes$id != "M169" & nodes$id != "M17"  & nodes$id != "M173" &
           nodes$id != "M178" & nodes$id != "M181" & nodes$id != "M183" & nodes$id != "M185" & nodes$id != "M186" &
           nodes$id != "M187" & nodes$id != "M189" & nodes$id != "M195" & nodes$id != "M198" & nodes$id != "M199" &
           nodes$id != "M202" & nodes$id != "M208" & nodes$id != "M210" & nodes$id != "M212" & nodes$id != "M213" &
           nodes$id != "M215" & nodes$id != "M216" & nodes$id != "M218" & nodes$id != "M22"  & nodes$id != "M220" &
           nodes$id != "M221" & nodes$id != "M226" & nodes$id != "M228" & nodes$id != "M235" & nodes$id != "M241" &
           nodes$id != "M246" & nodes$id != "M247" & nodes$id != "M25"  & nodes$id != "M255" & nodes$id != "M38"  &
           nodes$id != "M64"  & nodes$id != "M70"  & nodes$id != "M77"  & nodes$id != "M89"  & nodes$id != "M91"  &
           nodes$id != "M96"  & nodes$id != "M97"  & nodes$id != "U10"  & nodes$id != "U11"  & nodes$id != "U13"  &
           nodes$id != "U14"  & nodes$id != "U16"  & nodes$id != "U18"  & nodes$id != "U19"  & nodes$id != "U2"   &
           nodes$id != "U20"  & nodes$id != "U22"  & nodes$id != "U24"  & nodes$id != "U27"  & nodes$id != "U28"  &
           nodes$id != "U29"  & nodes$id != "U3"   & nodes$id != "U30"  & nodes$id != "U32"  & nodes$id != "U33"  &
           nodes$id != "U34"  & nodes$id != "U35"  & nodes$id != "U36"  & nodes$id != "U37"  & nodes$id != "U38"  &
           nodes$id != "U39"  & nodes$id != "U4"   & nodes$id != "U40"  & nodes$id != "U41"  & nodes$id !=  "U44" &
           nodes$id != "U45"  & nodes$id != "U46"  & nodes$id != "U47"  & nodes$id != "U49"  & nodes$id != "U50"  &
           nodes$id != "U51"  & nodes$id != "U52"  & nodes$id != "U53"  & nodes$id != "U56"  & nodes$id != "U57"  &
           nodes$id != "U59"  & nodes$id != "U6"   & nodes$id !=  "U64" & nodes$id != "U65"  & nodes$id != "U66"  &
           nodes$id != "U67"  & nodes$id != "U7",]
bh_eles <- sort(unique(bh$elephant))
bh_nodes <- sort(unique(nodes$id))
which(bh_nodes != bh_eles) # all the same individuals

length(unique(bh$encounter)) # 96 -- removing F157 (uncertain elephant) removed a unique sighting
# bh[which(bh$elephant == 'F157'),] # encounter = 587
# bh[which(bh$encounter == '587'),] # elephant = F157
str(nodes)

### create group-by-individual matrix
eles_asnipe <- bh[,c(3,4,2,5)]
eles_asnipe$Date <- as.integer(eles_asnipe$date)
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)                        # convert time values to seconds through day
which(is.na(eles_asnipe$Time))                                         # all filled

eles_asnipe <- eles_asnipe[,c(5,6,3,4)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
str(eles_asnipe)

# get_gbi generates a group by individual matrix. The function accepts a data.table with individual identifiers and a group column. The group by individual matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_')
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter))
max(eles_asnipe$group) # 96 -- number of different BH sightings for which elephants were identified (after excluding F157)

eles_asnipe <- eles_asnipe[,c(3,7)]
eles_asnipe <- data.table::setDT(eles_asnipe)
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')
sum(gbi_matrix) # 443

### convert gbi matrix format to dyadic data frame -- NOTE: this step takes several hours, run overnight
gbi_df <- data.frame(node_1 = numeric(), node_2 = numeric(), social_event = numeric(), obs_id = numeric())
for (obs_id in 1:nrow(gbi_matrix)) {
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
}
gbi_df
test <- gbi_df[gbi_df$social_event == 1,]
length(unique(test$obs_id)) # 92

write_delim(gbi_df, 'data_processed/motnp_bayesian_allpairwiseevents_breeding_22.01.13.csv', delim = ',')

### add elephant ID numbers to assigned index factors
str(gbi_df)
gbi_id <- data.frame(id_1 = colnames(gbi_matrix), node_1 = as.numeric(1:116))
gbi_check <- left_join(x = gbi_df, y = gbi_id, by = 'node_1')
gbi_id <- data.frame(id_2 = colnames(gbi_matrix), node_2 = as.numeric(1:116))
gbi_check <- left_join(x = gbi_check, y = gbi_id, by = 'node_2')
gbi_check[gbi_check$id_1 == 'F1',]  # looks fine
gbi_check[gbi_check$id_1 == 'M1',]  # nothing there -- good
gbi_check[gbi_check$id_1 == 'M7',]  # looks fine

### correct obs_id to match encounter numbers in other spreadsheets (some encounters have missing data: encounter numbers 1,2,3,5,6,8... where obs_id 1,2,3,4,5,6...). Warning: sometimes this step can exceed vector memory and not run.
eles$obs_id <- as.integer(as.factor(eles$encounter))
eles <- eles[,c(1,18,2:17)]
gbi_encounter <- left_join(x = gbi_check, y = eles, by = 'obs_id')
length(unique(gbi_encounter$encounter)) # 96 -- correct

### remove duplicate rows where an dyad is recording as being observed together twice during the same sighting
gbi_encounter$unique <- paste(gbi_encounter$node_1, gbi_encounter$node_2, gbi_encounter$social_event, gbi_encounter$obs_id, sep = '_')
gbi_distinct <- dplyr::distinct(gbi_encounter)
colnames(gbi_distinct) # "node_1","node_2","social_event","obs_id","id_1","id_2","encounter","date","time","gps_s","gps_e","total_elephants","total_elephants_numeric","total_elephants_uncert","total_id_dy","total_id_hkm","total_id_diff","perc_id_dy","perc_id_hkm","herd_type","number","elephant","location","unique"
gbi_distinct <- gbi_distinct[,c(1:14,16,19,20,23)]
colnames(gbi_distinct)[c(7,15:16)] <- c('encounter_id','total_id','perc_id')

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
  mutate(node_1_id = as.integer(as.factor(id_1)), node_2_id = as.integer(as.factor(id_2)))
head(df_agg) ; tail(df_agg)

### node_1_id and node_2_id wrong but unnecessary in these data sets -- only needed in one for all group types
df_agg <- df_agg[,c(1:4)]

### check values
summary(df_agg$event_count) # 0-18
df_agg[which(df_agg$id_1 == 'F52' & df_agg$id_2 == 'F60'),] # event_count = 10
df_agg[which(df_agg$event_count == 10),]                    # F52-F60 and F60-U17 -- makes sense

B7 <- eles[eles$elephant == 'F52' | eles$elephant == 'F60' | eles$elephant == 'U17' | eles$elephant == 'U21' | eles$elephant == 'F98',]
B7 <- B7[,c(1:9,11,14,15,17)]
head(B7)

length(which(B7$elephant == 'F52' & B7$herd_type == 'BH')) # 11
length(which(B7$elephant == 'U17' & B7$herd_type == 'BH')) # 12 -- one too many? was one too few before so likely this is the issue

length(which(B7$elephant == 'F60' & B7$herd_type == 'BH')) # 10 -- matches number of sightings of F52_F60 dyad
length(which(B7$elephant == 'U21' & B7$herd_type == 'BH')) # 9
df_agg[which(df_agg$event_count == 9),]   # F60 and U21, F52 and U21, U17 and U21 -- perfect

f52 <- sort(B7$encounter[which(B7$elephant == 'F52' & B7$herd_type == 'BH')])
f60 <- sort(B7$encounter[which(B7$elephant == 'F60' & B7$herd_type == 'BH')])
f52 # 106 158 163 244 356 487 532 553 557 615 778
f60 # 106 158  -  244 356 487 532 553 557 615 778

### write csv
readr::write_delim(df_agg, 'data_processed/motnp_bayesian_trimmedpairwiseevents_breeding_22.01.13.csv', delim = ',')

#### Join dyadic data sets together into a single data frame ####
all <- read_csv('data_processed/motnp_6thJan2022/motnp_bayesian_trimmedpairwiseevents_22.01.10.csv')
# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(all$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = all$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
all$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
all$dem_class_1 <- as.factor(all$dem_class_1)
all$dem_class_2 <- as.factor(all$dem_class_2)
all$dem_type <- as.factor(paste(all$dem_class_1, all$dem_class_2, sep = '_'))
all$age_type <- as.factor(paste(all$age_class_1, all$age_class_2, sep = '_'))
all$sex_type <- as.factor(paste(all$sex_1, all$sex_2, sep = '_'))
str(all)

mo <- read_csv('data_processed/motnp_bayesian_trimmedpairwiseevents_maleonly_22.01.13.csv')
str(mo)

mx <- read_csv('data_processed/motnp_bayesian_trimmedpairwiseevents_mixed_22.01.13.csv')
str(mx)

bh <- read_csv('data_processed/motnp_bayesian_trimmedpairwiseevents_breeding_22.01.13.csv')
str(bh)

all$dyad <- paste(all$id_1, all$id_2, sep = '_')
mo$dyad <- paste(mo$id_1, mo$id_2, sep = '_')
mx$dyad <- paste(mx$id_1, mx$id_2, sep = '_')
bh$dyad <- paste(bh$id_1, bh$id_2, sep = '_')

combine <- left_join(all, mo, by = 'dyad')
colnames(combine)
combine <- combine[,c(1:6,27,7:24)]
head(combine[combine$dem_type == 'AM_AM',])
colnames(combine)[c(1:7)] <- c('dyad_id','id_1','id_2','node_1','node_2','all_events','mo_events')

combine <- left_join(combine, mx, by = 'dyad')
colnames(combine)
combine <- combine[,c(1:7,28,8:25)]
head(combine[combine$dem_type == 'AF_AM',])
colnames(combine)[c(1:5,8)] <- c('dyad_id','id_1','id_2','node_1','node_2','mx_events')

combine <- left_join(combine, bh, by = 'dyad')
colnames(combine)
combine <- combine[,c(1:8,29,9:26)]
colnames(combine)[c(1:3,9:11)] <- c('dyad_id','id_1','id_2','bh_events','id_no_1','id_no_2')
head(combine)

### write to csv
write_delim(combine, 'data_processed/motnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',')
