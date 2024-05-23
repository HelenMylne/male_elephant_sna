######## Information ####
# Analysis script for paper on statistical methods for social network analysis using gambit-of-the-group data on male elephant population of Mosi-Oa-Tunya National Park, Zambia.
# Data collected: 19th May 2016-21st December 2017.
# Collected by: Mr David Youldon, Mr Dabwiso Sakala, Miss Helen Mylne and African Lion and Environmental Research Trust during this time
# Data supplied by: African Lion and Environmental Research Trust (ALERT)

######## Set up ####
# ## import packages needed
# install.packages(cmdstanr, lib = '../packages/')
# install.packages(tidyverse, lib = '../packages')
# install.packages(lubridate, lib = '../packages')
# install.packages(janitor, lib = '../packages')
# install.packages(hms, lib = '../packages')
# install.packages(readxl, lib = '../packages')
# install.packages(data.table, lib = '../packages')
# install.packages(spatsoc, lib = '../packages')
# install.packages(igraph, lib = '../packages/')

## load packages needed
library(cmdstanr, lib.loc = '../packages/')  # run models
library(tidyverse, lib.loc = '../packages')  # data manipulation
library(lubridate, lib.loc = '../packages')  # sort date columns out
library(janitor, lib.loc = '../packages')    # clean up data frame names
library(hms, lib.loc = '../packages')        # sort time columns out
library(readxl, lib.loc = '../packages')     # read in .xlxs files
library(data.table, lib.loc = '../packages') # convert data frames to data tables
library(spatsoc, lib.loc = '../packages')    # create social network data
library(igraph, lib.loc = '../packages/')    # extract eigenvector centralities

### set stan path
set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

### set seed
set.seed(12345)

# define plot theme
theme_set(theme_bw(base_size = 12))

#---------------- 1) Clean raw data                           ##########
######## IDs                                                  ####
### import data
id.raw <- readxl::read_excel("../data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 2) # read in raw data
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
ages <- readxl::read_excel("../data_raw/Raw_ALERT_ElephantAges_Saunders211019.xlsx", sheet = 1)  # read in age data
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
write_delim(id, "../data_processed/step1_dataprocessing/motnp_id.csv", delim = ",", col_names = T)

### Clear Environment
rm(ages, days.numbers, id, id.raw, names, nkash, days, i, months, weeks)

######## Encounters                                           ####
### import data
encounters <- readxl::read_excel("../data_raw/Raw_ALERT_ElephantDatabase_Youldon210811.xlsx", sheet = 3) # read in raw data
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

### write initial csv -- to be overwritten later on following generation of long format data
write_delim(d, '../data_processed/step1_dataprocessing/motnp_encounters.csv', na = 'NA', col_names = T, delim = ',')

### clear environment, leave d for mapping
rm(counts, e, plus, group.size, i, j)

######## Convert sightings data to long format for analysis   ####
### import encounter data
for(j in 2:ncol(d)){
  for(i in 1:nrow(d)){
    d[i,j] <- ifelse(d[i,j] == 'NA',NA,d[i,j]) # convert character "NA" to actual NA
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
max(first$date) - min(first$date)     # 561 days
str(first)

### make long, so every row is a sighting of an individual elephant
eles_long <- gather(data = first, key = number, value = elephant, e1:e76, factor_key = TRUE)
max(eles_long$date) - min(eles_long$date) # 561 days
eles_na <- eles_long[is.na(eles_long$elephant),] # filter out rows with no identified elephant
eles_long <- eles_long[!is.na(eles_long$elephant),] # filter out rows with no identified elephant
max(eles_long$date) - min(eles_long$date) # 504 days -- last 2 months no elephants identified
length(unique(eles_long$encounter))       # 574 -- number of independent sightings containing identified individuals

### add column for total number of males per sighting
eles_long$total_males <- NA
for(i in unique(eles_long$encounter)){
  x <- eles_long[eles_long$encounter == i,]
  eles_long$total_males[which(eles_long$encounter == i)] <- nrow(x)
}
eles_long <- eles_long[,c(1:7,17,8:16)]
total_eles <- eles_long[,c('encounter','total_males')] %>% distinct()
mean(total_eles$total_males) ; sd(total_eles$total_males)

min(eles_na$date)
max(eles_na$date)
min(eles_long$date)
max(eles_long$date)

sightings_na <- encounters[encounters$Date > '2017-10-05',]

### 2 incorrectly labelled elephants: no such code as N or D so N40 and D128 must be wrong. On keyboard, N next to M and D next to F --> assume that D128=F128 and N40=M40 (M40 one of the most commonly observed elephants).
eles_long$elephant[which(eles_long$elephant == 'D128')] <- 'F128'  # row 2297
eles_long$elephant[which(eles_long$elephant == 'N40')] <- 'M40'    # row 1859

### write to csv for later use
write_delim(eles_long, '../data_processed/step1_dataprocessing/motnp_eles_long.csv',col_names = T, delim = ',')

### clear environment
rm(d, typ1, typNA, first, i, j)

######## Create nodes data frame                              ####
### nodes data frame
ele_nodes <- read_delim('../data_processed/step1_dataprocessing/motnp_id.csv', delim = ',')  # read in ID data
ele_nodes <- ele_nodes[,c(1:7,11:12)]                                         # select desired columns

### make long (zero-padded) label for each elephant in eles_long to match ele_nodes
eles_long$time <- as.numeric(eles_long$time)   # fix time variable
eles_long$id_no <- eles_long$elephant          # add id_no variable to pad out
eles_long <- separate(eles_long, id_no, into = c('sex','num'), sep = 1)  # separate id_no variable for padding
eles_long$num <- sprintf("%04d", as.integer(eles_long$num))              # pad out
eles_long$id_no <- paste(eles_long$sex, eles_long$num, sep='')           # recombine
eles_long <- eles_long[,c(1:8,10,13,14,16,17,19)] # remove unneccesary columns

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
write_delim(ele_nodes, path = '../data_processed/step1_dataprocessing/motnp_elenodes.csv', delim = ',', col_names = T)

######## Correct herd type in encounter data                  ####
d <- read_csv('../data_processed/step1_dataprocessing/motnp_encounters.csv')
for(j in 2:ncol(d)){
  for(i in 1:nrow(d)){
    d[i,j] <- ifelse(d[i,j] == 'NA',NA,d[i,j]) # convert character "NA" to actual NA
  }
}
str(d)

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

#View(d[c(50,88,276,277,628,736,1006),]) # still an issue with 50 -- 1 elephant, but 8 identified, suggests that these may have got scrambled together somehow, but nothing to be done about it

bh <- d[d$type == 'BH',]
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

### overwrite csv with corrected herd types
write_delim(d, '../data_processed/step1_dataprocessing/motnp_encounters.csv', na = 'NA', col_names = T, delim = ',')

### clear environment
rm(list = ls()) ; gc()

#---------------- 2) Create dyadic data frame of sightings    ##########
######## Import data                                          ####
## elephant encounters
eles <- read_delim(file = '../data_processed/step1_dataprocessing/motnp_eles_long.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_') # make single variable for unique locations
eles <- eles[,c(1,16,2,3,17,4,5,14,7,8,10,13)]            # rearrange variables
str(eles)

## nodes
nodes <- read_delim(file = '../data_processed/step1_dataprocessing/motnp_elenodes.csv', delim = ',') # read in node data
colnames(nodes)
str(nodes)

######## Create group-by-individual matrix                    ####
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

## get_gbi generates a group-by-individual matrix. The function accepts a data.table with individual identifiers and a group column. The gbi matrix can then be used to build a network using asnipe::get_network.
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')  # 0-pad dates
eles_asnipe$encounter <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_') # unique value for each sighting
eles_asnipe$group <- as.integer(as.factor(eles_asnipe$encounter)) # unique factor for every sighting
max(eles_asnipe$group)                         # 574 -- agrees with number of different sightings for which elephants were identified
eles_asnipe <- eles_asnipe[,c(3,7)]            # create data table for gbi matrix
eles_asnipe <- data.table::setDT(eles_asnipe)  # create data table for gbi matrix
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe, group = 'group', id = 'ID')  # create gbi matrix

######## Convert to dyadic dataframe (warning -- slow)        ####
## code to convert gbi matrix format to dyadic data frame, available from @JHart96 GitHub repository (https://github.com/JHart96/bison_examples/blob/main/examples/convert_gbi.md) -- NOTE: this step can take several hours to run
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
  if(obs_id %% 10 == 0) {print(obs_id)}
  if(obs_id %% 10 == 0) {print(Sys.time())}
}
gbi_df # check structure of gbi_df
write_delim(gbi_df, '../data_processed/step1_dataprocessing/motnp_bayesian_bernoullipairwiseevents.csv', delim = ',')

######## Clean up dyadic data frame                           ####
## add elephant ID numbers to assigned index factors
str(gbi_df)
gbi_id <- data.frame(id_1 = colnames(gbi_matrix), node_1 = as.numeric(1:472))
gbi_check <- left_join(x = gbi_df, y = gbi_id, by = 'node_1')
gbi_id <- data.frame(id_2 = colnames(gbi_matrix), node_2 = as.numeric(1:472))
gbi_check <- left_join(x = gbi_check, y = gbi_id, by = 'node_2')

## correct obs_id to match encounter numbers in other spreadsheets (some encounters have missing data: encounter numbers 1,2,3,5,6,8... where obs_id 1,2,3,4,5,6...). Warning: this step can take quite a lot of memory - may need to clear some things from environment
eles$obs_id <- as.integer(as.factor(eles$encounter)) ; eles <- eles[,c(1,13,2:12)] # create factor variable per encounter
gbi_encounter <- left_join(x = gbi_check, y = eles, by = 'obs_id')  # add sighting information to gbi_df
length(unique(gbi_encounter$encounter)) # 574 -- correct

## remove duplicate rows where an dyad is recorded as being observed together twice during the same sighting
gbi_encounter$unique <- paste(gbi_encounter$node_1, gbi_encounter$node_2, gbi_encounter$social_event, gbi_encounter$obs_id, sep = '_') # create unique variable for every dyad and encounter
gbi_distinct <- dplyr::distinct(gbi_encounter) # 40708781 obs
colnames(gbi_distinct) # "node_1","node_2","social_event","obs_id","id_1","id_2","encounter","elephant","date","time","location","gps_s","gps_e",'herd_type","total_elephants_numeric","total_elephants_uncert","total_id_hkm","perc_id_hkm","unique"
head(gbi_distinct, 10)
gbi_distinct <- gbi_distinct[,c(1:7,9:18)]                                   # rearrange columns
colnames(gbi_distinct)[c(7,16:17)] <- c('encounter_id','total_id','perc_id') # rename columns

## convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power
gbi_distinct$dyad <- paste(gbi_distinct$id_1, gbi_distinct$id_2, sep = '_') # create character variable unique to each dyad
gbi_distinct$dyad_id <- as.integer(as.factor(gbi_distinct$dyad))            # create integer variable unique to each dyad
gbi_distinct$location_id <- as.integer(as.factor(gbi_distinct$location))    # create integer variable unique to each location
gbi_distinct <- dplyr::distinct(gbi_distinct)                               # remove any duplicate sightings of same individual at same place and time
head(gbi_distinct)

######## Aggregate sightings to gambit-of-the-group format    ####
## convert to Binomial model data format
df_agg <- gbi_distinct %>%
  group_by(id_1, id_2) %>%  # group by dyad pair
  summarise(event_count=sum(social_event), dyad_id=cur_group_id()) %>% # count total number of times dyad were seen together
  mutate(node_1_id=as.integer(as.factor(id_1)), node_2_id=as.integer(as.factor(id_2))) # record node IDs
length(df_agg$id_1) == cumsum(1:471)[471] # check have correct number of dyads -- number will be the (n-1)th value of the triangular number sequence in which n = total number of elephants in analysis (472). If TRUE, correct number of pairs.
head(df_agg) ; tail(df_agg)

## correct values in node_1_id and node_2_id using factor values.
df_agg <- df_agg[,c(1:4)]                              # select desired columns
df_agg$node_1 <- as.integer(as.factor(df_agg$id_1))    # F1 = 1, F10 = 2, F100 = 3, F101 = 4, F102 = 5... U99 = N
df_agg$node_2 <- as.integer(as.factor(df_agg$id_2))+1  # add 1 so starts at 2 (otherwise F10 is 2 in id_1, and 1 in id_2)
head(df_agg,10) ; tail(df_agg,10)

######## Add individual data and clean up data frame          ####
## add data about nodes
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

## remove any elephants from whom their is disagreement in the different data frames regarding their names or ID numbers
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

## correct sex_1, which has loaded in as a logical vector not a character/factor
unique(dyads$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = dyads$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
dyads$sex_1 <- as.character(sex_1$sex) ; rm(sex_1) ; gc()
str(dyads)  # sex_1 still comes up as logical, but now contains the right levels

## create variable for age difference
unique(dyads$age_category_1) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA --> convert anything <5 to cat 1, and anything 5-10 to cat 2
dyads$age_cat_id_1 <- ifelse(dyads$age_category_1 == '0-3', 1,
                                 ifelse(dyads$age_category_1 == '3-4', 1,
                                        ifelse(dyads$age_category_1 == '4-5', 1,
                                               ifelse(dyads$age_category_1 == '5-6', 2,
                                                      ifelse(dyads$age_category_1 == '6-7', 2,
                                                             ifelse(dyads$age_category_1 == '7-8', 2,
                                                                    ifelse(dyads$age_category_1 == '8-9', 2,
                                                                           ifelse(dyads$age_category_1 == '9-10', 2,
                                                                                  ifelse(dyads$age_category_1 == '10-15', 3,
                                                                                         ifelse(dyads$age_category_1 == '15-19', 4,
                                                                                                ifelse(dyads$age_category_1 == '20-25', 5,
                                                                                                       ifelse(dyads$age_category_1 == '20-35', 5,
                                                                                                              ifelse(dyads$age_category_1 == '25-40', 6,
                                                                                                                     ifelse(dyads$age_category_1 == '35-50', 6,
                                                                                                                            ifelse(dyads$age_category_1 == '40+', 7,
                                                                                                                                   ifelse(dyads$age_category_1 == '50+', 7, dyads$age_category_1))))))))))))))))
dyads$age_cat_id_1[which(is.na(dyads$age_cat_id_1))] <- 1 # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1

dyads$age_class_1 <- ifelse(dyads$age_cat_id_1 == 1, 'Calf',
                                ifelse(dyads$age_cat_id_1 == 2, 'Juvenile',
                                       ifelse(dyads$age_cat_id_1 > 4, 'Adult','Pubescent')))

unique(dyads$age_category_2) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
dyads$age_cat_id_2 <- ifelse(dyads$age_category_2 == '0-3', 1,
                                 ifelse(dyads$age_category_2 == '3-4', 1,
                                        ifelse(dyads$age_category_2 == '4-5', 1,
                                               ifelse(dyads$age_category_2 == '5-6', 2,
                                                      ifelse(dyads$age_category_2 == '6-7', 2,
                                                             ifelse(dyads$age_category_2 == '7-8', 2,
                                                                    ifelse(dyads$age_category_2 == '8-9', 2,
                                                                           ifelse(dyads$age_category_2 == '9-10', 2,
                                                                                  ifelse(dyads$age_category_2 == '10-15', 3,
                                                                                         ifelse(dyads$age_category_2 == '15-19', 4,
                                                                                                ifelse(dyads$age_category_2 == '20-25', 5,
                                                                                                       ifelse(dyads$age_category_2 == '20-35', 5,
                                                                                                              ifelse(dyads$age_category_2 == '25-40', 6,
                                                                                                                     ifelse(dyads$age_category_2 == '35-50', 6,
                                                                                                                            ifelse(dyads$age_category_2 == '40+', 7,
                                                                                                                                   ifelse(dyads$age_category_2 == '50+', 7, dyads$age_category_2))))))))))))))))
dyads$age_cat_id_2[which(is.na(dyads$age_cat_id_2))] <- 1   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1

dyads$age_class_2 <- ifelse(dyads$age_cat_id_2 == 1, 'Calf',
                                ifelse(dyads$age_cat_id_2 == 2, 'Juvenile',
                                       ifelse(dyads$age_cat_id_2 > 4, 'Adult','Pubescent')))

## correct dem_class with corrected age classes
dyads$dem_class_1 <- ifelse(dyads$age_class_1 == 'Adult', paste('A',dyads$sex_1, sep = ''),
                                ifelse(dyads$age_class_1 == 'Pubescent', paste('P',dyads$sex_1, sep = ''),
                                       ifelse(dyads$age_class_1 == 'Juvenile', paste('J',dyads$sex_1, sep = ''),
                                              paste('C',dyads$sex_1, sep = ''))))
dyads$dem_class_2 <- ifelse(dyads$age_class_2 == 'Adult', paste('A',dyads$sex_2, sep = ''),
                                ifelse(dyads$age_class_2 == 'Pubescent', paste('P',dyads$sex_2, sep = ''),
                                       ifelse(dyads$age_class_2 == 'Juvenile', paste('J',dyads$sex_2, sep = ''),
                                              paste('C',dyads$sex_2, sep = ''))))

## combined dem_class of dyad
dyads$age_class_id_1 <- ifelse(dyads$age_class_1 == 'Adult',4,
                                   ifelse(dyads$age_class_1 == 'Pubescent',3,
                                          ifelse(dyads$age_class_1 == 'Juvenile',2,1)))
dyads$age_class_id_2 <- ifelse(dyads$age_class_2 == 'Adult',4,
                                   ifelse(dyads$age_class_2 == 'Pubescent',3,
                                          ifelse(dyads$age_class_2 == 'Juvenile',2,1)))
dyads$dem_type <- ifelse(dyads$age_class_id_1 > dyads$age_class_id_2,
                             paste(dyads$dem_class_1, dyads$dem_class_2, sep = '_'), # when 1 is older: dc1_dc2
                             ifelse(dyads$age_class_id_1 < dyads$age_class_id_2,
                                    paste(dyads$dem_class_2, dyads$dem_class_1, sep = '_'), # when 2 is older: dc2_dc1
                                    ifelse(dyads$sex_1 == 'F',                                  # when age1 = age2...
                                           paste(dyads$dem_class_1, dyads$dem_class_2, sep = '_'), # ...when 1 is F: dc1_dc2
                                           ifelse(dyads$sex_2 == 'F',
                                                  paste(dyads$dem_class_2, dyads$dem_class_1, sep = '_'), # ...when 2 is F: dc2_dc1
                                                  ifelse(dyads$sex_1 == 'M', # ...when neither are F...
                                                         paste(dyads$dem_class_1, dyads$dem_class_2, sep = '_'), # ...when 1 is M: dc1_dc2
                                                         paste(dyads$dem_class_2, dyads$dem_class_1, sep = '_')))))) # ...when 2 is M or both are U: dc2_dc1
sort(unique(dyads$dem_type))

## add column for age difference between dyad
dyads$age_diff <- abs(as.numeric(dyads$age_cat_id_1) - as.numeric(dyads$age_cat_id_2))

## add column for total number of sightings per pair
dyads$count_dyad <- (dyads$count_1 + dyads$count_2) - dyads$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

## add column for total number of sightings per pair where they were NOT together
dyads$apart <- dyads$count_dyad - dyads$event_count

## reassign dyad numbers to remove gaps
dyads$node_1_nogaps <- as.integer(as.factor(dyads$node_1))
dyads$node_2_nogaps <- as.integer(as.factor(dyads$node_2))+1

## remove dead individuals (also removes all unknown sex calves)
dyads <- dyads %>% 
  filter(name_1 != 'Richard' & name_2 != 'Richard') %>% 
  filter(name_1 != 'Gabriel' & name_2 != 'Gabriel') %>% 
  filter(name_1 != 'Tori'    & name_2 != 'Tori')

######## Filter down to only males over 10 years old          ####
dyads <- dyads %>% 
  filter(sex_1 == 'M') %>% 
  filter(sex_2 == 'M') %>% 
  filter(as.numeric(age_cat_id_1) >= 3) %>% 
  filter(as.numeric(age_cat_id_2) >= 3)

### create nodes id variable which is 1:max(id) as currently covers all elephants, not just males
dyads$node_1_males <- as.integer(as.factor(dyads$node_1_nogaps))
dyads$node_2_males <- as.integer(as.factor(dyads$node_2_nogaps))+1

### standardise dyad_id for males only
dyads$dyad_males <- as.integer(as.factor(dyads$dyad_id))

### write out data as ready to go into model
write_csv(dyads, '../data_processed/step1_dataprocessing/motnp_binomialpairwiseevents.csv')

## clean up environment ready for next step
rm(list = ls()) ; gc()

### add time marker
print(paste0('Data wrangling completed at ', Sys.time(), '. Ready for analysis.'))

#---------------- 3) SRI      ##########
######## Import data          ####
counts_df <- read_csv('../data_processed/step1_dataprocessing/motnp_binomialpairwiseevents.csv') %>% 
  dplyr::select(-sex_1, -sex_2)

######## Calculate SRI values ####
# calculate sri
counts_df$sri <- counts_df$event_count / counts_df$count_dyad

# calculate percentages of 0s and 1s using SRI
( length(which(counts_df$sri==0))/length(counts_df$sri) ) * 100
( length(which(counts_df$sri==1))/length(counts_df$sri) ) * 100
length(unique(c(counts_df$id_1, counts_df$id_2)))

# calculate percentages of 0s and 1s using SRI, assuming a 5 sighting threshold per elephant
subset5 <- counts_df %>% filter(count_1 >= 5 & count_2 >= 5)
( length(which(subset5$sri==0))/length(subset5$sri) ) * 100
( length(which(subset5$sri==1))/length(subset5$sri) ) * 100
length(unique(c(subset5$id_1, subset5$id_2)))

# calculate percentages of 0s and 1s using SRI, assuming a 10 sighting threshold per elephant
subset10 <- counts_df %>% filter(count_1 >= 10 & count_2 >= 10)
( length(which(subset10$sri==0))/length(subset10$sri) ) * 100
length(unique(c(subset10$id_1, subset10$id_2)))
100-(( length(unique(c(subset10$id_1, subset10$id_2))) / length(unique(c(counts_df$id_1, counts_df$id_2))) )*100)

## clean up
rm(list = ls()[!ls() %in% c('counts_df')]) ; gc()

#---------------- 4) BISoN with Uniform prior ##########
######## Prep data for models                 ####
## create nodes data frame
nodes <- data.frame(id = sort(unique(c(dyads$id_1,dyads$id_2))),  # all unique individuals
                    node = NA, age = NA, sightings = NA)          # data needed on each
for(i in 1:nrow(nodes)){
  ## extract data about individual from dyads data frame
  if(nodes$id[i] %in% dyads$id_1) {
    x <- dyads[dyads$id_1 == nodes$id[i], c('id_1','node_1','count_1','age_category_1')] %>% 
      distinct()
  } else { 
    x <- dyads[dyads$id_2 == nodes$id[i], c('id_2','node_2','count_2','age_category_2')] %>%
      distinct() 
  }
  colnames(x) <- c('id','node','count','age')
  ## add individual data
  nodes$node[i] <- x$node
  nodes$age[i] <- x$age
  nodes$sightings[i] <- x$count
}
rm(x,i) ; gc()

## create data list
n_chains <- 4
n_samples <- 1000
n_dyads <- nrow(dyads)
counts_ls <- list(
  n_dyads    = n_dyads,                  # total number of times one or other of the dyad was observed
  dyad_ids   = dyads$dyad_id,            # identifier for each dyad
  together   = dyads$event_count,        # count number of sightings seen together
  count_dyad = dyads$count_dyad          # count total number of times seen
)

######## Compile edge model                   ####
edge_binary_uniform <- cmdstan_model('other/methods_paper/edge_binary_uniform.stan')

######## Run model                            ####
fit_edges_uniform <- edge_binary_uniform$sample(
  data = counts_ls, 
  chains = n_chains, parallel_chains = n_chains,
  iter_warmup = n_samples, iter_sampling = n_samples)

######## Extract edges                        ####
edges <- fit_edges_uniform$draws() %>% as.data.frame()
edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
edge_names <- data.frame(name = colnames(edges1)) %>%
  separate(name, into = c('chain','weight'), sep = '.edge_')
colnames(edges1) <- edge_names$weight
edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight   # select only chain 2
edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight   # select only chain 3
edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight   # select only chain 4
edges <- rbind(edges1, edges2, edges3, edges4)
n_samples <- nrow(edges)
edge_samples1 <- edges1 %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
rm(edges1, edges2, edges3, edges4) ; gc()

## save model for plotting
save.image('../outputs/sparse_network_methods_figures/model_run_uniform.RData')
rm(list = ls()[! ls() %in% c('counts_ls','nodes','counts_df')]) ; gc()

#---------------- 5) BISoN with Default prior ##########
######## Compile edge model                   ####
edge_binary_default <- cmdstan_model('other/methods_paper/edge_binary_gaussian.stan')

######## Run model                            ####
fit_edges_default <- edge_binary_default$sample(
  data = counts_ls, 
  chains = n_chains, parallel_chains = n_chains,
  iter_warmup = n_samples, iter_sampling = n_samples)

######## Extract edges                        ####
edges <- fit_edges_default$draws() %>% as.data.frame()
edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
edge_names <- data.frame(name = colnames(edges1)) %>%
  separate(name, into = c('chain','weight'), sep = '.edge_')
colnames(edges1) <- edge_names$weight
edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight   # select only chain 3
edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight   # select only chain 4
edges <- rbind(edges1, edges2, edges3, edges4)
n_samples <- nrow(edges)
edge_samples1 <- edges1 %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
rm(edges1, edges2, edges3, edges4) ; gc()

## save model for plotting
save.image('../outputs/sparse_network_methods_figures/model_run_default.RData')
rm(list = ls()[! ls() %in% c('counts_ls','nodes','counts_df')]) ; gc()

#---------------- 6) BISoN with Skewed prior  ##########
######## Compile edge model                   ####
edge_binary_skewed <- cmdstan_model('other/methods_paper/edge_binary_skewed.stan')

######## Run model                            ####
fit_edges_skewed <- edge_binary_skewed$sample(
  data = counts_ls, 
  chains = n_chains, parallel_chains = n_chains,
  iter_warmup = n_samples, iter_sampling = n_samples)

######## Extract edges                        ####
edges <- fit_edges_skewed$draws() %>% as.data.frame()
edges <- edges[,(n_chains+1):ncol(edges)]        # remove lp__ columns
edges1 <- edges[,seq(1,ncol(edges)-3, by = 4)]   # select only chain 1
edge_names <- data.frame(name = colnames(edges1)) %>%
  separate(name, into = c('chain','weight'), sep = '.edge_')
colnames(edges1) <- edge_names$weight
edges2 <- edges[,seq(2,ncol(edges)-2, by = 4)] ; colnames(edges2) <- edge_names$weight  # select only chain 2
edges3 <- edges[,seq(3,ncol(edges)-1, by = 4)] ; colnames(edges3) <- edge_names$weight   # select only chain 3
edges4 <- edges[,seq(4,ncol(edges)-0, by = 4)] ; colnames(edges4) <- edge_names$weight   # select only chain 4
edges <- rbind(edges1, edges2, edges3, edges4)
n_samples <- nrow(edges)
edge_samples1 <- edges1 %>%
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'weight')
rm(edges1, edges2, edges3, edges4) ; gc()

## save model for plotting
save.image('../outputs/sparse_network_methods_figures/model_run_skewedß.RData')
rm(list = ls()[! ls() %in% c('counts_ls','nodes','counts_df')]) ; gc()

#---------------- 7) BISoN with Conditional prior ##########
######## Compile edge model                       ####
edge_binary <- cmdstan_model("models/edge_binary_conditional.stan")

######## Run model                                ####
fit_edges_motnp <- edge_binary$sample(
  data = counts_ls, 
  chains = n_chains, parallel_chains = n_chains,
  iter_warmup = n_samples, iter_sampling = n_samples)

######## Extract edges -- different code to others -- check that this definitely does the same thing! ####
posterior_samples <- fit_edges_motnp$draws()

### break down into parameters
edge_weights_matrix <- posterior_samples[,,2:(nrow(dyads)+1)]
rm(posterior_samples) ; gc()

### save edge samples -- convert from being an array where each layer is a set of 4 chains, to a data frame where all chains are saved together in long format
edges <- as.data.frame(edge_weights_matrix[,,1])                                          # extract matrix for first dyad
colnames(edges) <- c('chain1','chain2','chain3','chain4')                                 # rename columns so know which chain is which
edges <- pivot_longer(edges, everything(), values_to = 'edge_draw', names_to = 'chain')   # convert to long format
edges$dyad <- counts_ls$dyad_ids[1]                                                       # add dyad ID for first dyad
edges$position <- rep(1:n_samples, each = n_chains)                                       # add chain position for each dyad
for(i in 2:n_dyads){                                                                      # repeat for all other dyads and append
  x <- as.data.frame(edge_weights_matrix[,,i])
  colnames(x) <- c('chain1','chain2','chain3','chain4')
  x <- pivot_longer(x, everything(), values_to = 'edge_draw', names_to = 'chain')
  x$dyad <- counts_ls$dyad_ids[i]
  x$position <- rep(1:n_samples, each = n_chains)
  edges <- rbind(edges, x)
}
save.image('../outputs/sparse_network_methods_figures/model_run_conditional.RData')
rm(list = ls()) ; gc()

#---------------- 8) Create SRI plots for paper ##########
pdf('../outputs/sparse_network_methods_figures/all_plots.pdf')
colours <- c('#21918c','#440154')

######## SRI Edge weights: 2 panel plot of SRI edges and edge_vs_sightings ####
#### Edge bar graph                                                        ####
# bar plot
(edges_sri <- ggplot()+
   geom_bar(data = counts_df, aes(x = round(sri, 3)), fill = rgb(33/255, 145/255, 140/255),
            colour = rgb(33/255, 145/255, 140/255), linewidth = 0.2)+
   #geom_bar(data = counts_df[counts_df$count_dyad > 10,], aes(x = round(sri, 3)), fill = '#fde725', colour = 'transparent')+
   scale_x_continuous(name = 'SRI value')+
   scale_y_continuous(name = 'number of dyads',
                      #limits = c(0,1000),
                      expand = c(0,0))+
   annotate('text', x = 0.44, y = 165,
            label = paste0('+ ', length(which(counts_df$sri == 0)), ' dyads with \nSRI = 0'),
            size = unit(4, 'pt'),
            colour = rgb(33/255, 145/255, 140/255))+
   coord_cartesian(ylim = c(0,200))
)
ggsave(filename = 'edges_sri.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = edges_sri, device = 'png', width = 700, height = 700, units = 'px')

# redo annotation on so that size isn't too large for multi plot
(edges_sri <- ggplot()+
    geom_bar(data = counts_df, aes(x = round(sri, 3)), fill = rgb(33/255, 145/255, 140/255),
             colour = rgb(33/255, 145/255, 140/255), linewidth = 0.2)+
    scale_x_continuous(name = 'SRI value')+
    scale_y_continuous(name = 'number of dyads',
                       expand = c(0,0))+
    annotate('text', x = 0.5, y = 165,
             label = paste0(length(which(counts_df$sri == 0)), ' dyads with \nSRI = 0'),
             size = unit(4, 'pt'),
             colour = rgb(33/255, 145/255, 140/255))+
    coord_cartesian(ylim = c(0,200))
)

#### Edges vs Sightings                                                    ####
## basic plot
(edgesightings_sri.1 <- ggplot(counts_df, aes(x = count_dyad, y = sri))+
   geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
              size = 0.5,
              shape = 19)+
   geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
   scale_x_continuous(name = 'total dyad sightings')+
   scale_y_continuous(name = 'SRI weight', limits = c(-0.02,1.02), expand = c(0,0))
)

## split lines by if elephants were ever seen in the same group
counts_df$together <- ifelse(counts_df$event_count == 0,
                             'never together', 'together at least once')
(edgesightings_sri.2 <- ggplot()+
    geom_point(data = counts_df,
               aes(x = count_dyad, y = sri),
               colour = rgb(94/255, 201/255, 98/255, 0.2),
               size = 1, shape = 19)+
    geom_smooth(data = counts_df[counts_df$together == 'together at least once',],
                aes(x = count_dyad, y = sri, colour = together))+
    geom_smooth(data = data.frame(x = 2:max(counts_df$count_dyad[counts_df$event_count == 0]),
                                  y = 0, together = 'never together'),
                aes(x = x, y = y, colour = together))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'SRI weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_colour_manual(values = colours,
                        breaks = c('never together','together at least once'),
                        name = NULL)+
    theme(legend.position = 'right',
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)))

#### Merge                                                                 ####
(combined <- (edges_sri + edgesightings_sri.2)+
   plot_annotation(tag_levels = 'a'))
ggsave(filename = 'outputs_sri.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 700, units = 'px')

######## SRI Eigenvector Centrality: 2 panel plot of SRI eigenvector against sighting count and number of non-associations ####
#### Individual sighting count                                                                                             ####
# calculate SRI and identify individuals
counts_df$sri <- counts_df$event_count / counts_df$count_dyad
elephants <- unique(c(counts_df$id_1, counts_df$id_2))

# calculate adjacency matrix
sri_adjmat <- matrix(NA, nrow = length(elephants), ncol = length(elephants), dimnames = list(x = elephants, y = elephants))
for(i in 1:nrow(sri_adjmat)){
  for(j in 1:ncol(sri_adjmat)){
    if(i <= j){
      sri_adjmat[i,j] <- ifelse(i == j, 0,
                                counts_df$sri[which(counts_df$id_1 == rownames(sri_adjmat)[i] &
                                                      counts_df$id_2 == colnames(sri_adjmat)[j])])
    } else {
      sri_adjmat[i,j] <- sri_adjmat[j,i]
    }
  }
}

# create network
sri_net <- igraph::graph_from_adjacency_matrix(sri_adjmat, weighted = T, mode = 'undirected', diag = FALSE)

# calculate eigenvector centrality
eigen_sri <- igraph::eigen_centrality(sri_net)$vector %>%
  as.data.frame()
colnames(eigen_sri) <- 'eigen'
eigen_sri$id <- rownames(eigen_sri)

# create data frame showing number of sightings per individual with other elephants
together0 <- counts_df %>% filter(event_count == 0)
eigen_sri$together0 <- NA ; for(i in 1:nrow(eigen_sri)) {
  x <- together0 %>%
    filter(id_1 == eigen_sri$id[i] | id_2 == eigen_sri$id[i])
  eigen_sri$together0[i] <- nrow(x)
}

# calculate number of sightings per elephant
counts <- data.frame(id = elephants, count = NA)
for(i in 1:nrow(counts)){
  counts$count[i] <- ifelse(i == nrow(counts),
                            unique(counts_df$count_2[which(counts_df$id_2 == counts$id[i])]),
                            unique(counts_df$count_1[which(counts_df$id_1 == counts$id[i])]))
}
eigen_sri <- eigen_sri %>%
  left_join(counts, by = 'id')

# plot
(eigensightings_sri <- ggplot(eigen_sri, aes(x = count, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255),
               size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

#### Total number of dyads where together = 0                                                                              ####
## no trend line
(eigen0_sri.1 <- ggplot(eigen_sri, aes(x = together0, y = eigen))+
   geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
   scale_x_continuous(name = 'dyads together = 0')+
   scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## with trend line
(eigen0_sri.2 <- ggplot(eigen_sri, aes(x = together0, y = eigen))+
    geom_point(colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

#### Merge                                                                                                                 ####
(eigen0_sri.2 + eigensightings_sri)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_sri.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 1600, height = 700, units = 'px')
rm(eigen0_sri.2,eigensightings_sri.2) ; gc()

save.image('../outputs/sparse_network_methods_figures/sri.RData')

#---------------- 9) Create Unconditional BISoN plots for paper ####
######## Unconditional Priors: 3 panel plot of uniform/default BISoN/right-skewed BISoN for priors ####
x <- seq(0, 1, length = 100)

#### Uniform, uniform(0,1)                                                                         ####
# calculate probability of all values of x assuming completely flat prior
flat_prior <- dunif(x = x, min = 0, max = 1)
data <- data.frame(x = x, density = flat_prior)

# plot
(priors_sri <- ggplot(data)+
    geom_line(aes(x = x, y = density),
              linewidth = 1.2,
              colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density', limits = c(0,1.05), expand = c(0,0))
)

# remove all information around it
inset_sri <- ggplot(data)+
  geom_line(aes(x = x, y = density),
            linewidth = 1.2,
            colour = rgb(33/255, 145/255, 140/255))+
  scale_x_continuous(#breaks = c(0,0.5,1),
    name = NULL)+
  scale_y_continuous(name = NULL,
                     #breaks = c(0,0.5,1),
                     limits = c(0,1.05),
                     expand = c(0,0))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())

#### Default, normal(0,2.5)                                                                        ####
# calculate probability of all values of x assuming default bisonR prior shape
default_prior <- get_default_priors('binary') # default edge prior = normal(0, 2.5)
default_prior <- c(0, 2.5)
default_prior <- dnorm(x = logit(x), mean = default_prior[1], sd = default_prior[2])
data <- data.frame(x = x, density = default_prior)

# plot
(priors_default <- ggplot(data)+
    geom_line(aes(x = x, y = density),
              linewidth = 1.2,
              colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density',
                       limits = c(0,0.2),
                       #limits = c(0,1),
                       expand = c(0,0))
)

inset_default <- ggplot(data)+
  geom_line(aes(x = x, y = density),
            linewidth = 1.2,
            colour = rgb(33/255, 145/255, 140/255))+
  scale_x_continuous(#breaks = c(0,0.5,1),
    name = NULL)+
  scale_y_continuous(name = NULL,
                     #breaks = c(0,0.1,0.2),
                     limits = c(0,0.2),
                     expand = c(0,0))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())

#### Skewed, beta(0.7,5)                                                                           ####
# calculate probability of all values of x assuming right skewed prior shape
skewed_prior <- c(1, 5)
skewed_prior <- dbeta(x = x, shape1 = skewed_prior[1], shape2 = skewed_prior[2])
data <- data.frame(x = x, density = skewed_prior)

# plot
(priors_skewed <- ggplot(data)+
    geom_line(aes(x = x, y = density),
              linewidth = 1.2,
              colour = rgb(33/255, 145/255, 140/255))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density',
                       #limits = c(0,0.3),
                       #limits = c(0,1),
                       expand = c(0,0))
)

inset_skewed <- ggplot(data)+
  geom_line(aes(x = x, y = density),
            linewidth = 1.2,
            colour = rgb(33/255, 145/255, 140/255))+
  scale_x_continuous(breaks = c(0,0.5,1),
                     name = NULL)+
  scale_y_continuous(name = NULL,
                     #breaks = c(0,1,2,3,4,5),
                     limits = c(-0.1,5),
                     expand = c(0,0))+
  theme(axis.text = element_blank(), axis.ticks = element_blank())

#### Merge                                                                                         ####
(priors_sri + priors_default + priors_skewed) +
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'priors_unconditional.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2100, height = 700, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('inset_sri','inset_default','inset_skewed','colours')]) ; gc()

######## Unconditional Outputs: 6 panel plot of uniform/default/skewed for posterior and edge_vs_sightings ####
#### Uniform, uniform(0,1)                                                                                 ####
load('../outputs/sparse_network_methods_figures/model_run_uniform.RData')

## plot posterior distribution
(figure_uniform_posterior <- ggplot(data = edge_samples1)+
    geom_density(aes(group = parameter, x = weight),
                 colour = rgb(33/255, 145/255, 140/255, 0.1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density')
)
ggsave(filename = 'posterior_uniform_alldyads.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure_uniform_posterior, device = 'png', width = 1600, height = 700, units = 'px')

## split colours by if elephants were ever seen grouping together
edge_samples1$dyad_males <- rep(1:nrow(counts_df), n_samples/4)
edge_samples1 <- edge_samples1 %>%
  left_join(counts_df[,c('dyad_males','id_1','id_2','count_1','count_2','count_dyad','event_count','sri')],
            by = 'dyad_males') %>%
  mutate(together0 = ifelse(event_count == 0, 'never together', 'together at least once'))

## plot a sample of 100 randomly selected dyads
set.seed(15) ; plot_dyads <- sample(1:nrow(counts_df), size = 100, replace = F)
plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]
(figure_uniform_posterior <- ggplot(data = edges_subset)+
    geom_density(aes(group = parameter, x = weight, colour = together0),
                 show.legend = F, linewidth = 0.2)+
    stat_density(aes(group = parameter, x = weight, colour = together0),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(0,1))+
    scale_y_continuous(name = 'density', limits = c(0,40))+
    scale_colour_manual(values = colours,
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
)
ggsave(filename = 'posterior_uniform_sampledyads.png',
       path = '../outputs/sparse_network_methods_figures/in',
       plot = figure_uniform_posterior, device = 'png',
       width = 1000, height = 1000, units = 'px')

## create dataframe to plot mean edge weight vs dyad sighting count
counts_df$together <- ifelse(counts_df$event_count == 0,
                             'never together',
                             'together at least once')

averages <- data.frame(dyad = colnames(edges),
                       mean = apply(edges, 2, mean),
                       median = apply(edges, 2, median))
averages$dyad_males <- 1:nrow(averages)
averages <- averages %>%
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','together','count_dyad','count_1','count_2','dyad_males')],
            by = 'dyad_males')

## all median points, no uncertainty
(figure_uniform_edgesightings.1 <- ggplot(averages, aes(x = count_dyad, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'median weight', limits = c(-0.02,1.02), expand = c(0,0))
)

## all median plots, no uncertainty, split lines by if they were ever seen together
(figure_uniform_edgesightings.2 <- ggplot()+
    geom_point(data = averages,
               aes(x = count_dyad, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = averages,
                aes(x = count_dyad, y = median,
                    linetype = as.factor(together)),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'median weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_linetype_manual(name = NULL,
                          values = c(1,6))+
    theme(legend.position = 'bottom',
          legend.key.height = unit(4, 'mm'),
          legend.text = element_text(size = 8)))

## all median points, showing uncertainty, split lines by if they were ever seen together
(figure_uniform_edgesightings.3 <- ggplot()+
    geom_point(data = edge_samples1,
               aes(x = count_dyad, y = weight),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages,
               aes(x = count_dyad, y = median),
               colour = rgb(94/255, 201/255, 98/255, 0.2),
               size = 1, shape = 19)+
    geom_smooth(data = averages,
                aes(x = count_dyad, y = median,
                    colour = as.factor(together)))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'edge weight',
                       limits = c(-0.02,1.02),
                       expand = c(0,0))+
    scale_colour_manual(values = colours)+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))
)

## save outputs
save.image('../outputs/sparse_network_methods_figures/plots_bisonuniform.RData')
rm(list = ls()[!ls() %in% c('counts_df','n_chains',
                            'inset_sri','inset_default','inset_skewed','inset_conditional',
                            'colours')]) ; gc()

#### Default, normal(0,2.5)                                                                                ####
load('../outputs/sparse_network_methods_figures/model_run_default.RData')

# plot posterior distribution
(figure_default_posterior <- ggplot(data = edge_samples1)+
    geom_density(aes(group = parameter, x = weight),
                 colour = rgb(33/255, 145/255, 140/255, 0.1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density')
)
ggsave(filename = 'posterior_default_alldyads.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure_default_posterior, device = 'png',
       width = 1600, height = 700, units = 'px')

## split colours by if elephants were ever seen grouping together
edge_samples1$dyad_males <- rep(1:nrow(counts_df), n_samples/4)
edge_samples1 <- edge_samples1 %>%
  left_join(counts_df[,c('dyad_males','id_1','id_2','count_1','count_2','count_dyad','event_count','sri')],
            by = 'dyad_males') %>%
  mutate(together0 = ifelse(event_count == 0, 'never together', 'together at least once'))

## plot a sample of 100 randomly selected dyads
set.seed(15) ; plot_dyads <- sample(1:nrow(counts_df), size = 100, replace = F)
plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]
(figure_default_posterior <- ggplot(data = edges_subset)+
    geom_density(aes(group = parameter, x = weight, colour = together0),
                 show.legend = F, linewidth = 0.2)+
    stat_density(aes(group = parameter, x = weight, colour = together0),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(0,1))+
    scale_y_continuous(name = 'density', limits = c(0,40))+
    scale_colour_manual(values = colours,
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
)

## create dataframe to plot mean edge weight vs dyad sighting count
counts_df$together <- ifelse(counts_df$event_count == 0,
                             'never together',
                             'together at least once')

averages <- data.frame(dyad = colnames(edges),
                       mean = apply(edges, 2, mean),
                       median = apply(edges, 2, median))
averages$dyad_males <- 1:nrow(averages)

averages <- averages %>%
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','together','count_dyad','count_1','count_2','dyad_males')],
            by = 'dyad_males')

## all median points, no uncertainty
(figure_default_edgesightings.1 <- ggplot(averages, aes(x = count_dyad, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)

## all median plots, no uncertainty, split lines by if they were ever seen together
figure_default_edgesightings.2 <- ggplot()+
    geom_point(data = averages,
               aes(x = count_dyad, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = averages,
                aes(x = count_dyad, y = median, linetype = as.factor(together)),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))+
    theme(legend.position = 'none')
(figure_default_edgesightings.2 +
  scale_linetype_manual(name = NULL,
                        values = c(1,6))+
  theme(legend.position = c(0.655,0.825),
        legend.background = element_rect(fill = 'white', colour = 'black'),
        legend.key.height = unit(4, 'mm'),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)))

## all median points, showing uncertainty, split lines by if they were ever seen together
(figure_default_edgesightings.3 <- ggplot()+
    geom_point(data = edge_samples1,
               aes(x = count_dyad, y = weight),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5, shape = 19)+
    geom_point(data = averages,
               aes(x = count_dyad, y = median),
               colour = rgb(94/255, 201/255, 98/255, 0.2),
               size = 1, shape = 19)+
    geom_smooth(data = averages,
                aes(x = count_dyad, y = median,
                    colour = as.factor(together)))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'edge weight',
                       limits = c(-0.02,1.02),
                       expand = c(0,0))+
    scale_colour_manual(values = colours)+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))
)

## save outputs
save.image('../outputs/sparse_network_methods_figures/plots_bisondefault.RData')
rm(list = ls()[!ls() %in% c('counts_df','n_chains',
                            'inset_sri','inset_default','inset_skewed','inset_conditional',
                            'colours')]) ; gc()

#### Skewed, beta(0.7,5)                                                                                   ####
load('../outputs/sparse_network_methods_figures/model_run_skewed.RData')

# plot posterior distribution
(figure_skewed_posterior <- ggplot(data = edge_samples1)+
    geom_density(aes(group = parameter, x = weight),
                 colour = rgb(33/255, 145/255, 140/255, 0.1))+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density')
)
ggsave(filename = 'posterior_skewed_alldyads.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure_skewed_posterior, device = 'png', width = 1600, height = 700, units = 'px')

## split colours by if elephants were ever seen grouping together
edge_samples1$dyad_males <- rep(1:nrow(counts_df), n_samples/4)
edge_samples1 <- edge_samples1 %>%
  left_join(counts_df[,c('dyad_males','id_1','id_2','count_1','count_2','count_dyad','event_count','sri')],
            by = 'dyad_males') %>%
  mutate(together0 = ifelse(event_count == 0, 'never together', 'together at least once'))

## plot a sample of 100 randomly selected dyads
set.seed(15) ; plot_dyads <- sample(1:nrow(counts_df), size = 100, replace = F)
plot_dyad_ids <- unique(edge_samples1$parameter)[plot_dyads]
edges_subset <- edge_samples1[edge_samples1$parameter %in% plot_dyad_ids,]
(figure_skewed_posterior <- ggplot(data = edges_subset)+
    geom_density(aes(group = parameter, x = weight, colour = together0),
                 show.legend = F, linewidth = 0.2)+
    stat_density(aes(group = parameter, x = weight, colour = together0),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(0,1))+
    scale_y_continuous(name = 'density', limits = c(0,40))+
    scale_colour_manual(values = colours,
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
)

# clean up and save
save.image('../outputs/sparse_network_methods_figures/model_run_skewed.RData')

## create dataframe to plot mean edge weight vs dyad sighting count
counts_df$together <- ifelse(counts_df$event_count == 0,
                             'never together',
                             'together at least once')

averages <- data.frame(dyad = colnames(edges),
                       mean = apply(edges, 2, mean),
                       median = apply(edges, 2, median))
averages$dyad_males <- 1:nrow(averages)
averages <- averages %>%
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','together','count_dyad','count_1','count_2','dyad_males')],
            by = 'dyad_males')

## all median points, no uncertainty
(figure_skewed_edgesightings.1 <- ggplot(averages,
                                         aes(x = count_dyad, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'median weight',
                       limits = c(-0.02,1.02), expand = c(0,0))
)

## all median plots, no uncertainty, split lines by if they were ever seen together
(figure_skewed_edgesightings.2 <- ggplot()+
    geom_point(data = averages,
               aes(x = count_dyad, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = averages,
                aes(x = count_dyad, y = median, linetype = as.factor(together)),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'median weight',
                       limits = c(-0.02,1.02), expand = c(0,0))+
    scale_linetype_manual(name = NULL,
                          values = c(1,6))+
    theme(legend.position = c(0.655,0.825),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
)

## all median points, showing uncertainty, split lines by if they were ever seen together
(figure_skewed_edgesightings.3 <- ggplot()+
    geom_point(data = edge_samples1,
               aes(x = count_dyad, y = weight),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages,
               aes(x = count_dyad, y = median),
               colour = rgb(94/255, 201/255, 98/255, 0.2),
               size = 1,
               shape = 19)+
    geom_smooth(data = averages,
                aes(x = count_dyad, y = median,
                    colour = as.factor(together))
    )+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'edge weight',
                       limits = c(-0.02,1.02),
                       expand = c(0,0))+
    scale_colour_manual(values = colours)+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))
)

## save outputs
save.image('../outputs/sparse_network_methods_figures/plots_bisonskewed.RData')

#### Merge                                                                                                 ####
#rm(list = ls()) ; gc()
# load('../outputs/sparse_network_methods_figures/plots_bisonskewed.RData')
rm(list = ls()[!ls() %in% c('figure_skewed_posterior','figure_skewed_edgesightings.3',
                            'inset_sri','inset_default','inset_skewed','inset_conditional',
                            'colours')]) ; gc()
load('../outputs/sparse_network_methods_figures/plots_bisondefault.RData')
rm(list = ls()[!ls() %in% c('figure_default_posterior','figure_default_edgesightings.3',
                            'figure_skewed_posterior','figure_skewed_edgesightings.3',
                            'inset_sri','inset_default','inset_skewed','inset_conditional',
                            'colours')]) ; gc()
load('../outputs/sparse_network_methods_figures/plots_bisonuniform.RData')
rm(list = ls()[!ls() %in% c('figure_uniform_posterior','figure_uniform_edgesightings.3',
                            'figure_default_posterior','figure_default_edgesightings.3',
                            'figure_skewed_posterior','figure_skewed_edgesightings.3',
                            'inset_sri','inset_default','inset_skewed','inset_conditional',
                            'colours')]) ; gc()
save.image('../outputs/sparse_network_methods_figures/plots_unconditional.RData')

(combined_top <- (figure_uniform_posterior + 
                    figure_default_posterior + 
                    figure_skewed_posterior) +
    plot_annotation(tag_levels = 'a') +
    plot_layout(guides = 'collect') &
    theme(legend.position = 'bottom'))
ggsave(filename = 'posterior_unconditional_noinset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2700, height = 700, units = 'px')

combined_bottom <- (figure_uniform_edgesightings.3 + figure_default_edgesightings.3 + figure_skewed_edgesightings.3)+
  plot_annotation(tag_levels = list(c('d','e','f')))
(combined_bottom <- combined_bottom + 
    plot_layout(guides = 'collect') & theme(legend.position = 'bottom'))
ggsave(filename = 'edgesightings_unconditional_noinset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2700, height = 700, units = 'px')

(combined_top / combined_bottom)
ggsave(filename = 'outputs_unconditional_noinset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2700, height = 1600, units = 'px')

uniform_posterior_inset <- wrap_elements(figure_uniform_posterior + 
                                           inset_element(inset_sri,
                                                         left = 0.5, right = 0.99,
                                                         bottom = 0.5, top = 0.99) +
                                           theme(legend.position = 'none'))
default_posterior_inset <- wrap_elements(figure_default_posterior + 
                                           inset_element(inset_default,
                                                         left = 0.5, right = 0.99,
                                                         bottom = 0.5, top = 0.99) +
                                           theme(legend.position = 'bottom', legend.title = element_blank()))
skewed_posterior_inset <- wrap_elements(figure_skewed_posterior + 
                                          inset_element(inset_skewed,
                                                        left = 0.5, right = 0.99, 
                                                        bottom = 0.5, top = 0.99) +
                                          theme(legend.position = 'none'))

combined_top <- (uniform_posterior_inset + default_posterior_inset + skewed_posterior_inset)+
  plot_annotation(tag_levels = 'a')
(combined_top <- combined_top + 
    theme(legend.position = 'bottom', legend.title = element_blank()) +
    plot_layout(guides = 'collect')
)
ggsave(filename = 'posterior_unconditional_inset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2700, height = 700, units = 'px')

(combined_top / combined_bottom)
ggsave(filename = 'outputs_unconditional_inset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2700, height = 1600, units = 'px')

rm(list = ls()[!ls() %in% c('inset_sri','inset_default','inset_skewed',
                            'colours')]) ; gc()

######## Unconditional Eigenvector Centrality: 6 panel plot of eigenvector against sighting count and number of non-associations ####
#### Uniform, uniform(0,1)                                                                                                       ####
load('../outputs/sparse_network_methods_figures/model_run_uniform.RData')

## create adjacency array from uniform edge posterior
uniform_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                        dimnames = list(elephants, elephants, NULL))
N <- nrow(counts_df)
for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  uniform_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edges[1:1000, i]
}

## make matrix to save eigenvector values to
eigen_uniform <- matrix(NA, nrow = 1000, ncol = length(elephants),
                        dimnames = list(1:1000, elephants))

## calculate eigenvector scores
for(i in 1:nrow(eigen_uniform)){
  uniform_net <- igraph::graph_from_adjacency_matrix(uniform_adjarr[,,i], weighted = T, mode = 'undirected', diag = FALSE)
  eigen <- as.data.frame(igraph::eigen_centrality(uniform_net)$vector)
  eigen_uniform[i,] <- eigen$`igraph::eigen_centrality(uniform_net)$vector`
}

## make into long format data frame
eigen_uniform <- eigen_uniform %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>%
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

## create data frame of average eigenvector centrality scores
averages_uniform <- eigen_uniform %>%
  select(id, mean, median, together0, count) %>%
  distinct()

## plot eigenvector centrality vs individual sighting count, with uncertainty, no trend line
(eigensightings_uniform.1 <- ggplot()+
    geom_point(data = eigen_uniform, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_uniform, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs individual sighting count, with uncertainty and trend line
(eigensightings_uniform.2 <- ggplot()+
    geom_point(data = eigen_uniform, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5, 
               shape = 19)+
    geom_point(data = averages_uniform, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = eigen_uniform, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs individual sighting count, no uncertainty, with trend line
(eigensightings_uniform.3 <- ggplot(data = averages_uniform, aes(x = count, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0,0))
)

## eigenvector centrality vs total number of dyads where together = 0, with uncertainty, no trend line
(eigen0_uniform.1 <- ggplot()+
   geom_point(data = eigen_uniform,
              aes(x = together0, y = eigen),
              colour = rgb(253/255, 231/255, 37/255, 0.01),
              size = 0.5,
              shape = 19)+
   geom_point(data = averages_uniform,
              aes(x = together0, y = median),
              colour = rgb(33/255, 145/255, 140/255),
              size = 0.5,
              shape = 19)+
   scale_x_continuous(name = 'dyads together = 0')+
   scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs total number of dyads where together = 0, with uncertainty and trend line
(eigen0_uniform.2 <- ggplot()+
    geom_point(data = eigen_uniform,
               aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_uniform,
               aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = eigen_uniform,
                aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs total number of dyads where together = 0, no uncertainty, with trend line
(eigen0_uniform.3 <- ggplot(averages_uniform)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5,
               shape = 19)+
    geom_smooth(aes(x = together0, y = median),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)

# save workspace
rm(list = ls()[!ls() %in% c('eigen_uniform','averages_uniform',
                            'eigensightings_uniform.1','eigensightings_uniform.2','eigensightings_uniform.3',
                            'eigen0_uniform.1','eigen0_uniform.2','eigen0_uniform.3',
                            'inset_sri','inset_default','inset_skewed','inset_conditional','colours')]) ; gc()
save.image('../outputs/sparse_network_methods_figures/eigen_checks_unconditional.RData')

#### Default, normal(0,2.5)                                                                                                      ####
## eigenvector centrality vs individual sighting count
load('../outputs/sparse_network_methods_figures/model_run_default.RData')

# get adjacency array for default prior
default_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                        dimnames = list(elephants, elephants, NULL))
for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  default_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edges[1:1000, i]
}

# make matrix to save eigenvector values to
eigen_default <- matrix(NA, nrow = 1000, ncol = length(elephants),
                        dimnames = list(1:1000, elephants))

# calculate eigenvector scores
for(i in 1:nrow(eigen_default)){
  default_net <- igraph::graph_from_adjacency_matrix(default_adjarr[,,i], weighted = T, mode = 'undirected', diag = FALSE)
  eigen <- as.data.frame(igraph::eigen_centrality(default_net)$vector)
  eigen_default[i,] <- eigen$`igraph::eigen_centrality(default_net)$vector`
}
rm(default_net, dyad_row, eigen); gc()

# make long format data frame
eigen_default <- eigen_default %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>%
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

# calculate averages
averages_default <- eigen_default %>%
  select(id, mean, median, together0, count) %>%
  distinct()

## plot eigenvector centrality vs individual sighting count, with uncertainty, no trend line
(eigensightings_default.1 <- ggplot()+
    geom_point(data = eigen_default, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_default, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs individual sighting count, with uncertainty and trend line
(eigensightings_default.2 <- ggplot()+
    geom_point(data = eigen_default, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_default, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    geom_smooth(data = eigen_default, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs individual sighting count, no uncertainty, with trend line
(eigensightings_default.3 <- ggplot(data = averages_default, aes(x = count, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255),
               size = 0.5,
               shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'median eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)

## eigenvector centrality vs total number of dyads where together = 0, with uncertainty, no trend line
(eigen0_default.1 <- ggplot()+
   geom_point(data = eigen_default, aes(x = together0, y = eigen),
              colour = rgb(253/255, 231/255, 37/255, 0.01),
              size = 0.5,
              shape = 19)+
   geom_point(data = averages_default, aes(x = together0, y = median, colour = count),
              colour = rgb(33/255, 145/255, 140/255),
              size = 0.5,
              shape = 19)+
   scale_x_continuous(name = 'dyads together = 0')+
   scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs total number of dyads where together = 0, with uncertainty and trend line
(eigen0_default.2 <- ggplot()+
    geom_point(data = eigen_default, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5,
               shape = 19)+
    geom_point(data = averages_default, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_default, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs total number of dyads where together = 0, no uncertainty, with trend line
(eigen0_default.3 <- ggplot(averages_default)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5,
               shape = 19)+
    geom_smooth(aes(x = together0, y = median),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)

# save workspace
rm(list = ls()[!ls() %in% c('eigen_uniform','averages_uniform',
                            'eigen_default','averages_default',
                            'eigensightings_uniform.1','eigensightings_uniform.2','eigensightings_uniform.3',
                            'eigensightings_default.1','eigensightings_default.2','eigensightings_default.3',
                            'eigen0_uniform.1','eigen0_uniform.2','eigen0_uniform.3',
                            'eigen0_default.1','eigen0_default.2','eigen0_default.3',
                            'inset_sri','inset_default','inset_skewed','inset_conditional','colours')]) ; gc()
save.image('../outputs/sparse_network_methods_figures/eigen_checks_unconditional.RData')

#### Skewed, beta(0.7,5)                                                                                                         ####
## eigenvector centrality vs individual sighting count
load('../outputs/sparse_network_methods_figures/model_run_skewed.RData')

# get adjacency array for skewed prior
skewed_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                       dimnames = list(elephants, elephants, NULL))
for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  skewed_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edges[1:1000, i]
}

# make matrix to save eigenvector values to
eigen_skewed <- matrix(NA, nrow = 1000, ncol = length(elephants),
                       dimnames = list(1:1000, elephants))

# calculate eigenvector scores
for(i in 1:nrow(eigen_skewed)){
  skewed_net <- igraph::graph_from_adjacency_matrix(skewed_adjarr[,,i], weighted = T, mode = 'undirected', diag = FALSE)
  eigen <- as.data.frame(igraph::eigen_centrality(skewed_net)$vector)
  eigen_skewed[i,] <- eigen$`igraph::eigen_centrality(skewed_net)$vector`
}

# make long format data frame
eigen_skewed <- eigen_skewed %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>%
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

# calculate averages
averages_skewed <- eigen_skewed %>%
  select(id, mean, median, together0, count) %>%
  distinct()

## plot eigenvector centrality vs individual sighting count, with uncertainty, no trend line
(eigensightings_skewed.1 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_skewed, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs individual sighting count, with uncertainty and trend line
(eigensightings_skewed.2 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_skewed, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_skewed, aes(x = count, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs individual sighting count, no uncertainty, with trend line
(eigensightings_skewed.3 <- ggplot(data = averages_skewed, aes(x = count, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean eigenvector centrality',
                       limits = c(0,1),
                       expand = c(0,0))
)

## eigenvector centrality vs total number of dyads where together = 0, with uncertainty, no trend line
(eigen0_skewed.1 <- ggplot()+
   geom_point(data = eigen_skewed, aes(x = together0, y = eigen),
              colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
   geom_point(data = averages_skewed, aes(x = together0, y = median),
              colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
   scale_x_continuous(name = 'dyads together = 0')+
   scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs total number of dyads where together = 0, with uncertainty and trend line
(eigen0_skewed.2 <- ggplot()+
    geom_point(data = eigen_skewed, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01), size = 0.5, shape = 19)+
    geom_point(data = averages_skewed, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255), size = 0.5, shape = 19)+
    geom_smooth(data = eigen_skewed, aes(x = together0, y = eigen),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

## eigenvector centrality vs total number of dyads where together = 0, no uncertainty, with trend line
(eigen0_skewed.3 <- ggplot(averages_skewed)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5, shape = 19)+
    geom_smooth(aes(x = together0, y = median),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')
)

# save workspace
rm(list = ls()[!ls() %in% c('eigen_uniform','averages_uniform',
                            'eigen_default','averages_default',
                            'eigen_skewed', 'averages_skewed',
                            'eigensightings_uniform.1','eigensightings_uniform.2','eigensightings_uniform.3',
                            'eigensightings_default.1','eigensightings_default.2','eigensightings_default.3',
                            'eigensightings_skewed.1', 'eigensightings_skewed.2', 'eigensightings_skewed.3',
                            'eigen0_uniform.1','eigen0_uniform.2','eigen0_uniform.3',
                            'eigen0_default.1','eigen0_default.2','eigen0_default.3',
                            'eigen0_skewed.1', 'eigen0_skewed.2', 'eigen0_skewed.3',
                            'inset_sri','inset_default','inset_skewed','inset_conditional','colours')]) ; gc()
save.image('../outputs/sparse_network_methods_figures/eigen_checks_unconditional.RData')

#### Merge                                                                                                                       ####
rm(list = ls()[!ls() %in% c('eigen0_uniform.2','eigensightings_uniform.2',
                            'eigen0_default.2','eigensightings_default.2',
                            'eigen0_skewed.2', 'eigensightings_skewed.2',
                            'inset_sri','inset_default','inset_skewed','colours')]) ; gc()

## plot without insets included
(eigen0_uniform.2 + eigen0_default.2 + eigen0_skewed.2) /
  (eigensightings_uniform.2 + eigensightings_default.2 + eigensightings_skewed.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_unconditional_noinset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2400, height = 1600, units = 'px')

## set up insets
uniform_eigen0_inset <- wrap_elements(eigen0_uniform.2 + 
                                        inset_element(inset_sri,
                                                      left = 0.01, right = 0.25,
                                                      bottom = 0.75, top = 0.99))
default_eigen0_inset <- wrap_elements(eigen0_default.2 + 
                                        inset_element(inset_default,
                                                      left = 0.01, right = 0.25,
                                                      bottom = 0.75, top = 0.99))
skewed_eigen0_inset <- wrap_elements(eigen0_skewed.2 + 
                                       inset_element(inset_skewed,
                                                     left = 0.01, right = 0.25,
                                                     bottom = 0.75, top = 0.99))

## plot with priors inset into top row
(uniform_eigen0_inset + default_eigen0_inset + skewed_eigen0_inset) /
  (eigensightings_uniform.2 + eigensightings_default.2 + eigensightings_skewed.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_unconditional_inset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 2400, height = 1600, units = 'px')

#---------------- 10) Create Conditional BISoN plots for paper ####
######## Conditional Prior: 1 panel plot of conditional prior  ####
# define sequence over which to plot
x <- seq(0, 1, length = 100)

# calculate probability of all values of x assuming completely flat prior
conditional1 <- dbeta(x = x, shape1 = 0.7, shape2 = 10)
conditional2 <- dbeta(x = x, shape1 = 1, shape2 = 5)
data <- data.frame(x = x,
                   density1 = conditional1,
                   density2 = conditional2) %>%
  pivot_longer(cols = c(density1, density2),
               names_to = 'density',
               values_to = 'y') %>%
  mutate(together = ifelse(density == 'density1', 'never together',
                           'together at least once'))

# plot
(prior_conditional <- ggplot(data)+
    geom_line(aes(x = x, y = y, colour = together),
              linetype = 1,
              linewidth = 1.2)+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density',
                       expand = c(0,0),
                       limits = c(0,15))+
    scale_colour_manual(name = NULL,
                        values = c(rgb(33/255, 145/255, 140/255),
                                   rgb(68/255,1/255,84/255,1)))+
    theme(legend.position = 'right')
)
ggsave(filename = 'prior_conditional.png',
       path = '../outputs/sparse_network_methods_figures/ines_suggestions/',
       plot = prior_conditional, device = 'png', width = 1200, height = 600, units = 'px')

inset_conditional <- ggplot(data)+
  geom_line(aes(x = x, y = y, colour = together, linetype = together),
            linewidth = 1.2)+
  scale_x_continuous(name = NULL)+
  scale_y_continuous(name = NULL,
                     expand = c(0,0),
                     limits = c(0,15))+
  scale_colour_manual(name = NULL,
                      values = c(rgb(33/255, 145/255, 140/255),
                                 rgb(68/255,1/255,84/255,1)))+
  theme(axis.text = element_blank(),
        legend.position = 'none')

# clean up
rm(list = ls()[!ls() %in% c('inset_sri','inset_default','inset_skewed','inset_conditional',
                            'colours')]) ; gc()


######## Conditional Outputs: 2 panel plot of conditional for posterior and edge_vs_sightings ####
#### Posterior distribution                                                                   ####
load('motnp_edgeweights_conditionalprior.RData')
counts_df$together <- ifelse(counts_df$event_count == 0, 'never together', 'together at least once')
counts <- counts_df %>%
  select(dyad_id, count_dyad, together, event_count) %>%
  rename(dyad = dyad_id)
edges <- edges %>%
  filter(chain == 'chain1') %>%
  left_join(counts, by = 'dyad')
rm(edge_binary, edge_samples, edgelist, motnp_ages, nodes, summary, x, n_dyads, n_samples, make_edgelist, plot_network_threshold, i) ; gc()

## plot posterior distribution
(posterior_conditional.1 <- ggplot(data = edges)+
    geom_density(aes(group = dyad, x = edge_draw, colour = together), show.legend = FALSE)+
    stat_density(aes(group = dyad, x = edge_draw, colour = together),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(-0.05, 1.05), expand = c(0,0))+
    scale_y_continuous(name = 'density', limits = c(-2, 72), expand = c(0,0))+
    scale_colour_manual(values = c(rgb(33/255, 145/255, 140/255, 0.05),
                                   rgb(68/255, 1/255, 84/255, 0.05)),
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom', #c(0.5,0.7),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))+
    labs(colour = 'sightings together')
)

## plot a sample of 100 randomly selected dyads
set.seed(15) ; plot_dyads <- sample(1:nrow(counts), size = 100, replace = F)
plot_dyad_ids <- unique(edges$dyad)[plot_dyads]
edges_subset <- edges[edges$dyad %in% plot_dyad_ids,]
edges_subset <- edges_subset[edges_subset$chain == 'chain1',]
(posterior_conditional.2 <- ggplot(data = edges_subset)+
    geom_density(aes(group = dyad, x = edge_draw, colour = together),
                 show.legend = F, linewidth = 0.2)+
    stat_density(aes(group = dyad, x = edge_draw, colour = together),
                 geom = "line", position = "identity", linewidth = 0)+
    scale_x_continuous(name = 'edge weight', limits = c(0,1))+
    scale_y_continuous(name = 'density', limits = c(0, 72))+
    scale_colour_manual(values = colours,
                        aesthetics = 'colour')+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))+
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
)

#### Mean edge weight vs dyad sightings                                                       ####
## create dataframe to plot mean edge weight vs dyad sighting count
averages_conditional <- edges %>%
  group_by(dyad) %>%
  mutate(mean = mean(edge_draw),
         median = median(edge_draw)) %>%
  select(dyad, mean, median) %>%
  distinct() %>%
  rename(dyad_id = dyad) %>%
  left_join(counts_df[,c('dyad_id','node_1_males','node_2_males','event_count','together','count_dyad','count_1','count_2')],
            by = 'dyad_id')

## all median points, no uncertainty
(edgesightings_conditional.1 <- ggplot(averages_conditional,
                                       aes(x = count_dyad, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255, 0.1), size = 0.5, shape = 19)+
    geom_smooth(colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))
)

## all median plots, no uncertainty, split lines by if they were ever seen together
(edgesightings_conditional.2 <- ggplot()+
    geom_point(data = averages_conditional,
               aes(x = count_dyad, y = median),
               colour = rgb(33/255, 145/255, 140/255, 0.1),
               size = 0.5, shape = 19)+
    geom_smooth(data = averages_conditional,
                aes(x = count_dyad, y = median,
                    linetype = as.factor(together)),
                colour = rgb(68/255, 1/255, 84/255))+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'mean weight', limits = c(-0.02,1.02), expand = c(0,0))+
    scale_linetype_manual(name = NULL,
                          values = c(1,6))+
    theme(legend.position = 'bottom',
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8))
)

## all median points, showing uncertainty, split lines by if they were ever seen together
(edgesightings_conditional.3 <- ggplot()+
    geom_point(data = edges,
               aes(x = count_dyad, y = edge_draw),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5, shape = 19)+
    geom_point(data = averages_conditional,
               aes(x = count_dyad, y = median),
               colour = rgb(94/255, 201/255, 98/255, 0.2),
               size = 1, shape = 19)+
    geom_smooth(data = averages_conditional,
                aes(x = count_dyad, y = median,
                    colour = as.factor(together))
    )+
    scale_x_continuous(name = 'total dyad sightings')+
    scale_y_continuous(name = 'edge weight',
                       limits = c(-0.02,1.02),
                       expand = c(0,0))+
    scale_colour_manual(values = colours)+
    theme(legend.position = 'bottom',
          legend.background = element_blank(),
          legend.key.height = unit(4, 'mm'),
          legend.title = element_blank(),
          legend.text = element_text(size = 8))
)

## save outputs
save.image('../outputs/sparse_network_methods_figures/plots_conditional.RData')
rm(list = ls()[ !ls() %in% c('counts_df','n_chains','colours',
                             'posterior_conditional.2','edgesightings_conditional.3',
                             'inset_conditional','fit_edges_motnp')]) ; gc()

#### Merge                                                                                    ####
## create posterior plot without legend to avoid problems with merging legends
post_nolegend <- posterior_conditional.2 + theme(legend.position = 'none')

## combine plots without prior inset
combined <- (post_nolegend + edgesightings_conditional.3) +
  plot_annotation(tag_levels = 'a')
(combined + plot_layout(guides = 'collect') & theme(legend.position = 'bottom'))
ggsave(filename = 'outputs_conditional_noinset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 1700, height = 1000, units = 'px')

## add prior as inset
conditional_posterior_inset <- wrap_elements(posterior_conditional.2 + 
                                               inset_element(inset_conditional,
                                                             left = 0.5, right = 0.99,
                                                             bottom = 0.5, top = 0.99)) +
  theme(legend.position = 'none')
combined <- (conditional_posterior_inset + edgesightings_conditional.3) +
  plot_annotation(tag_levels = 'a')
(combined + plot_layout(guides = 'collect') & theme(legend.position = 'bottom'))
ggsave(filename = 'outputs_conditional_inset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 1680, height = 700, units = 'px')

# clean up
rm(list = ls()[!ls() %in% c('inset_conditional','colours')]) ; gc()

#### Conditional Eigenvector Centrality: 2 panel plot of eigenvector against sighting count and number of non-associations ####
#### Total number of dyads where together = 0                                                                              ####
# get edge samples
edge_samples <- fit_edges_motnp$draws() %>% as.data.frame()
edge_samples <- edge_samples[,(n_chains+1):ncol(edge_samples)]
edge_samples <- edge_samples[,seq(1,ncol(edge_samples)-3, by = 4)]
colnames(edge_samples) <- counts_df$dyad_id

# get adjacency array for default prior
conditional_adjarr <- array(NA, dim = c(length(elephants), length(elephants), 1000),
                            dimnames = list(elephants, elephants, NULL))
for (i in 1:N) {
  dyad_row <- counts_df[i, ]
  conditional_adjarr[dyad_row$id_1, dyad_row$id_2, ] <- edge_samples[, i]
  conditional_adjarr[dyad_row$id_2, dyad_row$id_1, ] <- edge_samples[, i]
}

# make matrix to save eigenvector values to
eigen_conditional <- matrix(NA, nrow = 1000, ncol = length(elephants),
                            dimnames = list(1:1000, elephants))

# calculate eigenvector scores
for(i in 1:nrow(eigen_conditional)){
  conditional_net <- igraph::graph_from_adjacency_matrix(adjmatrix = conditional_adjarr[,,i], diag = FALSE, mode = 'undirected', weighted = T)
  eigen <- as.data.frame(igraph::eigen_centrality(conditional_net)$vector)
  eigen_conditional[i,] <- eigen$`igraph::eigen_centrality(conditional_net)$vector`
}

# make long format data frame
eigen_conditional <- eigen_conditional %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'eigen') %>%
  group_by(id) %>%
  mutate(mean = mean(eigen),
         median = median(eigen)) %>%
  left_join(eigen_sri[,c('id', 'together0', 'count')], by = 'id')

# calculate averages
averages_conditional <- eigen_conditional %>%
  select(id, mean, median, count, together0) %>%
  distinct()

# plot
(eigen0_conditional.1 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages_conditional, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

(eigen0_conditional.2 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = together0, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages_conditional, aes(x = together0, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    geom_smooth(data = eigen_conditional, aes(x = together0, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'dyads together = 0',
                       breaks = c(100,150,200))+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)


(eigen0_conditional.3 <- ggplot(averages_conditional)+
    geom_point(aes(x = together0, y = median, colour = count),
               size = 0.5, shape = 19)+
    geom_smooth(aes(x = together0, y = median),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'dyads together = 0')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))+
    labs(colour = 'number of\nsightings')+
    theme(legend.position = c(0.5,0.2),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.direction = 'horizontal',
          legend.key.height = unit(2.5, 'mm'),
          legend.title = element_text(size = 6),
          legend.text = element_text(size = 4))
)

#### Individual sightings                                                                                                  ####
(eigensightings_conditional.1 <- ggplot()+
   geom_point(data = eigen_conditional, aes(x = count, y = eigen),
              colour = rgb(253/255, 231/255, 37/255, 0.01))+
   geom_point(data = averages_conditional, aes(x = count, y = median),
              colour = rgb(33/255, 145/255, 140/255))+
   scale_x_continuous(name = 'total node sightings')+
   scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

(eigensightings_conditional.2 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01))+
    geom_point(data = averages_conditional, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(data = eigen_conditional, aes(x = count, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)

(eigensightings_conditional.2 <- ggplot()+
    geom_point(data = eigen_conditional, aes(x = count, y = eigen),
               colour = rgb(253/255, 231/255, 37/255, 0.01),
               size = 0.5)+
    geom_point(data = averages_conditional, aes(x = count, y = median),
               colour = rgb(33/255, 145/255, 140/255),
               size = 0.5)+
    geom_smooth(data = eigen_conditional, aes(x = count, y = eigen),
                colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'node centrality', limits = c(-0.02,1.02), expand = c(0,0))
)


(eigensightings_conditional.3 <- ggplot(data = averages_conditional,
                                        aes(x = count, y = median))+
    geom_point(colour = rgb(33/255, 145/255, 140/255))+
    geom_smooth(colour = rgb(68/255,1/255,84/255))+
    scale_x_continuous(name = 'total node sightings')+
    scale_y_continuous(name = 'mean centrality',
                       limits = c(0,1),
                       expand = c(0.02,0.02))
)

# save workspace
save.image('../outputs/sparse_network_methods_figures/eigen_checks_conditional.RData')

#### Merge                                                                                                                 ####
#load('../outputs/sparse_network_methods_figures/eigen_checks.RData')
rm(list = ls()[!ls() %in% c('eigen0_sri.2','eigensightings_sri.2',
                            'eigen0_uniform.2','eigen0_default.2','eigen0_skewed.2',
                            'eigensightings_uniform.2','eigensightings_default.2','eigensightings_skewed.2',
                            'eigen0_conditional.2','eigensightings_conditional.2',
                            'inset_sri','inset_default','inset_skewed','inset_conditional','colours')]) ; gc()

# conditional
(eigen0_conditional.2 + eigensightings_conditional.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_conditional_noinset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 1600, height = 700, units = 'px')

conditional_eigensightings_inset <- wrap_elements(eigensightings_conditional.2 + 
                                                    inset_element(inset_sri,
                                                                  left = 0.75, right = 0.99,
                                                                  bottom = 0.75, top = 0.99))

(conditional_eigensightings_inset + eigen0_conditional.2)+
  plot_annotation(tag_levels = 'a')
ggsave(filename = 'eigen_outputs_conditional_inset.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = last_plot(), device = 'png', width = 1600, height = 700, units = 'px')

#---------------- 11) Supplementary material: simulation -- 6 facet plot of 10 sightings vs 2 sightings, never together vs half together vs always together ####
rm(list = ls()) ; gc()

# define sequence over which to plot
x <- seq(0, 1, length = 100)

# identify probability of each value of x depending on shape of distributions
distributions <- data.frame(x = x,
                            rare2 = dbeta(x = x, shape1 = 0.1, shape2 = 1.9),
                            rare10 = dbeta(x = x, shape1 = 1, shape2 = 9),
                            sometimes2 = dbeta(x = x, shape1 = 1, shape2 = 1),
                            sometimes10 = dbeta(x = x, shape1 = 10, shape2 = 10),
                            usual2 = dbeta(x = x, shape1 = 1.9, shape2 = 0.1),
                            usual10 = dbeta(x = x, shape1 = 9, shape2 = 1)) %>%
  pivot_longer(names_to = 'category', values_to = 'dbeta',
               cols = c("rare2", "rare10",
                        "sometimes2", "sometimes10",
                        "usual2", "usual10")) %>%
  mutate(sightings = ifelse(category == 'rare2' |
                              category == 'sometimes2' |
                              category == 'usual2',
                            2, 10),
         together = ifelse(category == 'rare2' | category == 'rare10',
                           'rare',
                           ifelse(category == 'usual2' | category == 'usual10',
                                  'usual', 'sometimes')),
         labels_together = factor(ifelse(together == 'rare',
                                         'together 10% of the time',
                                         ifelse(together == 'sometimes',
                                                'together 50% of the time',
                                                'together 90% of the time')),
                                  levels = c('together 10% of the time',
                                             'together 50% of the time',
                                             'together 90% of the time')),
         labels_sightings = ifelse(sightings == 2, 'seen twice','seen 10 times')
  )

sri <- data.frame(x = rep(seq(0, 1, length.out = 11), 6),
                  category = rep(unique(distributions$category), each = 11)) %>%
  mutate(sightings = ifelse(category == 'rare2' |
                              category == 'sometimes2' |
                              category == 'usual2',
                            2, 10),
         observed = x*sightings,
         together = ifelse(category == 'rare2' | category == 'rare10',
                           'rare',
                           ifelse(category == 'usual2' | category == 'usual10',
                                  'usual', 'sometimes')),
         true_value = ifelse(together == 'rare', 0.1,
                             ifelse(together == 'sometimes', 0.5, 0.9)),
         labels_together = factor(ifelse(together == 'rare',
                                         'together 10% of the time',
                                         ifelse(together == 'sometimes',
                                                'together 50% of the time',
                                                'together 90% of the time')),
                                  levels = c('together 10% of the time',
                                             'together 50% of the time',
                                             'together 90% of the time')),
         labels_sightings = ifelse(sightings == 2, 'seen twice','seen 10 times')) %>%
  mutate(possible = ifelse(sightings == 10, 'yes',
                           ifelse(x == 0 | x == 0.5 | x == 1, 'yes', 'no'))) %>%
  filter(possible == 'yes') %>%
  select(-possible) %>%
  mutate(sri = dbinom(x = observed, size = sightings, prob = true_value),
         plot_y = 0) %>%
  pivot_longer(cols = c('sri', 'plot_y'), names_to = 'pairing', values_to = 'y_value') %>%
  mutate(pairing = as.integer(as.factor(paste0(x, category))))

true_values <- data.frame(x = rep(c(0.1,0.5,0.9), each = 2),
                          together = rep(unique(distributions$together), each = 2),
                          sightings = rep(c(2,10), 3)) %>%
  mutate(labels_together = factor(ifelse(together == 'rare',
                                         'together 10% of the time',
                                         ifelse(together == 'sometimes',
                                                'together 50% of the time',
                                                'together 90% of the time')),
                                  levels = c('together 10% of the time',
                                             'together 50% of the time',
                                             'together 90% of the time')),
         labels_sightings = ifelse(sightings == 2, 'seen twice','seen 10 times'))

(figure1 <- ggplot()+
    geom_vline(data = true_values, aes(xintercept = x,
                                       colour = 'True edge weight'),
               linewidth = 1.2)+
    geom_line(data = distributions, aes(x = x,
                                        y = dbeta/max(dbeta[dbeta != 'Inf']),
                                        colour = 'Bayesian posterior distribution'),
              linewidth = 1.2,
              key_glyph = 'rect')+
    geom_line(data = sri, aes(x = x, y = y_value, group = pairing,
                              colour = 'SRI probability'),
              linewidth = 1.2)+
    scale_x_continuous(name = 'edge weight')+
    scale_y_continuous(name = 'density', limits = c(0,1))+
    facet_grid(labels_sightings ~ labels_together)+
    scale_color_manual(name='Edge weight',
                       breaks=c('True edge weight', 'SRI probability', 'Bayesian posterior distribution'),
                       values=c('True edge weight'='#fde725',
                                'SRI probability'='#440154',
                                'Bayesian posterior distribution'='#21918c'))+
    theme(legend.position = 'bottom')
)
ggsave(filename = 'example_sri_vs_bison.png',
       path = '../outputs/sparse_network_methods_figures/',
       plot = figure1, device = 'png', width = 2100, height = 1600, units = 'px')

# clean up
rm(list = ls()) ; gc()
dev.off()
