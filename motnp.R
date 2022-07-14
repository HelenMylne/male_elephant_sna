library(tidyverse)
library(lubridate)
library(asnipe)
library(rethinking) # col.alpha
library(igraph)
#### Observation Sessions ####
sessions <- readxl::read_excel("Raw_ElephantDatabase_Youldon210811.xlsx", sheet = 1)
#View(sessions)
#s <- head(sessions, 20)
#s$start <- hms::as_hms(s$Start) ; s <- s[c(1:2,11,3:10)]
s <- sessions[,1:4]
s <- janitor::clean_names(s)
colnames(s)[4] <- 'mins'
s$start <- hms::as_hms(s$start)
s$end <- hms::as_hms(s$end)
hist(s$mins)
#s$dt.start <- paste(s$date, s$start)
#s$session <- ifelse(s$dt.start > lubridate::ymd_hms(paste(s$date, "11:59:00")), 'afternoon',
                    #ifelse(s$end < "14:00:00", 'morning', 'day'))
s$start[1]

s$session <- ifelse(s$start > as.POSIXct("12:00:00",format="%H:%M:%S"), 'afternoon', 'morning')

s$session <- ifelse(s$start > 12*60*60, 'afternoon', 'morning')
a <- s[s$session == 'afternoon',] ; m <- s[s$session == 'morning',]
min(a$start) ; max(m$start)

s$session <- ifelse(s$start >= 22*60*60, 'night', 
                    ifelse(s$start >= 18*60*60, 'evening', 
                           ifelse(s$start >= 12*60*60, 'afternoon',
                                  ifelse(s$end <= 5*60*60, 'night',
                                         ifelse(s$end <= 15*60*60, 'morning','day')))))

length(unique(s$date))
days <- aggregate(s$minutes, by = list(s$date), FUN = sum)
colnames(days) <- c('date','mins')
days$day <- weekdays(days$date)
table(days$day)

#### IDs ####
id.raw <- readxl::read_excel("Raw_ElephantDatabase_Youldon210811.xlsx", sheet = 2)
id <- id.raw[5:849,1:368] # trim off top and empty end columns
colnames(id)[1:10] <- id.raw[4,1:10]
#months <- c('May-2016','June-2016','July-2016','August-2016','September-2016')
#months <- c(lubridate::ym('2016-05', '2016-06','2016-07'))
#colnames(id)[11:13] <- months

#may16 <- dmy('1st May 2016') ; dec17 <- dmy('1st Dec 2017')
#months <- seq(from = may16, to = dec17, by = 'month')
#months

months <- seq(from = dmy('1st May 2016'), to = dmy('1st Dec 2017'), by = 'month')
colnames(id)[11:30] <- months
# still appears as numbers which is fine so long as those numbers are correct... 
dmy_hms(months)
date(months)
colnames(id)[11:30]
date(16922, origin = '1970-01-01')
library(zoo)
as.Date(c(16922,16953,16983,17014,17045,17075,17106,17136,
          17167,17198,17226,17257,17287,17318,17348,17379,17410,17440,17471,17501))

months <- paste('m',as.character(months), sep = '_')
colnames(id)[11:30] <- as.character(months)

weeks <- seq(from = dmy('19th May 2016'), to = dmy('15th Nov 2017'), by = 'weeks')
weeks
colnames(id)[31:108] <- weeks

weeks <- paste('w',as.character(weeks), sep = '_')
colnames(id)[31:108] <- as.character(weeks)

DE <- 26*4+5    # day start column number in Excel data = 109
NO <- 26*14+15  # day end column = 379
NO-DE           # length of days columns = 270

days <- seq(from = dmy('19th May 2016'), to = dmy('15th Nov 2017'), by = 'days')
length(days)    # 546 -- 276 days excluded

days <- id.raw[4,109:368]
View(days)
days.numbers <- data.frame(as.numeric(days)) # all fine except for last 3
colnames(days.numbers) <- 'nums'
days.numbers$dates <- as.Date(days.numbers$nums) # mostly 70 years and 2 days ahead
days.numbers$dates.new <- as.Date(days.numbers$nums - (70*365.25+1))
days.numbers$dates.new[258:260] <- c(days.numbers$dates.new[257]+2,
                                     days.numbers$dates.new[257]+5,
                                     days.numbers$dates.new[257]+6)

colnames(id)[109:368] <- days.numbers$dates.new-0.5
day.check <- as.numeric(colnames(id)[109:368])
days.numbers$dates.check <- as.Date(day.check)
days.numbers$check <- ifelse(days.numbers$dates.check == days.numbers$dates.new-0.5,1,0)
sum(days.numbers$check) == length(days.numbers$check)
colnames(id)[109:368] <- as.character(days.numbers$dates.new-0.5)

days <- paste('d',as.character(days.numbers$dates.new-0.5), sep = '_')
colnames(id)[109:368] <- as.character(days)

id <- janitor::clean_names(id)

#id[is.na(id[11:368])] <- 0
#id[13] <- imputeTS::na_replace(id[11],0)
#id[13] <- ifelse(id[11] == 'Yes', 1,0)
#for(i in 19:20){
#  id[i] <- imputeTS::na_replace(id[i],0)
#  id[i] <- ifelse(id[i] == 'Yes', 1,0)
#}
#sum(id[11:30])

str(id)
id2 <- id[c(1:10, rep(11:368, each=2))]
#id2[12] <- imputeTS::na_replace(id2[12],0)
#id2[12] <- ifelse(id2[12] == 'Yes',1,0)

#for(i in seq(12,40,2)){
#  id2[i] <- ifelse(id2[i] == 'Yes',1,0)
#}
#for(i in seq(46,48,2)){
#  id2[i] <- ifelse(id2[i] == 'Yes',1,0)
#}
#id2[42] <- ifelse(id2[41] == 'Yes',1,0)
#unique(id2[41])

id2[12] <- ifelse(is.na(id2[11]), 0, 1) # proof of concept

for(i in seq(14,100,2)){
  id2[i] <- ifelse(is.na(id2[i-1]), 0, 1) # check1
}
for(i in seq(102,726,2)){
  id2[i] <- ifelse(is.na(id2[i-1]), 0, 1) # check2
}

id2 <- id[c(1:10, rep(11:368, each=2))]
for(i in seq(12,726,2)){
  id2[i] <- ifelse(is.na(id2[i-1]), 0, 1) # all run together
}

for(i in 11:368){
  id[i] <- ifelse(is.na(id[i]), 0, 1) # run in actual data rather than doubled data
}

check <- ifelse(id[11:368] == id2[seq(12,726,2)], 1,0) # check that id matches eye-checked id2
sum(check) # all cells = 1 --> worked as intended! 

str(id)

# months = columns 11:30
id$months <- as.numeric(id$months)      # redo counting - check that R agrees with Excel
id$months2 <- rowSums(id[11:30])        # count up rows
id <- id[c(1:8,369,9:368)]              # inspect by eye -- look like they agree
id$months.check <- ifelse(id$months == id$months2, 1, 0)
sum(id$months.check) == length(id$months) # TRUE (845) --> counts of sightings for each individual per month agree.
max(id$months) # 12
id$months <- as.numeric(id$months)

# weeks = columns 31:108 (32:109 with addition of id$months2)
id$weeks <- as.numeric(id$weeks)        # redo counting - check that R agrees with Excel
id$weeks2 <- rowSums(id[32:109])        # count up rows
id <- id[c(1:10,371,11:370)]            # inspect by eye -- look like they agree
id$weeks.check <- ifelse(id$weeks == id$weeks2, 1, 0)
sum(id$weeks.check) == length(id$weeks) # TRUE (845) --> counts of sightings for each individual per week agree.
max(id$weeks) # 25
id$weeks <- as.numeric(id$weeks)

# days = columns 109:368 (111:370 with addition of id$months2 and id$weeks2)
id$days <- as.numeric(id$days)        # redo counting - check that R agrees with Excel
id$days2 <- rowSums(id[111:370])      # count up rows
id <- id[c(1:12,373,13:372)]          # inspect by eye -- look like they agree
id$days.check <- ifelse(id$days == id$days2, 1, 0)
sum(id$days.check) == length(id$days) # TRUE (845) --> counts of sightings for each individual per month agree.
max(id$days) # 36
id$days <- as.numeric(id$days)

# check that days >= weeks >= months --> can't see an individual in more months than days
id$dmy.check <- ifelse(id$days < id$weeks, 0,
                       ifelse(id$weeks < id$months, 0, 
                              ifelse(id$days < id$months, 0, 1)))
sum(id$dmy.check) == length(id$dmy.check)
which(id$dmy.check != 1)
# correct this one -- at some point will need to go through both IDs and encounters and check they agree
nkash <- id[434,] # seen on 2 weeks but only 1 day -- July 2017, weeks 13-19th July and 20th-26th, day = 20th July. Encounters sheet only listed for 20th July.
nkash[94] <- 0
# BUT -- this is with more columns than actual id sheet will have
id <- id.raw[5:849,1:368]
colnames(id)[1:10] <- id.raw[4,1:10]
colnames(id)[11:30] <- months
colnames(id)[31:108] <- weeks
colnames(id)[109:368] <- days
id <- janitor::clean_names(id)
for(i in 11:368){ id[i] <- ifelse(is.na(id[i]), 0, 1)}    # replace NA with 0 (not seen) and Yes with 1 (seen)
nkash$w_2017_07_13 <- 0
nkash$weeks <- 1
id[434,] <- nkash

id2 <- id.raw[5:849,1:368]
colnames(id2)[1:10] <- id.raw[4,1:10]
colnames(id2)[11:30] <- months
colnames(id2)[31:108] <- weeks
colnames(id2)[109:368] <- days
id2 <- janitor::clean_names(id)
for(i in 11:368){ id2[i] <- ifelse(is.na(id2[i]), 0, 1)}    # replace NA with 0 (not seen) and Yes with 1 (seen)
id2[434,] <- nkash

## replace Excel counts with R counts (though already proven the same)
id2$months <- rowSums(id2[11:30])
id2$weeks <- rowSums(id2[31:108])
id2$days <- rowSums(id2[109:368])

id2$months == id$months
id2$weeks == id$weeks
id2$days == id$days

mwd <- cbind(id2$months, id2$weeks, id2$days)
mwd.raw <- id.raw[5:849, 8:10]
mwd <- cbind(mwd, mwd.raw)
mwd$m <- ifelse(mwd$`1` == mwd$`468...8`, 1, 0)
mwd$w <- ifelse(mwd$`2` == mwd$`468...9`, 1, 0)
mwd$d <- ifelse(mwd$`3` == mwd$`469...10`, 1, 0)
sum(mwd$m) == length(id$id_no)
sum(mwd$d) == length(id$id_no)
sum(mwd$w) == length(id$id_no)
which(mwd$w == 0) # 434
id$name[which(mwd$w == 0)] # Nkash -- so that's the one already corrected elsewhere. 

id <- id[!is.na(id[,1]),]

# age category spreadsheet
ages <- id[,1:6]
ages$age_category <- NA
ages$age_category[1] <- '50+'
write_csv2(ages, na = '', path = '/Users/helen/Documents/ELES!! ❤️🐘🇧🇼/EleSNA/ALERTages.csv', col_names = T)

# new data frame for males
m <- id[id$sex == 'Male',]
m <- m[!is.na(m$id_no),]
m <- m[m$age_class != 'Calf' & m$age_class != 'Juvenile',]
str(m)
hist(m$days, breaks = 30)

#### Encounters ####
encounters <- readxl::read_excel("Raw_ElephantDatabase_Youldon210811.xlsx", sheet = 3)
e <- encounters[c(1:3,8:27)]
e <- janitor::clean_names(e)
colnames(e)[c(4,5,7,8,23)] <- c('gps_s','gps_e','total_id','perc_id','present')

## doesn't work to try and split by breeding herd -- only works for a single column
test <- e %>% separate(present,
                       as.character(c(1:4)),
                       sep = 'B') # splits breeding herds
test <- e %>% separate(present,
                       as.character(c(1:4)),
                       sep = 'B1')
test <- test %>% separate(present, # now multiple columns to split by this
                          as.character(c(1:4)),
                          sep = 'B2')

#unique(e$present) # 456 combinations + NA
d <- e %>% separate(present,
                    as.character(c(1:max(e$total_id, na.rm = T))))
group.size <- paste('e',1:76, sep = '')
colnames(d)[23:98] <- paste('e',1:76, sep = '')

unique(d$e1)
unique(id$herd)

eles <- rep(NA,76)
for(i in 23:98){
  eles[i] <- unique(d[,i])
}
eles[3]

which(d$e1 == 'B1')
d$e1.noB <- ifelse(d$e1 == 'B1', NA,
                   ifelse(d$e1 == 'B2', NA,
                          ifelse(d$e1 == 'B3', NA,
                                 ifelse(d$e1 == 'B4', NA,
                                        ifelse(d$e1 == 'B6', NA,
                                               ifelse(d$e1 == 'B7', NA, d$e1))))))
unique(d$e1)
unique(d$e1.noB)

eles <- d[23:98]
for(i in 1:ncol(eles)) {
  for(j in 1:nrow(eles)){
    eles[j,i] <- ifelse(eles[j,i] == 'B1', NA, eles[j,i])
  }
}
unique(eles$e1)

for(i in 1:ncol(eles)) {
  for(j in 1:nrow(eles)){
    eles[j,i] <- ifelse(eles[j,i] == 'B1', NA,
                        ifelse(eles[j,i] == 'B2', NA,
                               ifelse(eles[j,i] == 'B3', NA,
                                      ifelse(eles[j,i] == 'B4', NA,
                                             ifelse(eles[j,i] == 'B6', NA,
                                                    ifelse(eles[j,i] == 'B7', NA, eles[j,i]))))))
  }
}
unique(eles$e1)
d$e4[2] # B3
unique(eles$e4) # no B3 so it worked to get rid of that too, but there is a blank non-NA cell?
unique(d$e4) # blank non-NA exists before for loop
for(i in 1:ncol(eles)) {
  for(j in 1:nrow(eles)){
    eles[j,i] <- ifelse(eles[j,i] == '', NA, eles[j,i])
  }
}
unique(eles$e4) # no blank cells

eles <- d[23:98]
for(i in 1:ncol(eles)) {
  for(j in 1:nrow(eles)){
    eles[j,i] <- ifelse(eles[j,i] == 'B1', NA,
                        ifelse(eles[j,i] == 'B2', NA,
                               ifelse(eles[j,i] == 'B3', NA,
                                      ifelse(eles[j,i] == 'B4', NA,
                                             ifelse(eles[j,i] == 'B6', NA,
                                                    ifelse(eles[j,i] == 'B7', NA,
                                                           ifelse(eles[j,i] == '', NA, eles[j,i])))))))
  }
}

eles2 <- eles
unite(eles2, present, 1:76, remove = F, sep = '_', na.rm=T)
eles2 <- unite(data = eles2, col = present, 1:76, remove = F, sep = '_', na.rm=T)

e$present <- eles2$present
d <- e %>% separate(col = present, sep = '_',
                    into = as.character(group.size, na.rm = T))
colnames(d)[23:98] <- paste('e',1:76, sep = '')
d$e1 <- ifelse(d$e1 == '', NA, d$e1)

## In short:
encounters <- readxl::read_excel("Raw_ElephantDatabase_Youldon210811.xlsx", sheet = 3)
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

### write csv
### write csv
write_delim(d, 'data_processed/motnp_encounters.csv', na = 'NA', col_names = T, delim = ';')
#write_delim(d, 'data_processed/motnp_encounters2.csv', na = ' ', col_names = T, delim = ';') # doesn't help - still stores as na='NA' and reads back in as a character rather than a blank

#### map encounters -- 4 totally wrong GPS points corrected ####
encounters <- readxl::read_excel("data_raw/Raw_ElephantDatabase_Youldon210811.xlsx", sheet = 3)
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

plot(d$gps_e, d$gps_s) # 2 extreme outliers - one very north, one very east, 2 others outside park
which(d$gps_e > 5e6)   # 232
which(d$gps_s < 5e5)   # 672
which(d$gps_e < 2e6)   # 480 and 624
plot(d$gps_e, d$gps_s, col = col.alpha(rangi2, 0.2), pch = 19)
for(i in c(232,480,624,672)) points(d$gps_e[i], d$gps_s[i], col = 'red', pch = 19)

wrong.locations <- d[c(232,480,624,672),1:8]
d$gps_e[232] <- d$gps_e[232]/10 # value was 25470333 -- 2547033 puts it within cluster of all other points
d$gps_e[480] <- d$gps_e[480]*10 # value was 254665 -- 2546650 puts it within cluster of all other points
d$gps_e[624] <- d$gps_e[624]+800000 # value started 17 rather than 25 (1745269), but then matched format of other east coordinates -- 2545269 puts it within cluster of all other points
d$gps_s[672] <- d$gps_s[672]*10 # value was 175296 -- 1752960 puts it within cluster of all other points

par(xpd = F)
plot(d$gps_e, d$gps_s, col = col.alpha(rangi2, 0.3), pch = 19, xaxt = 'n', yaxt = 'n',
     xlab = 'east', ylab = 'south',
     xlim = c(min(d$gps_e, na.rm=T),max(d$gps_e, na.rm=T)),
     ylim = c(max(d$gps_s, na.rm=T),min(d$gps_s, na.rm=T))) # upside down so more south is bottom
for(i in c(232,480,624,672)) points(d$gps_e[i], d$gps_s[i], col = 'red', pch = 19) # fall within common sighting locations now.
max(d$gps_e, na.rm=T) ; min(d$gps_e, na.rm=T) # 2554583 ; 2542980
max(d$gps_s, na.rm=T) ; min(d$gps_s, na.rm=T) # 1755287 ; 1745598
abline(v = seq(2543000,2555000,1000), col = col.alpha('black',0.3))
abline(h = seq(1745000,1755000,1000), col = col.alpha('black',0.3))
axis(1, at = seq(2543000,2555000,1000), las = 2,
     labels = as.character(seq(25.43,25.55,0.01)))
axis(2, at = seq(1745000,1755000,1000), las = 1,
     labels = as.character(seq(-17.45,-17.55,-0.01)))

plot(d$gps_e, d$gps_s, col = 'white', pch = 19, xaxt = 'n', yaxt = 'n',
     xlab = 'east', ylab = 'south',
     xlim = c(min(d$gps_e, na.rm=T),max(d$gps_e, na.rm=T)),
     ylim = c(max(d$gps_s, na.rm=T),min(d$gps_s, na.rm=T)))
for(i in 1:length(d$date)) points(d$gps_e[i], d$gps_s[i], pch = 19,
                                  col = ifelse(d$type[i] == 'MO','blue',
                                               ifelse(d$type[i] == 'MX','darkgreen',
                                                      ifelse(d$type[i] == 'BH','magenta','black'))))
abline(v = seq(2543000,2555000,1000), col = col.alpha('black',0.3))
abline(h = seq(1745000,1755000,1000), col = col.alpha('black',0.3))
axis(1, at = seq(2543000,2555000,1000), las = 2,
     labels = as.character(seq(25.43,25.55,0.01)))
axis(2, at = seq(1745000,1755000,1000), las = 1,
     labels = as.character(seq(-17.45,-17.55,-0.01)))
legend(x = 25.51, y = -17.46, legend = c('MO','BH','MX','UK'), pch = 19, # doesn't appear
       col = c('blue','magenta','darkgreen','black'))
legend(locator(1), legend = c('MO','BH','MX','UK'), pch = 19,            # does but annoying
       col = c('blue','magenta','darkgreen','black'))

plot(d$gps_e, d$gps_s, col = 'white', pch = 19, xaxt = 'n', yaxt = 'n',
     xlab = 'east', ylab = 'south',
     xlim = c(min(d$gps_e, na.rm=T),max(d$gps_e, na.rm=T)),
     ylim = c(max(d$gps_s, na.rm=T),min(d$gps_s, na.rm=T)))
for(i in 1:length(d$date)) points(d$gps_e[i], d$gps_s[i], col = col.alpha(rangi2,0.4),
                                  pch = ifelse(d$type[i] == 'MO', 1,
                                               ifelse(d$type[i] == 'MX', 2,
                                                      ifelse(d$type[i] == 'BH', 4, 0))))
abline(v = seq(2543000,2555000,1000), col = col.alpha('black',0.3))
abline(h = seq(1745000,1755000,1000), col = col.alpha('black',0.3))
axis(1, at = seq(2543000,2555000,1000), las = 2,
     labels = as.character(seq(25.43,25.55,0.01)))
axis(2, at = seq(1745000,1755000,1000), las = 1,
     labels = as.character(seq(-17.45,-17.55,-0.01)))

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

plot(gps_s ~ gps_e, data = d[d$type == 'MO',],
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

#### asnipe ####
library(asnipe)
data("group_by_individual")
str(gbi)
View(gbi)
data("identified_individuals")

# load data d (Encounters) and id (IDs)
d <- read_csv('data_processed/motnp_encounters.csv')
d <- separate(d, col = 'date;typ;time;gps_s;gps_e;total_elephants;total_id;perc_id;type;calf_male;calf_female;calf_unknown;juv_male;juv_female;juv_unknown;pub_male;pub_female;pub_unknown;adult_male;adult_female;adult_unknown;unknown;e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20;e21;e22;e23;e24;e25;e26;e27;e28;e29;e30;e31;e32;e33;e34;e35;e36;e37;e38;e39;e40;e41;e42;e43;e44;e45;e46;e47;e48;e49;e50;e51;e52;e53;e54;e55;e56;e57;e58;e59;e60;e61;e62;e63;e64;e65;e66;e67;e68;e69;e70;e71;e72;e73;e74;e75;e76',
              into = c('date','typ','time','gps_s','gps_e','total_elephants','total_id','perc_id','type',
                       'calf_male','calf_female','calf_unknown','juv_male','juv_female','juv_unknown',
                       'pub_male','pub_female','pub_unknown','adult_male','adult_female','adult_unknown','unknown',
                       'e1','e2','e3','e4','e5','e6','e7','e8','e9','e10','e11','e12','e13','e14','e15',
                       'e16','e17','e18','e19','e20','e21','e22','e23','e24','e25','e26','e27','e28',
                       'e29','e30','e31','e32','e33','e34','e35','e36','e37','e38','e39','e40','e41',
                       'e42','e43','e44','e45','e46','e47','e48','e49','e50','e51','e52','e53','e54',
                       'e55','e56','e57','e58','e59','e60','e61','e62','e63','e64','e65','e66','e67',
                       'e68','e69','e70','e71','e72','e73','e74','e75','e76'),
              sep = ';')
str(d)
d$date <- as.Date(d$date)
d$time <- as.numeric(d$time)
#d$total_elephants <- as.numeric(d$total_elephants) -- don't do this - 50+ becomes NA
d$total_id <- as.numeric(d$total_id)
d$perc_id <- as.numeric(d$perc_id)

which(is.na(d$total_id))

d$typ <- ifelse(d$typ == 'NA',NA,d$typ)
for(i in 1:nrow(d)){
    d[i,10] <- ifelse(d[i,10] == 'NA',NA,d[i,10])
}
for(j in 11:12){
  for(i in 1:nrow(d)){
    d[i,j] <- ifelse(d[i,j] == 'NA',NA,d[i,j])
  }
}
for(j in 2:ncol(d)){
  for(i in 1:nrow(d)){
    d[i,j] <- ifelse(d[i,j] == 'NA',NA,d[i,j])
  }
}

test <- d
encounters <- readxl::read_excel("data_raw/Raw_ElephantDatabase_Youldon210811.xlsx", sheet = 3)
e <- encounters[c(1:3,8:27)] ; e <- janitor::clean_names(e)
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
d$gps_e[232] <- d$gps_e[232]/10 ; d$gps_e[480] <- d$gps_e[480]*10 ; d$gps_s[672] <- d$gps_s[672]*10 ; d$gps_e[624] <- d$gps_e[624]+800000

d.test <- d==test
which(d.test == 'FALSE') # nothing -- they are the same
which(d == 2549597)      # ok yes the method above should work!
which(d == 2)            # ok yes the method above should work!
# two versions of d (one from raw, one from processed) are now identical and can be used interchangeably

eles <- id[,1:6]

typ1 <- d[d$typ == 1,] ; typNA <- d[is.na(d$typ),]
first <- rbind(typNA, typ1) ; first <- first[!is.na(first$date),]
first <- first[,c(1,3:9,23:98)]
first$encounter <- 1:nrow(first)
first <- first[,c(85,1:84)]

str(first$total_id)

# trying to make a data frame in which there is a new row for each observation of an elephant, so need to repeat the encounter number for the number of identified elephants...
f <- rep(first$encounter, each = first$total_id)
head(f, 30)

which(is.na(first$total_id))
which(is.na(first$e1))

f <- rep(NA, sum(first$total_id, na.rm = T))
for(i in 1:length(f)){
  f[1] <- 1
  f[i] <- ifelse()
}

each <- first$total_id ; each <- each[!is.na(each)]
rep(first$encounter, each = each)
for(i in 1:length(each)){
  f <- rep(first$encounter[i], each = first$total_id[i])
}

#test <- rbind(rep(first, sum(first$total_id, na.rm = T)))
f1 <- first[,c(1:10)]
f2 <- first[,c(1:9,11)]

for(i in 10:85){
  test <- rbind(first[,c(1:9,i)]) # rewrites over test on each loop so just ends with e76 only
}


library(tidyr)

typ1 <- d[d$typ == 1,] ; typNA <- d[is.na(d$typ),]
first <- rbind(typNA, typ1) ; first <- first[!is.na(first$date),]
first <- first[,c(1,3:9,23:98)]
first$encounter <- 1:nrow(first)
first <- first[,c(85,1:84)]
head(first)
first$date <- lubridate::as_date(first$date)

# The arguments to gather():
# - data: Data object
# - key: Name of new key column (made from names of data columns)
# - value: Name of new value column
# - ...: Names of source columns that contain values
# - factor_key: Treat the new key column as a factor (instead of character vector)

eles_long <- gather(data = first, key = number, value = elephant, e1:e76, factor_key = TRUE)
View(eles_long)
eles_long <- eles_long[!is.na(eles_long$elephant),]

eles_asnipe <- eles_long[,c(1,3,4,5,11)]
eles_asnipe$location <- as.factor(paste(eles_asnipe$gps_s, '_', eles_asnipe$gps_e))
eles_asnipe <- eles_asnipe[,c(1,2,5,6)]
data("identified_individuals")
identified_individuals
str(identified_individuals)
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
str(eles_asnipe)
eles_asnipe$Date <- as.numeric(eles_asnipe$Date)
eles_asnipe$ID <- as.factor(eles_asnipe$ID)
str(eles_asnipe)

summary(eles_asnipe$Time)
which(eles_asnipe$Time > 1500)
eles_asnipe$Time <- ifelse(eles_asnipe$Time > 1500, NA, eles_asnipe$Time)
eles_asnipe$Time <- eles_asnipe$Time*86400
summary(eles_asnipe$Time)
hist(eles_asnipe$Time, breaks = 30)
rethinking::dens(eles_asnipe$Time)

a <- get_associations_points_tw(eles_asnipe,
                                time_window = 60, # all individuals seen together have same time, more than 1 minute apart = different groups
                                which_days = NULL,
                                which_locations = NULL)
gbi <- a[[1]]
View(gbi)

head(first)
ggbi <- get_group_by_individual(first, identities = NULL, location = NULL,
                                data_format = "groups") # wow, big nope!

get_network(gbi, data_format = 'GBI', association_index = 'SRI') # 474x474 matrix
length(unique(id$id_no)) # 479 -- am I missing some individuals?

n <- get_network(gbi, data_format = 'GBI', association_index = 'SRI') # 474x474 matrix

sn <- as.data.frame(n)
hist(sn$F9, breaks = 30)

sn$ID1 <- row.names(sn) ; sn <- sn[,c(475, 1:474)]
sna <- gather(data = sn, key = ID2, value = SRI, 2:475, factor_key = TRUE) # make long dyadic
sna$SRI.NA <- ifelse(sna$ID1==sna$ID2, NA, sna$SRI) # remove 0 values for self-association - not sure if this is necessary or not
# note, this has duplicates for all pairs (e.g. F1-F2 and F2-F1)

#dyads.doubled <- data.frame(dyad = paste(sna$ID1,'_',sna$ID2, sep=''),
#                            sri = sna$SRI.NA)
#dyads <- data.frame(dyad = unique(dyads.doubled$dyad)) # doesn't work to remove duplicates because F1_F2 doesn't read the same as F2_F1!

#sna$id1 <- NEED TO WORK OUT HOW TO CONVERT F1 TO F0001 ETC. IN ORDER TO MATCH DATA TO 

hist(sna$SRI.NA, xlab = 'SRI',
     main = 'SRI association values throughout MOTNP')

sna$sex1 <- sna$ID1 ; sna$sex2 <- sna$ID2
sna <- separate(sna, sex1, into = c('sex1','number1'), sep = 1)
sna <- separate(sna, sex2, into = c('sex2','number2'), sep = 1)
unique(sna$sex1); unique(sna$sex2)  # D and N??
sort(unique(sna$ID1)) # D128 and N40 -- D right next to F and N right next to M -- safe to assume that D128 is F128 and N40 is M40?

### start over:

typ1 <- d[d$typ == 1,] ; typNA <- d[is.na(d$typ),]
first <- rbind(typNA, typ1) ; first <- first[!is.na(first$date),]
first <- first[,c(1,3:9,23:98)]
first$encounter <- 1:nrow(first)
first <- first[,c(85,1:84)]
head(first)
first$date <- lubridate::as_date(first$date)
eles_long <- gather(data = first, key = number, value = elephant, e1:e76, factor_key = TRUE)
eles_long <- eles_long[!is.na(eles_long$elephant),]

which(eles_long$elephant == 'D128')
eles_long$elephant[which(eles_long$elephant == 'D128')] <- 'F128'  # row 2297
eles_long$elephant[which(eles_long$elephant == 'N40')] <- 'M40'    # row 1859

eles_asnipe <- eles_long[,c(1,3,4,5,11)]
eles_asnipe$location <- as.factor(paste(eles_asnipe$gps_s, '_', eles_asnipe$gps_e))
eles_asnipe <- eles_asnipe[,c(1,2,5,6)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$Date <- as.numeric(eles_asnipe$Date)
eles_asnipe$ID <- as.factor(eles_asnipe$ID)
eles_asnipe$Time <- ifelse(eles_asnipe$Time > 1500, NA, eles_asnipe$Time)
eles_asnipe$Time <- eles_asnipe$Time*86400
a <- get_associations_points_tw(eles_asnipe,
                                time_window = 60,
                                which_days = NULL,
                                which_locations = NULL)
gbi <- a[[1]]
n <- get_network(gbi, data_format = 'GBI', association_index = 'SRI') # 472x472 matrix
sn <- as.data.frame(n)
sn$ID1 <- row.names(sn) ; sn <- sn[,c(473, 1:472)]
sna <- gather(data = sn, key = ID2, value = SRI, 2:473, factor_key = TRUE) # make long dyadic
sna$SRI.NA <- ifelse(sna$ID1==sna$ID2, NA, sna$SRI)
sna$sex1 <- sna$ID1 ; sna$sex2 <- sna$ID2
sna <- separate(sna, sex1, into = c('sex1','number1'), sep = 1)
sna <- separate(sna, sex2, into = c('sex2','number2'), sep = 1)
unique(sna$sex1); unique(sna$sex2)  # Now just F, M and U
sna <- sna[,c(1:5,7)]
sna <- sna[!is.na(sna$SRI.NA),] ; sna <- sna[,c(1:3,5:6)]


write_delim(eles_long, 'data_processed/motnp_eles_long.csv',col_names = T, delim = ';')
test <- read_csv2('data_processed/motnp_eles_long.csv')
str(test)
test$time <- as.numeric(test$time)

eles_asnipe <- test[,c(1,3,4,5,11)]
eles_asnipe$location <- as.factor(paste(eles_asnipe$gps_s, '_', eles_asnipe$gps_e))
eles_asnipe <- eles_asnipe[,c(1,2,5,6)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$Date <- as.numeric(eles_asnipe$Date)
eles_asnipe$ID <- as.factor(eles_asnipe$ID)
eles_asnipe$Time <- ifelse(eles_asnipe$Time > 1500, NA, eles_asnipe$Time)
eles_asnipe$Time <- eles_asnipe$Time*86400
a <- get_associations_points_tw(eles_asnipe,
                                time_window = 60,
                                which_days = NULL,
                                which_locations = NULL)
gbi <- a[[1]]
n <- get_network(gbi, data_format = 'GBI', association_index = 'SRI') # 472x472 matrix
sn <- as.data.frame(n)
sn$ID1 <- row.names(sn) ; sn <- sn[,c(473, 1:472)]
sna <- gather(data = sn, key = ID2, value = SRI, 2:473, factor_key = TRUE) # make long dyadic
sna$SRI.NA <- ifelse(sna$ID1==sna$ID2, NA, sna$SRI)
sna$sex1 <- sna$ID1 ; sna$sex2 <- sna$ID2
sna <- separate(sna, sex1, into = c('sex1','number1'), sep = 1)
sna <- separate(sna, sex2, into = c('sex2','number2'), sep = 1)
unique(sna$sex1); unique(sna$sex2)  # Now just F, M and U
sna <- sna[,c(1:5,7)]
sna <- sna[!is.na(sna$SRI.NA),] ; sna <- sna[,c(1:3,5:6)]
str(sna)
sna$ID2 <- as.character(sna$ID2)

### code from Dan to remove duplicates
df = read.table(text = "
ID1 ID2
J12 L64
L64 J12
K7 D76
E37 T5
D76 K7
", header=T, stringsAsFactors=F)

no_duplicates_df <- df %>%
  rowwise() %>%      # for each row
  mutate(Composite_ID = paste(sort(c(ID1, ID2)), collapse = "")) %>%  # sort the teams alphabetically and then combine them separating with -
  ungroup() %>%        # forget the row grouping
  distinct_at(vars(Composite_ID))

no_duplicates_df

### remove duplicates
sna_no_duplicates <- sna %>% 
  rowwise() %>% 
  mutate(dyad = paste(sort(c(id1,id2)), collapse = '_')) %>% 
  ungroup() %>% 
  distinct_at(vars(dyad))
# sna$dyad <- paste(sort(c(sna$id1,sna$id2)), collapse = '') -- don't run this code again!
paste(sort(c(sna$id1[1],sna$id2[1])), collapse = '_')
sna$dyad <- NA
for(i in 1:10){
  sna$dyad[i] <- paste(sort(c(sna$id1[i],sna$id2[i])), collapse = '_')
}
for(i in 1:nrow(sna)){
  sna$dyad[i] <- paste(sort(c(sna$id1[i],sna$id2[i])), collapse = '_')
}
# still can't use this to remove the duplicates...?
sna_no_duplicates <- separate(sna_no_duplicates, dyad, into = c('id1','id2'), sep = '_', remove = F)
sna_no_duplicates$sri <- NA
sna_no_duplicates <- separate(sna_no_duplicates, id1, into = c('sex1','number1'), sep = 1, remove = F)
sna_no_duplicates <- separate(sna_no_duplicates, id2, into = c('sex2','number2'), sep = 1, remove = F)
sna_no_duplicates <- sna_no_duplicates[,c(1,2,5,3,6,8)]

match(x = sna_no_duplicates$dyad[1], table = sna$dyad)
match(x = sna_no_duplicates$dyad[10], table = sna$dyad)
match(x = sna_no_duplicates$dyad[100], table = sna$dyad)
match(x = sna_no_duplicates$dyad[1000], table = sna$dyad)
match(x = sna_no_duplicates$dyad[10000], table = sna$dyad)

test <- data.frame(y = 1:5, z = 6:10)
match(x = 7, table = test$z)
test$a <- NA
for(i in 1:5){
  test$a[i] <- test$z[match(x = i, table = test$y)]
}
test2 <- data.frame(y2 = 1:5, z2 = 16:20)
for(i in 1:5){
  test$a[i] <- test2$z2[match(x = i, table = test2$y)]
}

sna_no_duplicates$test <- NA
for(i in 1:10){
  sna_no_duplicates$test[i] <- sna$id1[match(x = sna_no_duplicates$dyad[i], table = sna$dyad)]
}
sna_no_duplicates$test2 <- NA
for(i in 1:10){
  sna_no_duplicates$test2[i] <- sna$id2[match(x = sna_no_duplicates$dyad[i], table = sna$dyad)]
}

for(i in 1:nrow(sna_no_duplicates)) {
  sna_no_duplicates$test[i] <- sna$id1[match(x = sna_no_duplicates$dyad[i], table = sna$dyad)]
}
for(i in 1:nrow(sna_no_duplicates)) {
  sna_no_duplicates$test2[i] <- sna$id2[match(x = sna_no_duplicates$dyad[i], table = sna$dyad)]
}

for(i in 1:nrow(sna_no_duplicates)) {
  sna_no_duplicates$sri[i] <- sna$sri[match(x = sna_no_duplicates$dyad[i], table = sna$dyad)]
}

# check that sna and sna_no_duplicates actually do have the same sri distribution
par(mfrow = c(2,1))
hist(sna$sri, breaks = 1000, ylim = c(0,500), las = 1)
hist(sna_no_duplicates$sri, breaks = 1000, ylim = c(0,250), las = 1)
summary(sna$sri)
summary(sna_no_duplicates$sri)
sum(sna$sri) == 2*sum(sna_no_duplicates$sri)

write_delim(sna, "ALERT_dyads2.csv", delim = " ", na = "NA", append = FALSE,
            col_names = T, quote_escape = "double")
write_delim(sna, "data_processed/motnp_dyads.csv", delim = ";", na = "", append = FALSE,
            col_names = T, quote_escape = "double")


associations <- sna[sna$SRI.NA > 0,]
associations <- associations[!is.na(associations$SRI.NA),]
associations <- sna[sna$SRI.NA > 0,]

par(mfrow = c(2,1))
hist(sna$SRI.NA, xlab = 'SRI', las = 1, ylab = '',
     main = 'SRI association values throughout MOTNP,\nincluding all dyads')
hist(associations$SRI.NA, xlab = 'SRI', las = 1, ylab = '',
     main = 'SRI association values throughout MOTNP,\nafter removing unassociated dyads')

par(mfrow = c(2,1))
hist(sna$SRI.NA, xlab = 'SRI', las = 1, ylab = '',
     main = 'Including all dyads')
hist(associations$SRI.NA, xlab = 'SRI', las = 1, ylab = '',
     main = 'After removing unassociated dyads')

id.counts <- id[,c(1,10)]
id.counts <- data.frame(id = unique(sna$ID1),
                        sightings = NA)

#### Plotting -- most of the working out is stored in the Network Visualisation project, this is mostly cleaned already ####
library(ggnetwork)
net <- intergraph::asNetwork(sn)  #ggnetwork requires a network object, but can't make this work
net <- ggnetwork(net)
ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(aes(color = as.factor(region),size=4)) +
  guides(color=FALSE,size=FALSE)+ #Remove legends
  theme_blank()

library(plotly)
nodes <- colnames(sn)
edges <- as.data.frame(get.edgelist(sn)) # also doesn't work - wrong input format but not sure what it wants

#### nodes data frame
id.raw <- readxl::read_excel("data_raw/Raw_ElephantDatabase_Youldon210811.xlsx", sheet = 2)
id <- id.raw[5:849,1:368] ; id <- id[!is.na(id[,1]),]
colnames(id)[1:10] <- id.raw[4,1:10]
months <- seq(from = dmy('1st May 2016'), to = dmy('1st Dec 2017'), by = 'month')
months <- paste('m',as.character(months), sep = '_') ; colnames(id)[11:30] <- months
weeks <- seq(from = dmy('19th May 2016'), to = dmy('15th Nov 2017'), by = 'weeks')
weeks <- paste('w',as.character(weeks), sep = '_') ; colnames(id)[31:108] <- weeks
days <- id.raw[4,109:368]
days.numbers <- data.frame(as.numeric(days)) ; colnames(days.numbers) <- 'nums'
days.numbers$dates <- as.Date(days.numbers$nums, origin = "1970-01-01")
days.numbers$dates.new <- as.Date(days.numbers$nums - (70*365.25+1), origin = "1970-01-01")
days.numbers$dates.new[258:260] <- c(days.numbers$dates.new[257]+2, days.numbers$dates.new[257]+5, days.numbers$dates.new[257]+6)
days <- paste('d',as.character(days.numbers$dates.new-0.5), sep = '_') ; colnames(id)[109:368] <- days
id <- janitor::clean_names(id)
for(i in 11:368){ 
  id[i] <- ifelse(is.na(id[i]), 0, 1)
}
nkash <- id[434,] ; nkash$w_2017_07_13 <- 0 ; nkash$weeks <- as.character(1) ; id[434,] <- nkash
id$months <- as.numeric(rowSums(id[11:30])) ; id$weeks <- as.numeric(rowSums(id[31:108])) ; id$days <- as.numeric(rowSums(id[109:368]))

ele_nodes <- as.data.frame(id[,c(1:6,10)])

id.test <- read_csv2('data_processed/motnp_id.csv')
test <- id == id.test
which(test == 'FALSE') # 2406
test[2406]
2406/479   # 5.022965 = 5th column, 
2406-5*479 # 11 = 11th row
# yes, this one is meant to be different! WOO! 

ele_nodes <- id.test ; rm(id.test)

#### edges data frame -- code worked out in Network Visualisation project
ele_links <- read.csv("data_processed/ALERT_dyads_separated.csv", header=T, as.is=T)
ele_links <- ele_links[!is.na(ele_links$SRI.NA),]
ele_links$type <- paste(ele_links$sex1,ele_links$sex2, sep = '')
ele_links$type <- ifelse(ele_links$type == 'FM', 'MF',
                         ifelse(ele_links$type == 'FU','UF',
                                ifelse(ele_links$type == 'MU', 'UM', ele_links$type)))
ele_links$num1 <- ele_links$ID1 ; ele_links$num2 <- ele_links$ID2
ele_links <- separate(ele_links, num1, into = c('s1','num1'), sep = 1)
ele_links <- separate(ele_links, num2, into = c('s2','num2'), sep = 1)
ele_links$num1 <- sprintf("%04d", as.integer(ele_links$num1))
ele_links$num2 <- sprintf("%04d", as.integer(ele_links$num2))
ele_links$from <- paste(ele_links$s1,ele_links$num1,sep='')
ele_links$to <- paste(ele_links$s2,ele_links$num2,sep='')
ele_links <- ele_links[,c(12,13,7,3,5,6)]
colnames(ele_links)[4] <- 'weight'

# OR -- all appears to work until reach the stage of plotting and no edge lines appear...
ele_links2 <- read.csv("data_processed/motnp_dyads.csv", header=T, as.is=T)
ele_links2 <- separate(ele_links2, 'id1.id2.sri.sex1.sex2', sep = ';',
                      into = c('id1','id2','sri','sex1','sex2'))
ele_links2$type <- paste(ele_links2$sex1,ele_links$sex2, sep = '')
ele_links2$type <- ifelse(ele_links2$type == 'FM', 'MF',
                         ifelse(ele_links2$type == 'FU','UF',
                                ifelse(ele_links2$type == 'MU', 'UM', ele_links2$type)))
ele_links2$num1 <- ele_links2$id1 ; ele_links2$num2 <- ele_links2$id2
ele_links2 <- separate(ele_links2, num1, into = c('s1','num1'), sep = 1)
ele_links2 <- separate(ele_links2, num2, into = c('s2','num2'), sep = 1)
ele_links2$num1 <- sprintf("%04d", as.integer(ele_links2$num1))
ele_links2$num2 <- sprintf("%04d", as.integer(ele_links2$num2))
ele_links2$from <- paste(ele_links2$s1,ele_links2$num1,sep='')
ele_links2$to <- paste(ele_links2$s2,ele_links2$num2,sep='')
#ele_links2 <- ele_links2[,c(12,13,7,3,5,6)]
ele_links2 <- ele_links2[,c(11,12,6,3:5)]
colnames(ele_links2)[4] <- 'weight'

test <- ele_links == ele_links2
which(test[,1] == 'FALSE') # none
which(test[,2] == 'FALSE') # none
which(test[,3] == 'FALSE') # none
which(test[,4] == 'FALSE') # 48348 FALSE!
which(test[,5] == 'FALSE') # none
which(test[,6] == 'FALSE') # none

ele_links2$weight <- as.numeric(ele_links2$weight)
test <- data.frame(ALERT = ele_links$weight,
                   motnp = ele_links2$weight)
test$check <- ifelse(test$ALERT == test$motnp, 1, 0)
test$r.chk <- ifelse(round(test$ALERT,5) == round(test$motnp,5), 1, 0)
sum(test$r.chk) == length(test$r.chk) # same length -- now they are the same -- motnp was just more precise

#### encounters
encounters <- readxl::read_excel("data_raw/Raw_ElephantDatabase_Youldon210811.xlsx", sheet = 3)
e <- encounters[c(1:3,8:27)] ; e <- janitor::clean_names(e)
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
d <- e %>% separate(col = present, sep = '_', into = as.character(group.size, na.rm = T))
colnames(d)[23:98] <- group.size
d$e1 <- ifelse(d$e1 == '', NA, d$e1)
d$gps_e[232] <- d$gps_e[232]/10 ; d$gps_e[480] <- d$gps_e[480]*10 ; d$gps_s[672] <- d$gps_s[672]*10 ; d$gps_e[624] <- d$gps_e[624]+800000
typ1 <- d[d$typ == 1,] ; typNA <- d[is.na(d$typ),]
first <- rbind(typNA, typ1) ; first <- first[!is.na(first$date),]
first <- first[,c(1,3:9,23:98)]   # remove columns counting number of each demographic group
first$encounter <- 1:nrow(first) ; first <- first[,c(85,1:84)] # index variable of sightings
first$date <- lubridate::as_date(first$date)
eles_long <- gather(data = first, key = number, value = elephant, e1:e76, factor_key = TRUE)
eles_long <- eles_long[!is.na(eles_long$elephant),]
eles_long$elephant[which(eles_long$elephant == 'D128')] <- 'F128'  # row 2297
eles_long$elephant[which(eles_long$elephant == 'N40')] <- 'M40'    # row 1859
eles_long$id_no <- eles_long$elephant
eles_long <- separate(eles_long, id_no, into = c('sex','num'), sep = 1)
eles_long$num <- sprintf("%04d", as.integer(eles_long$num))
eles_long$id_no <- paste(eles_long$sex, eles_long$num, sep='')
eles_long <- eles_long[,c(1:9,11,12,14)]

#### cleaning -- code worked out in Network Visualisation project
rm(days.numbers,id.raw,nkash,days,months,weeks,i,d,e,encounters,typ1,typNA,j,group.size)
ele_nodes <- ele_nodes[!is.na(ele_nodes$age_class),] # remove individuals which are just a number
ele_nodes <- ele_nodes[which(ele_nodes$id_no != 'M0175'),]
ele_links <- ele_links[which(ele_links$from != 'F0158'),]
ele_links <- ele_links[which(ele_links$from != 'F0176'),]
ele_links <- ele_links[which(ele_links$from != 'M0013'),]
ele_links <- ele_links[which(ele_links$from != 'M0125'),]
ele_links <- ele_links[which(ele_links$to != 'F0158'),]
ele_links <- ele_links[which(ele_links$to != 'F0176'),]
ele_links <- ele_links[which(ele_links$to != 'M0013'),]
ele_links <- ele_links[which(ele_links$to != 'M0125'),]

for(i in 1:nrow(ele_nodes)){
  ele_nodes$count[i] <- sum(eles_long$id_no == ele_nodes$id_no[i])
}
ele_nodes$id <- ele_nodes$id_no
ele_nodes <- separate(ele_nodes, id, into = c('sex','num'), sep = 1)
ele_nodes$num <- as.numeric(ele_nodes$num)
ele_nodes$id <- paste(ele_nodes$sex,ele_nodes$num,sep='')

ele_links$age1   <- NA
ele_links$age2   <- NA
ele_links$days1  <- NA
ele_links$days2  <- NA
ele_links$count1 <- NA
ele_links$count2 <- NA
for(i in 1:nrow(ele_links)){
  for(j in 1:nrow(ele_nodes)){
    if(ele_nodes$id_no[j] == ele_links$from[i]) {ele_links$age1[i]   <- ele_nodes$age_class[j]}
    if(ele_nodes$id_no[j] == ele_links$to[i])   {ele_links$age2[i]   <- ele_nodes$age_class[j]}
    if(ele_nodes$id_no[j] == ele_links$from[i]) {ele_links$days1[i]  <- ele_nodes$days[j]}
    if(ele_nodes$id_no[j] == ele_links$to[i])   {ele_links$days2[i]  <- ele_nodes$days[j]}
    if(ele_nodes$id_no[j] == ele_links$from[i]) {ele_links$count1[i] <- ele_nodes$count[j]}
    if(ele_nodes$id_no[j] == ele_links$to[i])   {ele_links$count2[i] <- ele_nodes$count[j]}
  }
}
ele_links$age.dyad <- ifelse(ele_links$age1 == 'Adult',
                             ifelse(ele_links$age2 == 'Adult', 'AA',
                                    ifelse(ele_links$age2 == 'Pubescent', 'AP',
                                           ifelse(ele_links$age2 == 'Juvenile', 'AJ', 'AC'))),
                             ifelse(ele_links$age1 == 'Pubescent',
                                    ifelse(ele_links$age2 == 'Adult', 'AP',
                                           ifelse(ele_links$age2 == 'Pubescent', 'PP',
                                                  ifelse(ele_links$age2 == 'Juvenile', 'PJ', 'PC'))),
                                    ifelse(ele_links$age1 == 'Juvenile',
                                           ifelse(ele_links$age2 == 'Adult', 'AJ',
                                                  ifelse(ele_links$age2 == 'Pubescent', 'PJ',
                                                         ifelse(ele_links$age2 == 'Juvenile', 'JJ', 'CJ'))), 'CC')))

ele_nodes$dem_class <- paste(sep = '',
                             ifelse(ele_nodes$age_class == 'Adult','A',
                                    ifelse(ele_nodes$age_class == 'Pubescent','P',
                                           ifelse(ele_nodes$age_class == 'Juvenile', 'J','C'))),
                             ele_nodes$sex)
ele_links$dem_class1 <- paste(sep = '',
                              ifelse(ele_links$age1 == 'Adult','A',
                                     ifelse(ele_links$age1 == 'Pubescent','P',
                                            ifelse(ele_links$age1 == 'Juvenile','J','C'))),
                              ele_links$sex1)
ele_links$dem_class2 <- paste(sep = '',
                              ifelse(ele_links$age2 == 'Adult','A',
                                     ifelse(ele_links$age2 == 'Pubescent','P',
                                            ifelse(ele_links$age2 == 'Juvenile','J','C'))),
                              ele_links$sex2)

ele_nodes$family[11] <- 'U23 F69'

## all elephants -- code worked out in Network Visualisation project
e_edges <- ele_links[,c(1:4)] ; e_edges <- e_edges[!is.na(e_edges$from),]
e_nodes <- ele_nodes ; e_nodes <- e_nodes[!is.na(e_nodes$id_no),]
e.net <- igraph::graph_from_data_frame(d = e_edges,         # data frame of network edges
                                       vertices = e_nodes,  # node ID and attributes of nodes
                                       directed = F)        # direction of edge strengths

plot(e.net)
graph_attr(e.net, 'layout') <- layout_with_lgl
plot(e.net, edge.width = e_edges$weight,
     #vertex.shape = ifelse(e_nodes$sex == 'Female', 'square',
     #                     ifelse(e_nodes$sex == 'Male','circle','raster')),
     vertex.color= ifelse(e_nodes$age_class == 'Adult','seagreen1',
                          ifelse(e_nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(e_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 4,
     vertex.label = e_nodes$sex)
legend(x = -1.1, y = 1, legend = c('Adult','Pubescent','Juvenile','Calf'), # takes ages to appear
       pch = 19, col = c('seagreen1','skyblue','yellow','magenta'))

## all elephants seen 10 or more times
e10_n <- ele_nodes[ele_nodes$count >9,]
e10_l <- ele_links[ele_links$count1 >9,]
e10_l <- e10_l[e10_l$count2 >9,]
e10.net <- igraph::graph_from_data_frame(d = e10_l, vertices = e10_n, directed = F)
graph_attr(e10.net, 'layout') <- layout_in_circle
plot(e10.net, edge.width = e10_l$weight, edge.curved = 0,
     vertex.color= ifelse(e10_n$age_class == 'Adult','seagreen1',
                          ifelse(e10_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e10_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 4,
     vertex.label = e10_n$sex,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.2))

AP <- e10_l[e10_l$age.dyad == 'AA' | e10_l$age.dyad == 'AP' | e10_l$age.dyad == 'PP',]

MM <- AP[AP$type == 'MM',]
summary(MM$weight)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.06884 0.08211 0.12749 0.45816
FF <- AP[AP$type == 'FF',]
summary(FF$weight)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.05761 0.08904 0.12956 0.93407
MF <- AP[AP$type == 'MF',]
summary(MF$weight)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.00000 0.05854 0.07268 0.11663 0.79888
par(mfrow = c(3,1))
hist(MM$weight, breaks = 30, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Male-Male Interactions')
hist(MF$weight, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Male-Female Interactions')
hist(FF$weight, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Female-Female Interactions')

MM_not0 <- MM[MM$weight > 0,]
MF_not0 <- MF[MF$weight > 0,]
FF_not0 <- FF[FF$weight > 0,]
hist(MM_not0$weight, breaks = 30, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Male-Male Interactions')
hist(MF_not0$weight, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Male-Female Interactions')
hist(FF_not0$weight, breaks = 50, xlab = 'weight', xlim = c(0,1),
     main = 'Adult and Pubescent Female-Female Interactions')

e10_l_plot <- e10_l[e10_l$type == 'MM' | e10_l$type == 'MF' | e10_l$type == 'FF',]
e10_l_plot <- e10_l_plot[e10_l_plot$age.dyad == 'AA' | e10_l_plot$age.dyad == 'AP' |
                           e10_l_plot$age.dyad == 'PP',]
ggplot(data = e10_l_plot, aes(x = weight))+
  geom_histogram(aes(fill = type), colour = 'black')+
  facet_grid(type ~ age.dyad)+
  scale_fill_viridis_d()+
  theme_light()
ggplot(data = e10_l_plot, aes(x = weight))+
  geom_histogram(aes(fill = type), colour = 'black')+
  facet_free(type ~ age.dyad)+
  scale_fill_viridis_d()+
  theme_light()

## all elephants seen 15 or more times
e15_n <- ele_nodes[ele_nodes$count >14,]
e15_l <- ele_links[ele_links$count1 >14,]
e15_l <- e15_l[e15_l$count2 >14,]
e15.net <- igraph::graph_from_data_frame(d = e15_l, vertices = e15_n, directed = F)
graph_attr(e15.net, 'layout') <- layout_in_circle
plot(e15.net, edge.width = e15_l$weight, edge.curved = 0,
     vertex.color= ifelse(e15_n$age_class == 'Adult','seagreen1',
                          ifelse(e15_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e15_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 6,
     vertex.label = e15_n$sex,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.2))

## all elephants seen 20 or more times
e20_n <- ele_nodes[ele_nodes$count >19,]
e20_l <- ele_links[ele_links$count1 >19,]
e20_l <- e20_l[e20_l$count2 >19,]
e20.net <- igraph::graph_from_data_frame(d = e20_l, vertices = e20_n, directed = F)
graph_attr(e20.net, 'layout') <- layout_in_circle
plot(e20.net, edge.width = e20_l$weight, edge.curved = 0,
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$dem_class,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.2))
# mostly see no patterns from AM nodes, but one thick line from AM to CU
plot(e20.net, edge.width = e20_l$weight, edge.curved = 0,
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$id,
     edge.color = ifelse(e20_l$from == 'U0009'|e20_l$to == 'U0009', 'red',
                         rgb(0.1,0.1,0.1,alpha=0.2)))
plot(e20.net, edge.width = e20_l$weight, edge.curved = 0,
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$id,
     edge.color = ifelse(e20_l$from == 'M0005'|e20_l$to == 'M0005', 'red',
                         rgb(0.1,0.1,0.1,alpha=0.2)))
plot(e20.net, edge.width = e20_l$weight, edge.curved = 0,
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$id,
     edge.color = ifelse(e20_l$from == 'U0005'|e20_l$to == 'U0005', 'red',
                         rgb(0.1,0.1,0.1,alpha=0.2)))
plot(e20.net, edge.width = e20_l$weight, edge.curved = 0,
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$id,
     edge.color = ifelse(e20_l$from == 'M0006'|e20_l$to == 'M0006', 'red',
                         rgb(0.1,0.1,0.1,alpha=0.2)))
# never mind -- it's going between U5 and M6 (both Marianne's family), not U9 and M5!

plot(e20.net, edge.width = e20_l$weight, edge.curved = 0,
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$dem_class,
     edge.color = ifelse(e20_l$age.dyad == 'AA','seagreen1',
                         ifelse(e20_l$age.dyad == 'AP','seagreen3',
                                ifelse(e20_l$age.dyad == 'AJ','green',
                                       ifelse(e20_l$age.dyad == 'AC','green',
                  ifelse(e20_l$age.dyad == 'PP','skyblue',
                         ifelse(e20_l$age.dyad == 'PJ','purple',
                                ifelse(e20_l$age.dyad == 'PC', 'purple', 'magenta'))))))))
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
                         rgb(1,0,0,alpha=0.5),
                         rgb(0.1,0.1,0.1,alpha=0.2)))

df <- data.frame(v1 = rep(1:5,2), v2 = rep(1:2, each = 5))
plot(df, ylim = c(0,3), cex = 2, pch = 19, col = 'blue')
points(df, cex = 5, pch = 19,
     col = ifelse(df$v2 == 1, rgb(1,0,0,alpha=0.5),rgb(0.1,0.1,0.1,alpha=0.2) ))
plot(df, ylim = c(0,3), cex = 2, pch = 19, col = 'blue')
points(df, cex = 5, pch = 19,
       col = ifelse(df$v2 == 1, rgb(1,0,0,alpha=0.2),rgb(0.1,0.1,0.1,alpha=0.2) ))

## males only -- code worked out in Network Visualisation project
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
plot(m.net, edge.width = m_edges$weight*10)
plot(m.net, edge.width = m_edges$weight*10,
     #margin = c(-2,-2,-2,-2),
     vertex.color = ifelse(m_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.shape = ifelse(m_nodes$count < 5, 'square', 'circle'),
     vertex.label = m_nodes$id,
     vertex.size = ifelse(m_nodes$count < 5, 3, 8))

## males with 10 or more sightings only
males10 <- males[males$count > 9,]
males10 <- males10[!is.na(males10$id_no),]

m_links <- ele_links[ele_links$type == 'MM',]
m_links <- m_links[m_links$age1 == 'Adult' | m_links$age1 == 'Pubescent',]
m_links <- m_links[m_links$age2 == 'Adult' | m_links$age2 == 'Pubescent',]
m10_links <- m_links[m_links$count1 > 9,]
m10_links <- m10_links[m10_links$count2 > 9,]

m10_edges <- m10_links[,c(1:4)] ; m10_edges <- m10_edges[!is.na(m10_edges$from),]
m10_nodes <- males10 ; m10_nodes <- m10_nodes[!is.na(m10_nodes$id_no),]

m10.net <- igraph::graph_from_data_frame(d = m10_edges,         # data frame of network edges
                                         vertices = m10_nodes,  # node ID and attributes of nodes
                                         directed = F)          # direction of edge strengths
graph_attr(m10.net, 'layout') <- layout_with_lgl
set.seed(11) ; plot(m10.net, #edge.width = m10_edges$weight,
                    #margin = c(-2,-2,-2,-2),
                    vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
                    vertex.label = m10_nodes$id,
                    vertex.size = m10_nodes$count/2,
                    edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
                    edge.curved = 0)

graph_attr(m10.net, 'layout') <- layout_in_circle
set.seed(11) ; plot(m10.net, edge.width = m10_edges$weight,
                    vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
                    vertex.label = m10_nodes$id,
                    vertex.size = m10_nodes$count/2,
                    #vertex.label.degree = ifelse(m10_nodes$num <= 18 | m10_nodes$num >= 160, 0,
                    #                             ifelse(m10_nodes$num > 18 & m10_nodes$num <= 47,
                    #                                    '-pi/2',
                    #                                    ifelse(m10_nodes$num > 47 & m10_nodes$num <= 102,
                    #                                           'pi', 'pi/2'))),
                    edge.color = ifelse(m10_edges$from == 'M0078',
                                        ifelse(m10_edges$to == 'M0009',
                                               'red', rgb(0.1,0.1,0.1,alpha=0.5)),
                                        ifelse(m10_edges$from == 'M0009',
                                               ifelse(m10_edges$to == 'M0078',
                                                      'red', rgb(0.1,0.1,0.1,alpha=0.5)),
                                               rgb(0.1,0.1,0.1,alpha=0.5))))
set.seed(11) ; plot(m10.net, edge.width = m10_edges$weight,
                    vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
                    vertex.label = m10_nodes$id,
                    vertex.size = m10_nodes$count/2,
                    edge.color = ifelse(m10_edges$from == 'M0078' & m10_edges$to == 'M0009', 'red', 
                                        ifelse(m10_edges$from == 'M0009' & m10_edges$to == 'M0078', 'red',
                                               rgb(0.1,0.1,0.1,alpha=0.5))))

set.seed(11) ; plot(m10.net, edge.width = m10_edges$weight,
                    vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
                    vertex.label = m10_nodes$id,
                    vertex.size = m10_nodes$count/2,
                    edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
                    edge.curved = 0)

## ignore this
graph_attr(m10.net, 'layout') <- layout_with_kk(m10.net, dim = 3) # dim = 2 gives circular layout, dim = 3 gives a virus but with overlapping nodes
plot(m10.net, edge.width = m10_edges$weight,
     vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.label = m10_nodes$id,vertex.size = m10_nodes$count/2,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5), edge.curved = 0)

length(m10_nodes$id_no)/5 # 15
(m10.mat <- matrix(data = m10_nodes$id_no, nrow = 15, ncol = 5, byrow = T))
(m10.mat <- matrix(data = 1:75, nrow = 15, ncol = 5, byrow = T))
(m10.mat <- matrix(data = m10_nodes$id_no, nrow = 25, ncol = 3, byrow = T))
graph_attr(m10.net, 'layout') <- layout_with_kk(m10.net, dim = 3,
                                                coords = m10.mat)

## females -- just to see if the stronger ties are more obvious (10 or more sightings only) -- code worked out in Network Visualisation project
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
plot(f.net, edge.width = f_edges$weight*2)
plot(f.net, edge.width = f_edges$weight*2,
     vertex.color = ifelse(f_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f_nodes$id,
     vertex.size = f_nodes$count)

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
plot(f10.net, edge.width = f10_edges$weight*2)
plot(f10.net, edge.width = f10_edges$weight*2,
     vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f10_nodes$id,
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


### B7
b7.n <- ele_nodes[ele_nodes$herd == "B7",] ; b7.n <- b7.n[!is.na(b7.n$id_no),]
b7.e <- ele_links[ele_links$from == 'F0052' | ele_links$from == 'F0060' | ele_links$from == 'F0098' | 
                    ele_links$from == 'U0017' | ele_links$from == 'U0021',]
b7.e <- b7.e[b7.e$to == 'F0052' | b7.e$to == 'F0060' | b7.e$to == 'F0098' | 
               b7.e$to == 'U0017' | b7.e$to == 'U0021',]
b7 <- graph_from_data_frame(d = b7.e, vertices = b7.n, directed = F)
plot(b7, edge.curved = 0, vertex.label = b7.n$id,
     edge.width = b7.e$weight*8)

### Kamada-Kawai layout (plots using code from ALERT_clean, after cutting down #edges shown)
#graph_attr(b7, 'layout') <- layout_with_kk(b7, dim = 2) # puts individuals with the thickest lines on the outside and lowest weights in the middle
#plot(b7, edge.curved = 0.3, vertex.label = b7.n$id, edge.width = b7.e$weight*8)

#(b7.mat <- matrix(data = c(0,2, 2,2, 1,1, 0,0, 2,0),
#                  nrow = 5, ncol = 2, byrow = T))
#graph_attr(m10.net, 'layout') <- layout_with_kk(m10.net, dim = 2, coords = b7.mat) # nope, thought this might tell the nodes where to sit. says it requires 2 columns and n rows.

#graph_attr(m10.net, 'layout') <- g.layout_kamada_kawai(seed = 5)
l <- layout_with_kk(b7)
plot(b7, layout = l,
     edge.curved = 0, vertex.label = b7.n$id, edge.width = b7.e$weight*3)

plot(e20.net, layout = layout_with_kk(e20.net), # goes circular again
     edge.curved = 0,
     edge.width = ifelse(e20_l$weight < 0.2, 0, e20_l$weight),
     edge.color = ifelse(e20_l$weight < 0.2, 0, rgb(0.1,0.1,0.1,alpha=0.4)),
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 8,
     vertex.label = e20_n$dem_class)
plot(e20.net, layout = layout_with_kk(e20.net, dim = 3), # spherical, some overlapping
     edge.curved = 0,
     edge.width = ifelse(e20_l$weight < 0.2, 0, e20_l$weight),
     edge.color = ifelse(e20_l$weight < 0.2, 0, rgb(0.1,0.1,0.1,alpha=0.4)),
     vertex.color= ifelse(e20_n$age_class == 'Adult','seagreen1',
                          ifelse(e20_n$age_class == 'Pubescent','skyblue',
                                 ifelse(e20_n$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 10,
     vertex.label = e20_n$dem_class)

plot(m10.net, layout = layout_with_kk(m10.net, dim = 2),
     edge.width = ifelse(m10_edges$weight < 0.19, 0, m10_edges$weight),
     edge.color = ifelse(m10_edges$weight < 0.19, 0, rgb(0.1,0.1,0.1,alpha=0.6)),
     vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.label = m10_nodes$id,
     vertex.size = m10_nodes$count/2,
     edge.curved = 0)
plot(m10.net, layout = layout_with_kk(m10.net, dim = 3),
     edge.width = ifelse(m10_edges$weight < 0.19, 0, m10_edges$weight),
     edge.color = ifelse(m10_edges$weight < 0.19, 0, rgb(0.1,0.1,0.1,alpha=0.6)),
     vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.label = m10_nodes$id,
     vertex.size = m10_nodes$count/2,
     edge.curved = 0)

plot(f10.net, layout = layout_with_kk(f10.net, dim = 2),
     edge.width = ifelse(f10_edges$weight < 0.2, 0, f10_edges$weight),
     edge.color = ifelse(f10_edges$weight < 0.2, 0, rgb(0.1,0.1,0.1,alpha=0.5)),
     vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f10_nodes$id,
     vertex.size = f10_nodes$count/2,
     edge.curved = 0)
plot(f10.net, layout = layout_with_kk(f10.net, dim = 3),
     edge.width = ifelse(f10_edges$weight < 0.2, 0, f10_edges$weight),
     edge.color = ifelse(f10_edges$weight < 0.2, 0, rgb(0.1,0.1,0.1,alpha=0.5)),
     vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f10_nodes$id,
     vertex.size = f10_nodes$count/2,
     edge.curved = 0)


### write out so don't have to recreate the ele_links data set every time
write_delim(ele_links, path = 'data_processed/ALERT_elelinks.csv', delim = ';', col_names = T)
write_delim(ele_nodes, path = 'data_processed/ALERT_elenodes.csv', delim = ';', col_names = T)
rm(e_edges, e_nodes, e.net, eles_long, f_edges, f_links, f_nodes, f10_edges, f10_links, f10_nodes, f10.net, females, females10, first, id, m_edges, m_links, m_nodes, m.net, m10_edges, m10_links, m10_nodes, m10.net, males, males, males10, i, j)

# test usable
test.links <- read_csv('data_processed/ALERT_elelinks.csv')
test.links <- separate(test.links, col = 'from;to;type;weight;sex1;sex2;age1;age2;days1;days2;count1;count2;age.dyad;dem_class1;dem_class2', into = c('from','to','type','weight','sex1','sex2','age1','age2','days1','days2','count1','count2','age.dyad','dem_class1','dem_class2'), sep = ';')

str(test.links)
test.links$weight <- as.numeric(test.links$weight)
test.links$count1 <- as.numeric(test.links$count1)
test.links$count2 <- as.numeric(test.links$count2)
test.links$days1 <- as.numeric(test.links$days1)
test.links$days2 <- as.numeric(test.links$days2)
test.links$age.dyad <- as.factor(test.links$age.dyad)
test.links$dem_class1 <- as.factor(test.links$dem_class1)
test.links$dem_class2 <- as.factor(test.links$dem_class2)
str(test.links)

test.nodes <- read_csv('data_processed/ALERT_elenodes.csv')
test.nodes <- separate(test.nodes, col = 'id_no;herd;name;age_class;family;days;count;sex;num;id;dem_class', into = c('id_no','herd','name','age_class','family','days','count','sex','num','id', 'dem_class'), sep = ';')

str(test.nodes)
test.nodes$days <- as.numeric(test.nodes$days)
test.nodes$count <- as.numeric(test.nodes$count)

# remove F157
test.links <- test.links[test.links$from  != 'F0157',]
test.links <- test.links[test.links$to    != 'F0157',]
test.nodes <- test.nodes[test.nodes$id_no != 'F0157',]

e_edges <- test.links[,c(1:4)] ; e_edges <- e_edges[!is.na(e_edges$from),]
e_nodes <- test.nodes ; e_nodes <- e_nodes[!is.na(e_nodes$id_no),]
e.net <- igraph::graph_from_data_frame(d = e_edges,         # data frame of network edges
                                       vertices = e_nodes,  # node ID and attributes of nodes
                                       directed = F)        # direction of edge strengths
plot(e.net)

## males only
males <- test.nodes[test.nodes$sex=='M'&test.nodes$age_class!='Calf'&test.nodes$age_class!='Juvenile',]
males <- males[!is.na(males$id_no),]
m_links <- test.links[test.links$type == 'MM',]
m_links <- m_links[m_links$age1 == 'Adult' | m_links$age1 == 'Pubescent',]
m_links <- m_links[m_links$age2 == 'Adult' | m_links$age2 == 'Pubescent',]
m_edges <- m_links[,c(1:4)] ; m_edges <- m_edges[!is.na(m_edges$from),]
m_nodes <- males ; m_nodes <- m_nodes[!is.na(m_nodes$id_no),]
m.net <- igraph::graph_from_data_frame(d = m_edges,         # data frame of network edges
                                       vertices = m_nodes,  # node ID and attributes of nodes
                                       directed = F)        # direction of edge strengths
graph_attr(m.net, 'layout') <- layout_with_lgl
plot(m.net, edge.width = m_edges$weight*3,
     vertex.color = ifelse(m_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
     vertex.shape = ifelse(m_nodes$count < 5, 'square', 'circle'),
     vertex.label = m_nodes$id,
     vertex.size = ifelse(m_nodes$count < 5, 3, 8),
     edge.curved = 0)

## males with 10 or more sightings only
males10 <- males[males$count > 9,] ; males10 <- males10[!is.na(males10$id_no),]
m10_links <- m_links[m_links$count1 > 9,] ; m10_links <- m10_links[m10_links$count2 > 9,]
m10_edges <- m10_links[,c(1:4)] ; m10_edges <- m10_edges[!is.na(m10_edges$from),]
m10_nodes <- males10 ; m10_nodes <- m10_nodes[!is.na(m10_nodes$id_no),]
m10.net <- igraph::graph_from_data_frame(d = m10_edges,         # data frame of network edges
                                         vertices = m10_nodes,  # node ID and attributes of nodes
                                         directed = F)          # direction of edge strengths
graph_attr(m10.net, 'layout') <- layout_with_lgl
set.seed(11) ; plot(m10.net, edge.width = m10_edges$weight,
                    vertex.color = ifelse(m10_nodes$age_class == 'Adult', 'seagreen1', 'skyblue1'),
                    vertex.label = m10_nodes$id,
                    vertex.size = m10_nodes$count/2,
                    edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
                    edge.curved = 0)

## females -- just to see if the stronger ties are more obvious (10 or more sightings only)
females <- anti_join(test.nodes, males)
females <- females[!is.na(females$id_no),]
f_links <- test.links[test.links$dem_class1 != 'AM' & test.links$dem_class2 != 'AM',]
f_links <- f_links[f_links$dem_class1 != 'PM' & f_links$dem_class2 != 'PM',]
f_links <- f_links[f_links$from != 'F0157' | f_links$to != 'F0157',]
f_edges <- f_links[,c(1:4)] ; f_edges <- f_edges[!is.na(f_edges$from),]
f_nodes <- females[females$id_no != 'F0157',] ; f_nodes <- f_nodes[!is.na(f_nodes$id_no),]
f.net <- igraph::graph_from_data_frame(d = f_edges,         # data frame of network edges
                                       vertices = f_nodes,  # node ID and attributes of nodes
                                       directed = F)        # direction of edge strengths

length(unique(f_edges$from))  # 246
length(unique(f_edges$to))    # 246
length(unique(f_nodes$id_no)) # 245
f <- sort(unique(f_edges$from))  # 246
n <- sort(unique(f_nodes$id_no)) # 245
which(n!=f) # 152 onwards
f[152] # F157 has snuck her way back in again...
which(ele_links$from == 'F0157')
which(ele_nodes$id_no == 'F0157')

graph_attr(f.net, 'layout') <- layout_with_lgl
plot(f.net, edge.width = f_edges$weight*2,
     vertex.color = ifelse(f_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f_nodes$id,
     vertex.size = f_nodes$count,
     edge.curved = 0)

females10 <- females[females$count > 9,]  ; females10 <- females10[!is.na(females10$id_no),]
f10_links <- f_links[f_links$count1 > 9,] ; f10_links <- f10_links[f10_links$count2 > 9,]
f10_edges <- f10_links[,c(1:4)] ; f10_edges <- f10_edges[!is.na(f10_edges$from),]
f10_nodes <- females10 ; f10_nodes <- f10_nodes[!is.na(f10_nodes$id_no),]

f10.net <- igraph::graph_from_data_frame(d = f10_edges,         # data frame of network edges
                                         vertices = f10_nodes,  # node ID and attributes of nodes
                                         directed = F)          # direction of edge strengths
graph_attr(f10.net, 'layout') <- layout_with_lgl
plot(f10.net, edge.width = f10_edges$weight*2,
     vertex.color = ifelse(f10_nodes$age_class == 'Adult','seagreen1',
                           ifelse(f10_nodes$age_class == 'Pubescent','skyblue',
                                  ifelse(f10_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.label = f10_nodes$id,
     vertex.size = f10_nodes$count/1.5,
     edge.curved = 0,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5))

#### Summaries ####
table(ele_nodes$dem_class)
summary(males$count)
hist(males$count, breaks = 46)
ele_nodes$dem_class2 <- ifelse(ele_nodes$age_class == 'Calf','C',
                               ifelse(ele_nodes$age_class == 'Juvenile','J', ele_nodes$dem_class))
ele_nodes$dem_class2 <- factor(ele_nodes$dem_class2, levels = c('AM','AF','PM','PF','PU','J','C'))
boxplot(ele_nodes$count ~ ele_nodes$dem_class2, las = 1, pch = 16,
        xlab = 'demographic group',
        ylab = 'count of sightings',
        col = c('seagreen1','seagreen1','skyblue','skyblue','skyblue','yellow','magenta'),
        outcol = c('seagreen1','seagreen1','skyblue','skyblue','skyblue','yellow','magenta'))
points(pch = 1, ele_nodes$count ~ ele_nodes$dem_class2)

ele_nodes$count[which(ele_nodes$name == 'Emily-Ruth'),]

















#### Markdown TAP 1 ####
tinytex::install_tinytex()
?anti_join
layout_in_circle