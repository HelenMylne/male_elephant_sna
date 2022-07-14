#### Information ####
# Script to process and plot association data from Makgadikgadi Pans National Park, Bostwana
# Data collected: 1st November 2019 to 5th August 2021
# Collected by: 
# Data supplied by: Dr Kate Evans, Elephants for Africa, 19th October 2021

#### Set up ####
library(tidyverse)  # data manipulation
library(lubridate)  # sort date columns out
library(zoo)        # sort date columns out
library(asnipe)     # generating networks
library(rethinking) # col.alpha
library(igraph)     # plotting networks

#### Data import ####
d <- readxl::read_excel("data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx", sheet = 1)
str(d)
length(unique((d$ID_Elephant_Visual)))   # 13429 -- all unique = unique number per elephant sighting
length(unique((d$Elephant_Sighting_ID))) # 5162 -- unique number per group sighting
length(unique((d$Elephant_ID)))          # 6646 identified individuals
table(d$Elephant_ID)                     # 5062 unidentified, rest mostly seen once (max 22). T=temporaryID, B=checked.
e <- d[d$Elephant_ID != '-',]
barplot(table(e$Elephant_ID),horiz = T, las = 1, xlim = c(0,25))
table(d$In_The_Field_Elephant_Code)      # Unique code per elephant per day, max=FFFF=84?
length(unique((d$Date)))                 # 985 days
str(d$Time)                              # datetime rather than just date
d$time_cat <- lubridate::hour(d$Time)
e <- separate(d, Time, into = c('wrong_date','time'), sep = ' ')
head(e)
e$time <- hms::as_hms(e$time)
d <- e[,c(1:5,7,33,8:32)]
str(d)
d$Age_Range_ID <- as.factor(d$Age_Range_ID)
summary(d$Age_Range_ID)     # 2=1-4y, 3=5-9y, 4=10-15y, 5=16-20y, 6=21-25y, 7=26-35y, 8=>35y, 10=??
summary(d$Age)              # all NA
d$Activity_ID <- as.factor(d$Activity_ID)   
summary(d$Activity_ID)      # 0-10, not integer, 3 and 6 most common, no ethogram in doc from Kate, 12298 NA
table(d$In_musth)           # 118 Y, 757 UK, and then I'm assuming N=n=M=0 for not-musth = 12554
summary(d$Reaction_Indices) # -4=run away, 1=don't notice, 2=unbothered, 3=agitated, 4=charge, range 1-2
d$Distance_To_Observer <- as.numeric(d$Distance_To_Observer)
summary(d$social_id)        # range 1-3, mostly 2 or 3
barplot(table(d$social_id), las = 1, xlab = 'social ID value', ylab = 'count') ; abline(h = 0)
length(which(is.na(d$social_id) == TRUE))
table(d$elephant_sighting_id, d$social_id)
summary(d$Distance_To_Observer)  # 1-5252, use a range finder (GPS elephants = range+GPS_car+bearing)
table(d$Habitat_ID)         # 0-8 but not integer, UK=4, NA=12346
table(d$Vegetal_Species_ID) # only 31 filled
table(d$Feeding_Type_ID)    # only 57 filled (37 browse, 14 graze, 6 fruit)
summary(d$Pictures_Taken)   # 4181 no pictures, 9248 with pictures
table(d$Sick)               # 38 yes, 12602 no (+17*0+1*M), 771 unknown
d$Physical_Condition_ID <- as.factor(d$Physical_Condition_ID)
summary(d$Physical_Condition_ID) # 1=emaciated, 2=very thin, 3=thin, 4=good, 5=fat. 971 UK, 829 NA, max 5, mode 3, min 1.5

# L/R_TGS_**_ID -- TGS = temporal gland secretion, L/R = Left/Right
table(d$L_TGS_Length_ID)    # 0=4125, 1=378, 2=387, 3=443, 4=515, UK=7148
table(d$L_TGS_Width_ID)     # t,m,w=thin,medium,wide?? 0=2455, m=863, t=692, UK=4234, w=170
table(d$L_TGS_Age_ID)       # n=new=1249, o=old=463, 0=2461, UK=4231

table(d$R_TGS_Length_ID)    # 0=4171, 1=361, 2=364, 3=467, 4=497, UK=7136
table(d$R_TGS_Width_ID)     # t=631, m=882, w=186, 0=2480, UK=4277
table(d$R_TGS_Age_ID)       # n=1180, o=515, 0=2478, UK=4275

table(d$TG_Swelling_ID)     # 0=8557 no swelling of temporal glands, 1=57 (slightly swollen), 2=42 (swollen), 3=36 (very swollen), 93=unknown
summary(d$Hind_Foot_Length_Average)          # 154 T
summary(d$Front_Foot_Circumference_Average)  # 3T
summary(d$Height)           # all NA
table(d$Sex_ID)             # 1=13359 (males), 2=9 (?), 3=5 (?)

str(d)
d <- d[,c(1:15,18:27,31:32)]
head(d)
d$Date <- lubridate::as_date(d$Date)
d <- janitor::clean_names(d)

#### Clean data import ####
d <- readxl::read_excel("data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx", sheet = 1)
d$time_cat <- lubridate::hour(d$Time)
d <- separate(d, Time, into = c('wrong_date','time'), sep = ' ')
d$time <- hms::as_hms(d$time)
d <- d[,c(1:5,7,33,8:32)]
d$Age_Range_ID <- as.factor(d$Age_Range_ID)
d$Activity_ID <- as.factor(d$Activity_ID)   
d$Distance_To_Observer <- as.numeric(d$Distance_To_Observer)
d$Physical_Condition_ID <- as.factor(d$Physical_Condition_ID)
d <- d[,c(1:15,18:27,31:32)]
d$Date <- lubridate::as_date(d$Date)
d <- janitor::clean_names(d)

#### taking a look ####
table(d$elephant_id)
hist(table(d$elephant_sighting_id), breaks = 30, las = 1)
length(which(table(d$elephant_sighting_id) == 1)) # 2400 elephants seen singly
length(which(table(d$elephant_sighting_id) == 2)) # 1017 elephants seen as pairs

table(d$age_range_id)

summary(d$reaction_indices)
d$reaction_indices <- factor(d$reaction_indices, levels = c(1:3))
table(d$reaction_indices)
barplot(table(d$reaction_indices), xlab = 'reaction index score', ylab = 'count',
        las = 1, ylim = c(0,500), xlim = c(0,4))
abline(h = 0, lwd = 2) ; axis(2, at = seq(0,500,50), las = 1)
text('No response\n= 349', x = 0.7, y = 480)
text('Calm response\n= 448', x = 1.9, y = 480)
text('Agitated, charge\nor flee = 0', x = 3.1, y = 480)

range(d$date) # "2012-05-04" "2021-09-24"

#### restructuring ####
View(d)
motnp_long <- read_delim("data_processed/motnp_eles_long.csv", delim = ',')  # read in ALERT data for structure of how information should be.
View(motnp_long)
str(motnp_long)
#$ encounter      : num [1:4263] 1 2 3 5 6 8 9 10 12 14 ...
#$ date           : Date[1:4263], format: "2016-05-19" "2016-05-19" "2016-05-19" "2016-05-20" ...
#$ time           : num [1:4263] 0.421 0.485 0.491 0.342 0.59 ...
#$ gps_s          : num [1:4263] 1752912 1752031 1751973 1750466 1753007 ...
#$ gps_e          : num [1:4263] 2549597 2547860 2547883 2545378 2550000 ...
#$ total_elephants: chr [1:4263] "2" "9" "5" "9" ...
#$ total_id       : num [1:4263] 2 9 1 8 11 1 1 6 1 3 ...
#$ perc_id        : num [1:4263] 100 100 20 88.9 50 ...
#$ type           : chr [1:4263] "MO" "BH" "MO" "BH" ...
#$ number         : chr [1:4263] "e1" "e1" "e1" "e1" ...
#$ elephant       : chr [1:4263] "M1" "F9" "M8" "F1" ...

str(d)
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

# motnp_long$encountner = d$elephant_sighting_id
# motnp_long$date = d$date
# motnp_long$time = d$time
# motnp_long$gps_s = d$ DOESN'T EXIST
# motnp_long$gps_e = d$ DOESN'T EXIST
# motnp_long$total_elephants = d$ COUNT ROWS PER SIGHTING ID
# motnp_long$total_id = d$ COUNT ELEPHANT ID PER SIGHTING ID
# motnp_long$perc_id = d$ CALCULATE PERCENTAGE
# motnp_long$type = d$ sex_id IDENTIFY WHICH ARE FEMALES (~1-2%) AND MARK ALL OTHERS AS MO
# motnp_long$number = d$ IGNORE
# motnp_long$elephant = d$elephant_id

prop.table(table(d$sex_id))
#           1            2            3 
#0.9989531145 0.0006729978 0.0003738877   <-- 99.9% 1 = male, rest = BH/MX

d_long <- data_frame(encounter = d$elephant_sighting_id,
                     date = d$date,
                     time = d$time,
                     gps_s = NA,
                     gps_e = NA,
                     total_elephants = NA,
                     total_id = NA,
                     perc_id = NA,
                     type = ifelse(d$sex_id == 2, 'BH/MX', ifelse(d$sex_id == 3, 'BH/MX', 'MO')),
                     #number = NA,
                     elephant = ifelse(d$elephant_id == '-', NA, d$elephant_id))

which(d$sex_id == 2) ; which(d$sex_id == 3) ; which(d_long$type == 'BH/MX')

for(i in 1:length(d_long$total_elephants)) {
  d_long$total_elephants[i] <- length(which(d$elephant_sighting_id == d_long$encounter[i]))
  }

d_id <- d_long[!is.na(d_long$elephant),]
for(i in 1:length(d_long$total_id)) {
  d_long$total_id[i] <- length(which(d_id$encounter == d_long$encounter[i]))
}

d_long$perc_id <- 100*(d_long$total_id/d_long$total_elephants)

write_delim(d_long, 'data_processed/mpnp_eles_long_nogps.csv', delim = ',')

#### Nodes data frame ####
View(d)

motnp_nodes <- read_delim('data_processed/motnp_elenodes.csv', delim = ',')
colnames(motnp_nodes) # "id_no" "herd" "name" "age_class" "age_category" "family" "days" "comments" "sex" "num" "id" "count" "dem_class" -- carry over "id_no" "age_class" "age_category" "sex" "count" "dem_class"

ele_nodes <- data.frame(id_no = sort(unique(d_id$elephant)),
                        age_class = NA,
                        age_category = NA,
                        sex = NA,
                        count = NA,
                        dem_class = NA)
#for(i in 1:length(ele_nodes$id_no)){
#  ele_nodes$age_class[i] <- d$age_range_id[which(d$elephant_id == ele_nodes$id_no[i])]
#}

d$age_range_id[range(which(d$elephant_id == ele_nodes$id_no[1]))]                                   # 7 and 8, levels 2-10
d_age <- d[c(2,3,8)]
d_age$age_range_id <- as.numeric(d$age_range_id)
str(d_age)
d_age$age_range_id[range(which(d_age$elephant_id == ele_nodes$id_no[1]))]                           # 7 and 8, no other levels
ele_nodes$age_class[1] <- mean(d_age$age_range_id[which(d_age$elephant_id == ele_nodes$id_no[1])])  # mean = 7.055556
test <- d[d$elephant_id == 'B1000',] ; mean(as.numeric(test$age_range_id))                          # mean = 7.055556

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

ele_nodes$age_class[1] <- getmode(d_age$age_range_id[which(d_age$elephant_id == ele_nodes$id_no[1])])
test <- d[d$elephant_id == 'B1000',] ; getmode(as.numeric(test$age_range_id))

for(i in 1:length(ele_nodes$age_class)) {
  ele_nodes$age_class[i] <- getmode(d_age$age_range_id[which(d_age$elephant_id == ele_nodes$id_no[i])])
  ele_nodes$sex[i] <- d$sex_id[which(d$elephant_id == ele_nodes$id_no[i])]
  ele_nodes$count[i] <- length(which(d$elephant_id == ele_nodes$id_no[i]))
}

ele_nodes$sex <- ifelse(ele_nodes$sex == 2, 'F', ifelse(ele_nodes$sex == 3, 'F', 'M'))

ele_nodes$age_category <- ifelse(ele_nodes$age_class < 3, 'Calf',
                                 ifelse(ele_nodes$age_class == 3, 'Juvenile',
                                        ifelse(ele_nodes$age_class > 3 & ele_nodes$age_class <= 5, 'Pubescent',
                                               ifelse(ele_nodes$age_class > 5 & ele_nodes$age_class <= 8, 'Adult', 'UNK'))))
ele_nodes$age_class <- as.factor(ele_nodes$age_class)

ele_nodes$age <- ifelse(ele_nodes$age_category == 'Adult', 'A',
                        ifelse(ele_nodes$age_category == 'Pubescent', 'P',
                               ifelse(ele_nodes$age_category == 'Juvenile', 'J', 'C')))
ele_nodes$dem_class <- paste(ele_nodes$age, ele_nodes$sex, sep = '')
ele_nodes <- ele_nodes[c(1:6)]

write_delim(ele_nodes, 'data_processed/mpnp_elenodes.csv')

