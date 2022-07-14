#### Information ####
# Script to process and plot association data from Makgadikgadi Pans National Park, Bostwana
# Data collected: 1st November 2019 to 5th August 2021
# Collected by: 
# Data supplied by: Dr Kate Evans, Elephants for Africa, 2nd September 2021

#### Set up ####
library(tidyverse)  # data manipulation
library(lubridate)  # sort date columns out
library(zoo)        # sort date columns out
library(asnipe)     # generating networks
library(rethinking) # col.alpha
library(igraph)     # plotting networks

#### Data import ####
d <- readxl::read_excel("data_raw/Raw_EfA_ElephantDatabase_2019.2021_Evans210902.xlsx", sheet = 1)
str(d)
length(unique((d$ID_Elephant_Visual)))   # 1018 -- all unique = unique number per elephant sighting
length(unique((d$Elephant_Sighting_ID))) # 367 -- unique number per group sighting
length(unique((d$Elephant_ID)))          # 608 identified individuals
table(d$Elephant_ID)                     # 409 unidentified, rest all only seen once. T=temporaryID, B=checked.
table(d$In_The_Field_Elephant_Code)      # Unique code per elephant per day, max=BR=73 elephants seen in one day
length(unique((d$Date)))                 # 68 days
str(d$Time)                              # datetime rather than just date
d$time_cat <- lubridate::hour(d$Time)
e <- separate(d, Time, into = c('wrong_date','time'), sep = ' ')
head(e)
e$time2 <- hms::as_hms(e$time)
d <- e[,c(1:5,7,38,37,8:36)]
str(d$time2)
d$Age_Range_ID <- as.factor(d$Age_Range_ID)
summary(d$Age_Range_ID) # 1=<1y, 2=1-4y, 3=5-9y, 4=10-15y, 5=16-20y, 6=21-25y, 7=26-35y, 8=>35y, 9=UK
summary(d$Age)          # all NA -- 4=10-15y, 5=16-20y, 6=21-25y, 7=26-35y, 8=>35y
d$Activity_ID <- as.factor(d$Activity_ID)   
summary(d$Activity_ID)      # 2-8, not integer, 6 most common, no ethogram in doc from Kate
table(d$In_musth)           # 7 yes, 955 no
summary(d$Reaction_Index)   # -4=run away, 1=don't notice, 2=unbothered, 3=agitated, 4=charge, range 1-3
summary(d$Number_In_The_Group)   # 1-26
d$Distance_To_Observer <- as.numeric(d$Distance_To_Observer)
summary(d$Distance_To_Observer)  # 1-600, use a range finder (GPS elephants = range+GPS_car+bearing)
summary(d$Habitat_ID)       # all NA
table(d$Vegetal_Species_ID) # only 3 filled
summary(d$Feeding_Type_ID)  # all NA
summary(d$Pictures_Taken)   # 404 no pictures, 613 with pictures, 1 NA
table(d$Sick)               # 1 sick, 963 healthy, 53 unknown
d$Physical_Condition_ID <- as.factor(d$Physical_Condition_ID)
summary(d$Physical_Condition_ID) # 1=emaciated, 2=very thin, 3=thin, 4=good, 5=fat

# L/R_TGS_**_ID -- TGS = temporal gland secretion, L/R = Left/Right
table(d$L_TGS_Length_ID)    # vast majority are 0 or unknown. 0=340, 1=16, 2=18, 3=15, 4=16
table(d$L_TGS_Width_ID)     # t,m,w=thin,medium,wide?? t=16,m=46,w=3, 0=343
table(d$L_TGS_Age_ID)       # n=new=60, o=old=8, 0=336

table(d$R_TGS_Length_ID)    # 0=383, 1=19, 2=25, 3=22, 4=15
table(d$R_TGS_Width_ID)     # t=16, m=63, w=3, 0=386
table(d$R_TGS_Age_ID)       # n=70, o=12, 0=382

table(d$TG_Swelling_ID)     # 0=919 no swelling of temporal glands, 1=5 swelled glands, 93=unknown
summary(d$Hind_Foot_Length_Average)          # all NA
summary(d$Front_Foot_Circumference_Average)  # all NA
summary(d$Height)           # all NA
table(d$Sex_ID)             # 618 x 1, 400 x NA --> not sure if this indicates 618 females or males. EfA data states to record gender if they are female, which would suggest that NA=M/U and 1=F.
summary(d$...32) ; summary(d$...33) ; summary(d$...34) ; summary(d$...35) # all NA

str(d)
d <- d[,c(1:5,7:9,11:15,19:28,32:33)] # NOTE -- ANY FURTHER DATA FROM KATE MAY BE ABLE TO USE SOME OF THE DROPPED COLUMNS
head(d)
d$date <- lubridate::as_date(d$date)
d %>% rename(time = time2)
d <- janitor::clean_names(d)

#### Clean data import ####
d <- readxl::read_excel("data_raw/Raw_ElephantDatabase_Evans210902.xlsx", sheet = 1)
d$Date <- lubridate::as_date(d$Date)
d$time_cat <- lubridate::hour(d$Time)
d <- separate(d, Time, into = c('wrong_date','time'), sep = ' ')
d$time <- hms::as_hms(d$time)
d <- d[,c(1:5,7,37,8:36)]
d$Age_Range_ID <- as.factor(d$Age_Range_ID)
d$Activity_ID <- as.factor(d$Activity_ID)   
d$Distance_To_Observer <- as.numeric(d$Distance_To_Observer)
d$Physical_Condition_ID <- as.factor(d$Physical_Condition_ID)
d <- d[,c(1:8,10:14,18:27,31:32)] # ANY FURTHER DATA FROM KATE MAY BE ABLE TO USE SOME OF THE DROPPED COLUMNS
d <- janitor::clean_names(d)

#### taking a look ####
table(d$Elephant_ID)
hist(table(d$elephant_sighting_id))
hist(table(d$elephant_sighting_id), breaks = 30, las = 1)
length(which(table(d$elephant_sighting_id) == 1 | table(d$elephant_sighting_id) == 2)) # 237
length(which(table(d$elephant_sighting_id) == 1)) # 166

summary(d$reaction_index)
d$reaction_index <- factor(d$reaction_index, levels = c(1:4))
barplot(table(d$reaction_index), xlab = 'reaction index score', ylab = 'count',
        las = 1, ylim = c(0,160), xlim = c(0,5))
abline(h = 0, lwd = 2) ; axis(2, at = seq(0,150,10), las = 1)
#abline(v = seq(0,5,0.5)) ; abline(v = seq(0,5,0.1), lty = 3)
text('No response\n= 141', x = 0.7, y = 90)
text('Calm response\n= 147', x = 1.9, y = 90)
text('Agitated response\n= 1', x = 3.1, y = 90)
text('Charge/flee\n= 0', x = 4.3, y = 90)






