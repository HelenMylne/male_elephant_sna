# EfA data model prep
#### Information ####
# Data collected by Elephants for Africa (EfA) 2012-2021
# Data supplied by Dr Kate Evans
#### Set up ####
library(tidyverse)
library(rstan)
library(rethinking)
library(cmdstanr)
library(lubridate)

#### read in processed data files ####
aa <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1.250_22.03.10.csv') %>% distinct()
ab <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings251.387_22.03.08.csv') %>% distinct()
ac <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings388.400_22.03.08.csv') %>% distinct()
ad <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings401.450_22.03.08.csv') %>% distinct()
ae <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings451.500_22.03.08.csv') %>% distinct()
af <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings501.550_22.03.08.csv') %>% distinct()
ag <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings551.600_22.03.08.csv') %>% distinct()
ah <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings601.650_22.03.08.csv') %>% distinct()
ai <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings651.700_22.03.08.csv') %>% distinct()
aj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings701.750_22.03.08.csv') %>% distinct()
ak <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings751.800_22.03.08.csv') %>% distinct()
al <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings801.850_22.03.08.csv') %>% distinct()
am <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings851.900_22.03.08.csv') %>% distinct()
an <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings901.950_22.03.08.csv') %>% distinct()
ao <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings951.1000_22.03.08.csv') %>% distinct()
ap <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1001.1050_22.03.08.csv') %>% distinct()
aq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1051.1100_22.03.08.csv') %>% distinct()
ar <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1101.1150_22.03.08.csv') %>% distinct()
as <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1151.1200_22.03.08.csv') %>% distinct()
at <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1201.1250_22.03.08.csv') %>% distinct()
au <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1251.1300_22.03.08.csv') %>% distinct()
av <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1301.1350_22.03.08.csv') %>% distinct()
aw <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1351.1400_22.03.08.csv') %>% distinct()
ax <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1401.1450_22.03.08.csv') %>% distinct()
ay <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1451.1500_22.03.08.csv') %>% distinct()
az <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1501.1550_22.03.08.csv') %>% distinct()

ba <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1551.1600_22.03.08.csv') %>% distinct()
bb <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1601.1650_22.03.08.csv') %>% distinct()
bc <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1651.1700_22.03.08.csv') %>% distinct()
bd <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1701.1750_22.03.08.csv') %>% distinct()
be <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1751.1800_22.03.08.csv') %>% distinct()
bf <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1801.1850_22.03.08.csv') %>% distinct()
bg <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1851.1900_22.03.08.csv') %>% distinct()
bh <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1901.1950_22.03.08.csv') %>% distinct()
bi <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings1951.2000_22.03.08.csv') %>% distinct()
bj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2001.2050_22.03.08.csv') %>% distinct()
bk <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2051.2100_22.03.08.csv') %>% distinct()
bl <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2101.2150_22.03.08.csv') %>% distinct()
bm <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2151.2200_22.03.08.csv') %>% distinct()
bn <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2201.2250_22.03.08.csv') %>% distinct()
bo <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2251.2300_22.03.08.csv') %>% distinct()
bp <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2301.2350_22.03.08.csv') %>% distinct()
bq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2351.2400_22.03.08.csv') %>% distinct()
br <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2401.2450_22.03.08.csv') %>% distinct()
bs <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2451.2500_22.03.08.csv') %>% distinct()
bt <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2501.2550_22.03.08.csv') %>% distinct()
bu <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2551.2600_22.03.08.csv') %>% distinct()
bv <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2601.2650_22.03.08.csv') %>% distinct()
bw <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2651.2700_22.03.08.csv') %>% distinct()
bx <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2701.2750_22.03.08.csv') %>% distinct()
by <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2751.2800_22.03.08.csv') %>% distinct()
bz <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2801.2850_22.03.08.csv') %>% distinct()

ca <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2851.2900_22.03.08.csv') %>% distinct()
cb <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2901.2950_22.03.08.csv') %>% distinct()
cc <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings2951.3000_22.03.08.csv') %>% distinct()
cd <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3001.3050_22.03.08.csv') %>% distinct()
ce <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3051.3100_22.03.08.csv') %>% distinct()
cf <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3101.3150_22.03.08.csv') %>% distinct()
cg <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3151.3200_22.03.08.csv') %>% distinct()
ch <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3201.3250_22.03.08.csv') %>% distinct()
ci <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3251.3300_22.03.08.csv') %>% distinct()
cj <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3301.3350_22.03.08.csv') %>% distinct()
ck <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3351.3400_22.03.08.csv') %>% distinct()
cl <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3401.3450_22.03.08.csv') %>% distinct()
cm <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3451.3500_22.03.08.csv') %>% distinct()
cn <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3501.3550_22.03.08.csv') %>% distinct()
co <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3551.3600_22.03.08.csv') %>% distinct()
cp <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3601.3650_22.03.08.csv') %>% distinct()
cq <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3651.3700_22.03.08.csv') %>% distinct()
cr <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3701.3750_22.03.08.csv') %>% distinct()
cs <- read_csv('data_processed/mpnp_pairwiseevents/mpnp_bayesian_allpairwiseevents_sightings3751.3765_22.03.08.csv') %>% distinct()

# merge
all <- rbind(aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,ay,az,
             ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk,bl,bm,bn,bo,bp,bq,br,bs,bt,bu,bv,bw,bx,by,bz,
             ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl,cm,cn,co,cp,cq,cr,cs)
length(unique(all$obs_id))          # 3765

# clean environment
rm(aa,ab,ac,ad,ae,af,ag,ah,ai,aj,ak,al,am,an,ao,ap,aq,ar,as,at,au,av,aw,ax,ay,az,ba,bb,bc,bd,be,bf,bg,bh,bi,bj,bk,bl,bm,bn,bo,bp,bq,br,bs,bt,bu,bv,bw,bx,by,bz,ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl,cm,cn,co,cp,cq,cr,cs)

#### Recreate how data were produced to match up obs_id to actual encounters ####
# sightings data
s <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211214.xlsx')
str(s)
colnames(s)[c(1:23,57)] <- s[2,c(1:23,57)]
colnames(s)[24:56] <- c('CM','CF','CU','CM','CF','CU','JM','JF','JU','YPM','YPF','YPU','OPM','OPF','OPU',
                        'YAM','YAF','YAU','MAM','MAF','MAU','OAM','OAF','OAU','UM','UF','UU','SM','SF','SU',
                        'AM','AF','AU')
s <- s[3:nrow(s),]
s <- janitor::clean_names(s)

# individual data
efa <- readxl::read_excel('data_raw/Raw_EfA_ElephantVisuals_IndividualsGroups_Evans211019.xlsx')
str(efa)
efa$time_cat <- lubridate::hour(efa$Time)
efa <- separate(efa, Time, into = c('wrong_date','time'), sep = ' ')
efa$time <- hms::as_hms(efa$time)
efa$Date <- lubridate::as_date(efa$Date)
efa <- efa[,c(1:5,7,33,8:15,18:27,31:32)]
efa <- janitor::clean_names(efa)

efa_long <- data.frame(encounter = efa$elephant_sighting_id,
                       date = efa$date,
                       time = efa$time,
                       gps_s = NA,
                       gps_e = NA,
                       total_elephants = NA,
                       total_id = NA,
                       perc_id = NA,
                       type = ifelse(efa$sex_id == 1, 'MO', 'BH/MX'),
                       elephant = ifelse(efa$elephant_id == '-', NA, efa$elephant_id),
                       sex = ifelse(efa$sex_id == 1, 'M','F/U'),
                       age_range = efa$age_range_id)

for(i in 1:length(efa_long$total_elephants)) {
  efa_long$total_elephants[i] <- length(which(efa$elephant_sighting_id == efa_long$encounter[i]))
}

efa_id <- efa_long[!is.na(efa_long$elephant),]
for(i in 1:length(efa_id$total_id)) {
  efa_id$total_id[i] <- length(which(efa_id$encounter == efa_long$encounter[i]))
}

efa_id$perc_id <- 100*(efa_id$total_id/efa_id$total_elephants)

colnames(s)[1] <- 'encounter'
s$encounter <- as.numeric(s$encounter)
efa_gps <- left_join(x = efa_id, y = s, by = 'encounter') # first few rows are blank because 2 datasets formed at slightly different times -- group sightings only goes as far as 29th January 2021, individual sightings go until 24th September 2021
efa_gps <- efa_gps[,c(1,13,10,2,3,18,19,6:9,11,12)]
colnames(efa_gps)[4] <- c('date')

length(unique(efa_long$encounter)) # 5162
length(unique(efa_id$encounter))   # 3736
length(unique(efa_gps$encounter))  # 3736

efa_gps$location <- paste(efa_gps$latitude, efa_gps$longitude, sep = '_')
efa_gps$longitude[which(efa_gps$longitude < 20)] <- NA
efa_gps$longitude[which(efa_gps$longitude < 21)] # 20.512689999999
efa_gps$longitude[which(efa_gps$latitude < -20)] # 24.750029999999999 24.707709999999999 24.707709999999999 24.768260000000001 24.764379999999999 24.765647999999999 24.528880000000001 24.528880000000001

### create group-by-individual matrix
eles_asnipe <- efa_gps[,c(1,4,5,3,14)]  # encounter, date, time, elephant, location
eles_asnipe$Date <- as.integer(eles_asnipe$date)
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)         # start from 1, not 1st January 1970
eles_asnipe$Time <- hour(eles_asnipe$time)*60*60 + minute(eles_asnipe$time)*60 + second(eles_asnipe$time) # convert time values to seconds through day
eles_asnipe <- eles_asnipe[,c(1,6,7,4,5)]
colnames(eles_asnipe)[2:5] <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
eles_asnipe$group <- with(eles_asnipe, paste(Date, Time, Location, sep = '_'))

length(unique(eles_asnipe$encounter)) # 3736
length(unique(eles_asnipe$group))     # 3765
eles_asnipe$d_pad <- str_pad(eles_asnipe$Date, 3, pad = '0')
eles_asnipe$group <- paste(eles_asnipe$d_pad, eles_asnipe$Time, eles_asnipe$Location, sep = '_')
eles_asnipe$group_id <- as.integer(as.factor(eles_asnipe$group))
length(unique(eles_asnipe$encounter)) # 3736
max(eles_asnipe$group_id)             # 3765
eles_asnipe2 <- eles_asnipe[,c(4,8)]
colnames(eles_asnipe2)[2] <- 'group'
eles_asnipe2 <- data.table::setDT(eles_asnipe2)
gbi_matrix <- spatsoc::get_gbi(DT = eles_asnipe2, group = 'group', id = 'ID') # 3765 rows --  ROW 1 IN gbi_matrix = obs_id 1 IN FINAL OUTPUT

obs_id <- rownames(gbi_matrix)
asnipe <- unique(sort(as.integer(as.factor(eles_asnipe$group))))
which(obs_id != asnipe) # none -- both 1:3765

## clean up
rm(efa, efa_gps, efa_id, efa_long, eles_asnipe2, gbi_matrix, s, asnipe, i, obs_id)

#### create new variable to match up obs_id to groups ####
### SO: Create a new variable using eles_asnipe group -- do a for loop that adds 1 to the value every time a group changes and retains it whenever it is the same -- this will create a variable which is unique for 3765 encounters but which goes in the same order as the obs_id
### use observations dataframe which is unique by GROUP (NOT ENCOUNTER ID NUMBER) and then add number in order of group
eles_asnipe3 <- eles_asnipe %>% arrange(desc(encounter)) %>% arrange(desc(Time)) %>% arrange(desc(Date)) # convert to an order where group is genuinely going in order
which(eles_asnipe$encounter != eles_asnipe3$encounter)  # check to see which ones have changed -- row 9 = day 3421 = now in descending order as with date
eles_asnipe3$in_order_date <- c('yes',rep(NA,nrow(eles_asnipe3)-1))
eles_asnipe3$in_order_time <- c('yes',rep(NA,nrow(eles_asnipe3)-1))
for(i in 2:nrow(eles_asnipe3)){
  eles_asnipe3$in_order_date[i] <- ifelse(eles_asnipe3$Date[i] > eles_asnipe3$Date[i-1], 'no', 'yes')
  eles_asnipe3$in_order_time[i] <- ifelse(eles_asnipe3$Date[i] == eles_asnipe3$Date[i-1], 
                                         ifelse(eles_asnipe3$Time[i] <= eles_asnipe3$Time[i-1], 'yes', 'no'),
                                         'yes')
}
table(eles_asnipe3$in_order_date) ; table(eles_asnipe3$in_order_time) # sightings now in reverse chronological order, both by date and time

(N <- length(unique(eles_asnipe3$group)))   # 3765
length(unique(eles_asnipe3$encounter)) - N  # -29 -- encounters were not reliably numbered throughout dataset
eles_asnipe3$group_id_new <- c(N, rep(NA, nrow(eles_asnipe3)-1))
for(i in 2:nrow(eles_asnipe3)){
  eles_asnipe3$group_id_new[i] <- ifelse(eles_asnipe3$group[i] == eles_asnipe3$group[i-1],
                                         eles_asnipe3$group_id_new[i-1],
                                         eles_asnipe3$group_id_new[i-1]-1)
}
summary(eles_asnipe3$group_id_new)

observations <- eles_asnipe3[,c('Date','Time','Location','group_id_new')] %>% distinct() # 3765 observations

## SO... to identify which sightings in main data frame are actually represented in binomial data (all): merge 'observations' data with Binomial social interactions
rm(eles_asnipe, eles_asnipe3, N)
#### merge sightings information into dyad data #####
colnames(observations) <- c('date', 'time', 'location', 'obs_id')
all <- left_join(x = all, y = observations, by = 'obs_id')
head(all,10)
tail(all,10)
rm(observations)

#### convert to Bernoulli model data format -- can't actually use Bernoulli as would require too much computing power ####
all$dyad <- paste(all$node_1, all$node_2, sep = '_')
all$dyad_id <- as.integer(as.factor(all$dyad))
all$location_id <- as.integer(as.factor(all$location))
head(all)

max(all$date) - min(all$date) # 3430 days. Time period for ALERT data = 581 days --> split into 6 time windows (7 boundaries)
periods <- seq(from = min(all$date), to = max(all$date), length.out = 7)
events <- all[all$social_event == 1,] ; events <- events[!is.na(events$node_1),] # cuts down to 55490933 dyad pairs from 13043 -- assigning to time window doesn't take forever
events$period <- NA
for(i in 1:nrow(events)){
  events$period[i] <- which(periods <= events$date[i])[length(which(periods <= events$date[i]))] # take last value in vector
  if(i %% 1000 == 0) {print(i)}
}
range(events$obs_id)

periods ; View(events[c(sample(x = 1:nrow(events), size = 20, replace = F)),]) # visual check that periods have come out right

# check elephants all match up
length(unique(all$node_1))
length(unique(all$node_2))
length(unique(events$node_1))
length(unique(events$node_2))
#rm(all)

#### convert to Binomial model data format -- aggregate all sightings of each dyad together into a count ####
df_split <- events %>%
  group_by(node_1, node_2, period) %>%
  summarise(event_count=sum(social_event),
            dyad_id=cur_group_id())
head(df_split)
rm(events)

### add ID numbers
eles <- read_csv('data_processed/mpnp_eles_long_22.03.08.csv') %>% 
  select(elephant) %>% 
  distinct()
eles$node_1 <- as.integer(as.factor(eles$elephant))
colnames(eles)[1] <- c('id_1')
df <- left_join(df_split, eles, by = 'node_1') %>% distinct()
colnames(eles) <- c('id_2','node_2')
df <- left_join(df, eles, by = 'node_2') %>% distinct()
head(df) ; tail(df)
rm(df_split)

df$dyad <- paste(df$id_1, df$id_2, sep = '_')
df$dyad_id_period <- df$dyad_id              # every dyad has it's own ID number, including if same dyad in a different time window
df$dyad_id <- as.integer(as.factor(df$dyad)) # every dyad has it's own ID number, but same dyad in different windows share ID number

#### create dyad row for all pairs per period ####
dyads <- data.frame(id_1 = rep(sort(eles$id_2), each = nrow(eles)),
                    id_2 = rep(sort(eles$id_2), nrow(eles)))

colnames(eles) <- c('id_1','node_1')
dyads <- left_join(dyads, eles, by = 'id_1')
colnames(eles) <- c('id_2','node_2')
dyads <- left_join(dyads, eles, by = 'id_2')
dyads <- dyads[dyads$node_1 < dyads$node_2,]
rm(all, eles, i, periods)

dyads <- data.frame(id_1 = rep(dyads$id_1, length(unique(df$period))),
                    id_2 = rep(dyads$id_2, length(unique(df$period))),
                    node_1 = rep(dyads$node_1, length(unique(df$period))),
                    node_2 = rep(dyads$node_2, length(unique(df$period))),
                    period = rep(sort(unique(df$period)), each = nrow(dyads)))

head(df) ; head(dyads)
data <- left_join(x = dyads, y = df, by = c('id_1','id_2','period'))
data <- data[,c(1:5,8:10)]
colnames(data)[3:4] <- c('node_1','node_2')
data$event_count <- ifelse(is.na(data$event_count) == TRUE, 0, data$event_count)
table(data$event_count)
data$dyad <- paste(data$id_1, data$id_2, sep = '_')
data$dyad_id <- as.integer(as.factor(data$dyad))
head(data, 20)

periods <- data.frame(period_start = seq(from = min(sightings$obs_date), to = max(sightings$obs_date), length.out = 32)[1:31],
                      period = 1:31)
data <- left_join(x = data, y = periods, by = 'period')
head(data)

#### add data about nodes ####
eles <- read_csv('data_processed/mpnp_eles_long_22.03.08.csv') %>% 
  select(elephant, sex, age_range) %>% 
  distinct()
nrow(eles) - length(unique(eles$elephant)) # 472 elephants age/sex reclassified
unique(eles$sex)
eles <- eles[,c(1,3)] %>% distinct()
nrow(eles) - length(unique(eles$elephant)) # 450 elephants age reclassified

table(eles$age_range) # no 1, 9 is meant to be UK and 10 is not classified, but no 9s and 148 10s I'm thinking 10 is unknown age
eles$age_range_NA <- ifelse(eles$age_range == 10, NA, eles$age_range)

eles$age_unsure <- NA ; eles$age_min <- NA ; eles$age_max <- NA ; eles$age_maxmin <- NA ; eles$age_median <- NA
for (i in 1:nrow(eles)) {
  individual <- eles[eles$elephant == eles$elephant[i],]
  summary <- summary(individual$age_range_NA, na.rm = T)
  eles$age_unsure[i] <- nrow(individual)
  eles$age_min[i] <- summary[1]
  eles$age_max[i] <- summary[6]
  eles$age_maxmin[i] <- eles$age_max[i] - eles$age_min[i]
  eles$age_median[i] <- summary[3]
}
summary(eles$age_median)

colnames(eles)[1] <- 'id'
mpnp <- left_join(x = data, y = eles, by = 'id')

#### write to file ####
readr::write_delim(mpnp, 'data_processed/mpnp_bayesian_pairwiseevents_22.04.14.csv', delim = ',')

## clean environment
rm(list = ls())
