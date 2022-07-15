## First attempt to roughly build a Bayesian model of the ALERT data
##### set up and data import ####
library(tidyverse)
library(rstan)
library(asnipe)
library(dplyr)
library(igraph)
## sightings
sightings <- read_delim(file = 'data_processed/motnp_encounters_correctid.csv', delim = ',')
sightings <- s[1:22]
str(sightings)

## elephant encounters
eles <- read_delim(file = 'data_processed/motnp_eles_long_correctid.csv', delim = ',')
colnames(eles)
eles <- eles[c(1:12,14)]
str(eles)

## nodes
nodes <- read_delim(file = 'data_processed/motnp_elenodes.csv', delim = ',')
colnames(nodes)
str(nodes)

## links
links <- read_delim('data_processed/motnp_elelinks.csv', delim = ',')

##### plan for analysis #####
## ideas for how to "Bayesianise" SRI
# SRI ~ Binomial(p_ij, N_ij)
# logit(p_ij) = beta_s[sex_ij] * beta_a[age_ij] * beta_f[family_ij] ;
# N ~ Normal(mu, sigma)
# mu = ya + yb + yab + x
# sigma = constant or Exponential(1)
# (beta_s, beta_a) ~ 0,1 bounded Normal distribution
# beta_f ~ Exponential(1)
# Where sex_ij: 1=MM, 2=MF, 3=FF, 4=MU, 5=FU, 6=UU;  age_ij: 1=AA, 2=AP, 3=AJ, 4=AC, 5=PP, 6=PJ...10=CC; family_ij: 1=mother-calf, 2=same breeding herd, 3=different breeding herds). Could also form a composite of sex and age dyads, so AM-CU is definitely different to AF-CU.
# also good to incorporate uncertainty over misidentification or failure to identify/count individuals.

##### Get data into aggregated dyad form #####
# columns: encounter number, node_id1, node_id2, type_1, type_2, event_number, location_id
# aggregate to: node_id1, node_id2, count_interactions, dyad_id, node_id1_factor, node_id2_factor
colnames(e)
# encounter = event, gps_s&gps_e = location, node_id1&node_id2 = elephant, 
colnames(n)
# type_1&type_2 = dem_class

test <- e[e$encounter < 11,]
View(test)
# sighting 1 = M1 and M11 --> event 1 = M1-M11
# sighting 2 = F9, U8, F8, F13, M26, M87, F19, M19, M24 --> events 2-32(?) = F9-U8, F9-F8, F9-F13...F19-M19, F19-M24, M19-M24
# sighting 3 = M8 and 4 unknown others --> events 33-43? or no event?
# sighting 4 = 17 MX -- no elephants identified
# sighting 5 = F1, F2, F3, U12, F6, F7, M6, M7, U5 --> events 44-74, total ID recorded as 8 but 9 IDs given?
# sighting 6 = F8, F13, F21, M15, M26, M87, M3, M5, M8, M9, M52 + 11 unknown others --> what to do about all of the missing IDs?
# sighting 7 = 5 MO -- no elephants identified
# sighting 8 = F40 + 9 others --> missing IDs?
# sighting 9 = M3 alone --> no dyad so can't be an event, but also can't just drop it as it is a sighting?
# sighting 10 = M4, M5, M13, M18, M113, M118 + 2 others --> missing IDs? 

test2 <- e[e$encounter < 3,]
test2
## CURRENT FORMAT:
#    encounter       date  time   gps_s   gps_e total_elephants total_id_dy total_id_hkm total_id_diff perc_id_dy perc_id_hkm type  elephant
# 1          1 2016-05-19 0.421 1752912 2549597               2           2            2             0        100         100   MO        M1      
# 2          2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH        F9      
# 3          1 2016-05-19 0.421 1752912 2549597               2           2            2             0        100         100   MO       M11     
# 4          2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH        U8      
# 5          2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH        F8      
# 6          2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH       F13     
# 7          2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH       M26     
# 8          2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH       M87     
# 9          2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH       F19     
# 10         2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH       M19     
# 11         2 2016-05-19 0.485 1752031 2547860               9           9            9             0        100         100   BH       M24     

## TARGET FORMAT:
#   node_1 node_2   type_1   type_2 event location
# 1    Rey   Leia Lifeform Lifeform     1        E
# 2    Rey   Leia Lifeform Lifeform     1        C
# 3    Rey   Leia Lifeform Lifeform     1        D
# 4    Rey   Leia Lifeform Lifeform     1        D
# 5    Rey   Leia Lifeform Lifeform     1        C
# 6    Rey   Leia Lifeform Lifeform     1        F

test2$location <- paste(test2$gps_s, test2$gps_e, sep = '_')
test2 <- test2[c(1:3,14,6:13)]
head(test2)

el <- read_delim('data_processed/motnp_elelinks.csv', delim = ',')
head(el)
test3 <- el[el$id1 == 'M1' | el$id1 == 'F9' | el$id1 == 'M11' | el$id1 == 'U8' | el$id1 == 'F8' | el$id1 == 'F13' | 
              el$id1 == 'M26' | el$id1 == 'M87' | el$id1 == 'F19' | el$id1 == 'M19' | el$id1 == 'M24',]
test3 <- test3[test3$id2 == 'M1' | test3$id2 == 'F9' | test3$id2 == 'M11' | test3$id2 == 'U8' | test3$id2 == 'F8' | test3$id2 == 'F13' | 
                 test3$id2 == 'M26' | test3$id2 == 'M87' | test3$id2 == 'F19' | test3$id2 == 'M19' | test3$id2 == 'M24',]
head(test3)
#   from  to    type  weight sex1  sex2  id1   id2   dyad   ageclass1 ageclass2 age1  age2  days1 days2 count1 count2 age.dyad dem_class1 dem_class2
# 1 F0008 F0009 FF    0.0213 F     F     F8    F9    F8_F9  Adult     Adult     20-25 20-25     9    10     11     11 AA       AF         AF        
# 2 F0009 M0011 MF    0.136  F     M     F9    M11   F9_M11 Adult     Adult     20-25 40+      10    16     11     18 AA       AF         AM        
# 3 F0013 F0009 FF    0.0247 F     F     F13   F9    F13_F9 Pubescent Adult     15-19 20-25    10    10     11     11 AP       PF         AF        
# 4 F0009 M0001 MF    0.116  F     M     F9    M1    F9_M1  Adult     Adult     20-25 25-40    10    29     11     44 AA       AF         AM        
# 5 F0019 F0009 FF    0.0911 F     F     F19   F9    F19_F9 Adult     Adult     35-50 20-25    16    10     21     11 AA       AF         AF        
# 6 F0009 M0019 MF    0.338  F     M     F9    M19   F9_M19 Adult     Pubescent 20-25 10-15    10    13     11     13 AP       AF         PM

# 55 dyads formed of 11 elephants from 2 events.
dyads <- data.frame(node_1 = test3$id1, node_2 = test3$id2,
                    type_1 = test3$dem_class1, type_2 = test3$dem_class2)
events <- data.frame(event = test2$encounter, location = test2$location, elephant = test2$elephant)
## DYADS                           |   ## EVENTS
#    node_1 node_2 type_1 type_2   |   #   event        location elephant
# 1  F8     F9     AF     AF       |   # 1     1 1752912_2549597       M1
# 2  F9     M11    AF     AM       |   # 2     2 1752031_2547860       F9
# 3  F13    F9     PF     AF       |   # 3     1 1752912_2549597      M11
# 4  F9     M1     AF     AM       |   # 4     2 1752031_2547860       U8
# 5  F19    F9     AF     AF       |   # 5     2 1752031_2547860       F8
# 6  F9     M19    AF     PM       |   # 6     2 1752031_2547860      F13

## TARGET FORMAT:
#   node_1 node_2   type_1   type_2 event location
# 1    Rey   Leia Lifeform Lifeform     1        E
# 2    Rey   Leia Lifeform Lifeform     1        C
# 3    Rey   Leia Lifeform Lifeform     1        D
# 4    Rey   Leia Lifeform Lifeform     1        D
# 5    Rey   Leia Lifeform Lifeform     1        C
# 6    Rey   Leia Lifeform Lifeform     1        F

test4 <- e[e$perc_id_hkm == 100,]
head(test4)
test4$location <- paste(test4$gps_s,test4$gps_e,sep='_')
eles100 <- as.data.frame(unique(test4$elephant))
colnames(eles100) <- 'from'
test5 <- right_join(el, eles100, by = 'from')
colnames(eles100) <- 'to'
test5 <- right_join(test5, eles100, by = 'to')
length(unique(test5$from))

dyads100 <- data.frame(node_1 = test5$id1, node_2 = test5$id2,
                       type_1 = test5$dem_class1, type_2 = test5$dem_class2)
events100 <- data.frame(event = test4$encounter, location = test4$location, elephant = test4$elephant)
head(dyads100)
head(events100)
which(events100$elephant == 'M11')
events100 <- rename(events100, node_1 = elephant)
test100 <- full_join(dyads100, events100, by = 'node_1')

View(test100)
length(unique(test100$node_1)) # 469 -- should only be 263?


e$location <- paste(e$gps_s, e$gps_e, sep = '_')
e <- e[c(1:3,14,6:13)]

dyads_all <- data.frame(node_1 = el$id1, node_2 = el$id2,
                        type_1 = el$dem_class1, type_2 = el$dem_class2)
events_all <- data.frame(event = e$encounter, location = e$location, node_1 = e$elephant)

test_all <- full_join(dyads_all, events_all, by = 'node_1')
f1 <- test_all[test_all$node_1 == 'F1',]
table(f1$node_2)
test_all <- left_join(x = dyads_all, y = events_all, by = 'node_1')
f1 <- test_all[test_all$node_1 == 'F1',]
table(f1$node_2)
test_all <- left_join(y = dyads_all, x = events_all, by = 'node_1')
f1 <- test_all[test_all$node_1 == 'F1',]
table(f1$node_2)
test_all <- inner_join(x = dyads_all, y = events_all, by = 'node_1')
f1 <- test_all[test_all$node_1 == 'F1',]
table(f1$node_2)


(count_sightings_individuals <- as.data.frame(table(events_all$node_1)))
summary(count_sightings_individuals)
length(count_sightings_individuals$Var1) # 472
# join is producing a new line for every dyad whenever there is a sighting of node 1

head(test_all, 50)

test_all_noduplicates <- test_all %>% 
  rowwise() %>% 
  mutate(dyad = paste(sort(c(node_1,node_2)), collapse = '_')) %>% 
  ungroup() %>% 
  distinct_at(vars(dyad))
length(test_all_noduplicates$dyad) # 108352

for(i in 1:nrow(test_all)){
  test_all$dyad[i] <- paste(sort(c(test_all$node_1[i], test_all$node_2[i])), collapse = '_')
}
length(unique(test_all$dyad))      # 108352
which(sort(unique(test_all$dyad)) != sort(unique(test_all_noduplicates$dyad))) # all match

test_all$dyad_fct <- as.factor(as.numeric(as.factor(test_all$dyad)))

#test_all$id1_fct <- as.factor(as.numeric(as.factor(test_all$node_1)))
#test_all$id2_fct <- as.factor(as.numeric(as.factor(test_all$node_2)))
head(test_all, 50) # NOPE - THIS MEANS I HAVE DIFFERENT NUMBERS FOR THE SAME ELEPHANT DEPENDING ON THE NODE
# test_all <- test_all[,1:8]

eles_all <- data.frame(node_1 = sort(unique(test_all$node_1)),
                       node_1_fct = 1:472)
data_analyse <- left_join(x = test_all, y = eles_all, by = 'node_1')
eles_all <- data.frame(node_2 = sort(unique(test_all$node_1)),
                       node_2_fct = 1:472)
data_analyse <- left_join(x = data_analyse, y = eles_all, by = 'node_2')
head(data_analyse, 50)
str(data_analyse)
data_analyse$dyad_fct <- as.integer(data_analyse$dyad_fct)

#data_agg <- data_analyse %>%
#  group_by(node_1, node_2) %>%
#  summarise(event_count=sum(event), dyad_fct=cur_group_id())
#head(data_agg, 20) -- produces crazy high event counts -- summing event number rather than counting them

## TARGET:
##   node_1 node_2  event_count dyad_id node_1_id node_2_id
## 1 Rey    Leia             44       1         1         2
## 2 Rey    Obi-Wan          20       2         1         3
## 3 Rey    Luke             28       3         1         4
## 4 Rey    C-3PO             4       4         1         5
## 5 Rey    BB-8              0       5         1         6
## 6 Rey    R2-D2             0       6         1         7

#data_agg <- data_analyse %>%
#  group_by(dyad_fct) %>%
#  summarise(event_count=sum(event), dyad_fct=cur_group_id())
#head(data_agg, 20)  # same issue as above

event_count <- as.data.frame(table(data_analyse$dyad))
head(event_count, 20) ; tail(event_count,20)
colnames(event_count) <- c('dyad','event_count')
data_agg <- left_join(x = event_count, y = data_analyse, by = 'dyad')
data_agg_short <- data_agg[,c(1:6,9:11)]
data_agg_short <- distinct(.data = data_agg_short, .keep_all = T)
# well yes this worked, all apart from the very first step which I didn't notice before now gave every dyad the same score when it shared a node_1....

f1 <- data_analyse[data_analyse$node_1 == 'F1',]
table(f1$node_2)


##### STARTING OVER #####
sightings <- read_delim('data_processed/motnp_eles_long_correctid.csv', delim = ',') ; sightings <- sightings[,c(1:12,14)]
dyads <- read_delim('data_processed/motnp_elelinks.csv', delim = ',') ; dyads <- dyads[,c(1:3,5:13,18:20)]
target1 <- data.frame(node_1 = rep('M1',6),
                      node_2 = c(rep('M2',4), rep('M3',2)),
                      type_1 = rep('AM',6),
                      type_2 = c(rep('PM',4),rep('AM',2)),
                      event = c(1:6),
                      location = c('A','E','B','C','B','D'))
head(sightings) ; head(dyads) ; head(target1)

sightings$event_id <- as.factor(sightings$encounter)
head(sightings, 50)

df <- data.frame(node_1 = dyads$id1,
                 node_2 = dyads$id2,
                 type_1 = dyads$dem_class1,
                 type_2 = dyads$dem_class2)
s <- data.frame(event = sightings$event_id, 
                location = paste(sightings$gps_s, sightings$gps_e, sep='_'),
                individual = sightings$elephant)   # this is not node_1, this is both nodes
head(df) ; head(s) ; target1
# subset s by event --> if both individuals are present then include in event column
# for each event, find all individuals --> form dyads of those individuals --> add a new row per dyad --> fill in event details

encounters <- read_delim('data_processed/motnp_encounters_correctid.csv', delim = ',')
typ1 <- encounters[encounters$typ == 1,] ; typNA <- encounters[is.na(encounters$typ),]
first <- rbind(typNA, typ1) ; first <- first[!is.na(first$date),]
first$encounter <- 1:nrow(first)
first <- first[,c(102,8,26:101)]

myTriangular <- function(n) {  # number of dyads per group = triangular of group size
  k <- cumsum(seq(n)-1)
  return(k[n])
}
myTriangular(c(1:10))
myTriangular(first10$total_id_hkm[1])
myTriangular(first10$total_id_hkm[2])

first10 <- first[first$encounter < 11,]
first10$no_dyads <- NA
first10 <- first10[c(1,2,79,3:78)]
for(i in 1:nrow(first10)){
  first10$no_dyads[i] <- ifelse(first10$total_id_hkm[i] == 0, 0, myTriangular(first10$total_id_hkm[i]) )
}

encounter_number <- 1
for(i in 2:10){
  encounter_number <- paste(encounter_number, rep(first10$encounter[i], each = first10$no_dyads[i]))
}
encounter_number

s10 <- s[as.numeric(s$event) < 11,]
#####
m10 <- matrix(data = NA, nrow = 11, ncol = length(unique(s10$individual))+1)
colnames(m10) <- c('event', sort(unique(s10$individual)))
m10[1,] <- colnames(m10)
m10[c(2:11),1] <- 1:10

for(i in 2:nrow(m10)){
  for(j in 2:ncol(m10)){
    for(k in 1:nrow(s10)){
      m10[i,j] <- ifelse(m10[1,j] == s10$individual[k] & m10[i,1] == s10$event[k], 1, 0)
    }
  }
}
m10[1,15] # M1
m10[2,1]  # event 1
m10[1,15] == s10$individual[1] # TRUE
m10[2,1] == s10$event[1]       # TRUE
m10[2,15] <- ifelse(m10[1,15] == s10$individual[1] & m10[2,1] == s10$event[1], 1, 0) # 1
m10[2,15] <- NA

for(k in 1:nrow(s10)){
  m10[2,15] <- ifelse(m10[1,15] == s10$individual[k] & m10[2,1] == s10$event[k], 1, 0)
}

for(k in 1:nrow(s10)){
  for(i in 2:nrow(m10)){
    for(j in 2:ncol(m10)){
      m10[i,j] <- ifelse(m10[1,j] == s10$individual[k] & m10[i,1] == s10$event[k], 1, 0) 
    }
  }
}
#####
sightings <- read_delim('data_processed/motnp_eles_long_correctid.csv', delim = ',') ; sightings <- sightings[,c(1:12,14)]
target1 <- data.frame(node_1 = rep('M1',6),
                      node_2 = c(rep('M2',4), rep('M3',2)),
                      type_1 = rep('AM',6),
                      type_2 = c(rep('PM',4),rep('AM',2)),
                      event = c(1:6),
                      location = c('A','E','B','C','B','D'))
sightings$location <- paste(sightings$gps_s,sightings$gps_e,sep='_')
s <- sightings[,c(1,13,14)]
head(s) ; head(target1)
s10 <- s[s$encounter < 11,]
s10 <- s10[,c(1,2)]

#test <- as.data.frame(apply(s10, 2, FUN = combn, m = 2))
#length(test[test$encounter == ' 1',]$encounter)
#test <- lapply(seq_along(s10$elephant), combn, m = 2, simplify = FALSE)

#apply(x, MARGIN, FUN)	Apply a function to the rows or columns or both	                            Data frame or matrix	vector, list, array
#lapply(X, FUN)      	Apply a function to all the elements of the input	List, vector or data frame	list
#sapply(X, FUN)      	Apply a function to all the elements of the input	List, vector or data frame	vector or matrix

e2 <- s[s$encounter == 2,]
e2.eles <- e2$elephant
#test <- sapply(e2.eles, FUN = combn, m = 2, simplify = FALSE)
test <- do.call(c, lapply(seq_along(e2.eles), combn, x = e2.eles, simplify = FALSE))
View(test)

multi_combn_l <- function(dat, n) {
  unlist(lapply(1:n, function(x) combn(dat, x, simplify=F)), recursive=F)
}
test <- multi_combn_l(dat = e2.eles, n = 2) # produces all singles and doubles
df <- as.data.frame(test)
df2 <- data.frame(encounter = 2,
                  dyad = rep(NA, length(test)),
                  node_1 = rep(NA,length(test)),
                  node_2 = rep(NA,length(test)))
for(i in 1:nrow(df2)){
  df2$dyad[i]   <- as.character(test[[i]])
  df2$node_1[i] <- test[[i]]
  df2$node_2[i] <- test[[i]]
}

#####

library(tidyverse)
df1 <- read_delim('data_processed/mpnp_gbi_to_aggregated_dyad_test_firsthalf.csv', delim = ';')
df2 <- read_delim('data_processed/mpnp_gbi_to_aggregated_dyad_test_secondhalf.csv', delim = ';')

df <- rbind(df1,df2)

write_csv(df, 'data_processed/test.csv')
readr::write_csv(df, 'data_processed/test.csv') # may have tried to use the wrong package before when it wouldn't go in














### make data frame that matches required structure for asnipe::get_associations_points_tw(): "Input data must be of the following form: Date is an integer for day (usually starting at 1 on the first day). Time are the number of seconds elapsed from the start (continuous across all dates). ID is a unique character string for each individual. Location is a unique character string for each location."
eles_asnipe <- read_delim('data_processed/motnp_eles_long_correctid.csv', delim = ',')
eles_asnipe$location <- as.character(paste(eles_asnipe$gps_s, eles_asnipe$gps_e, sep = '_'))
eles_asnipe <- eles_asnipe[,c(2:3,14:15)]

eles_asnipe$time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong. 1 visit is recorded at time 1599 -- cannot be corrected, so convert to missing value.
eles_asnipe$time <- eles_asnipe$time*86400   # convert time values to seconds through day
which(is.na(eles_asnipe$time))      # 161  698 1122 1469 1770
eles_asnipe$time[c(161,698,1122,1469,1770)] <- 0

eles_asnipe$Date <- as.integer( eles_asnipe$date-min(eles_asnipe$date-1) )
eles_asnipe$Time <- eles_asnipe$time + (eles_asnipe$Date-1)*86400
eles_asnipe$Date <- as.integer(eles_asnipe$Date)
hist(eles_asnipe$Date, breaks = 50)
hist(eles_asnipe$Time, breaks = 50)

eles_asnipe <- eles_asnipe[,c(5,6,3,4)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
str(eles_asnipe)

which(is.na(eles_asnipe$Date))      # none
which(is.na(eles_asnipe$Time))      # none
which(is.na(eles_asnipe$ID))        # none
which(is.na(eles_asnipe$Location))  # none

a <- get_associations_points_tw(eles_asnipe,            # data to input
                                time_window = 60,       # >1 minute apart = different groups
                                which_days = NULL,      # all days
                                which_locations = NULL) # all locations
gbi <- a[[1]] # group-by-individual = first part of a

## THIS ISN'T WORKING -- TRY WITH TIME PER DAY RATHER THAN CONTINUOUS:
eles_asnipe <- read_delim('data_processed/motnp_eles_long_correctid.csv', delim = ',')
eles_asnipe$location <- as.character(paste(eles_asnipe$gps_s, eles_asnipe$gps_e, sep = '_'))
eles_asnipe <- eles_asnipe[,c(2:3,14:15)]

eles_asnipe$time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong. 1 visit is recorded at time 1599 -- cannot be corrected, so convert to missing value.
eles_asnipe$time <- eles_asnipe$time*86400   # convert time values to seconds through day
which(is.na(eles_asnipe$time))      # 161  698 1122 1469 1770
eles_asnipe$time[c(161,698,1122,1469,1770)] <- 0

eles_asnipe$Date <- as.integer( eles_asnipe$date-min(eles_asnipe$date-1) )

eles_asnipe <- eles_asnipe[,c(5,2:4)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)

a <- get_associations_points_tw(eles_asnipe,            # data to input
                                time_window = 60,       # >1 minute apart = different groups
                                which_days = NULL,      # all days
                                which_locations = NULL) # all locations
gbi <- a[[1]]                   # group-by-individual matrix = first part of a
times_per_group <- a[[2]]       # vector of times for each group = second part of a
dates_per_group <- a[[3]]       # vector of dates for each group
locat_per_group <- a[[4]]       # vector of locations for each group

## code from Dan to convert to usable format
obs <- t(sapply(1:20, function(x) rbinom(8, 1, 0.5)))
obs              # 8*20 matrix = 160 events
myTriangular(8)  # 28 possible dyads amongst 8 pairs
sum(obs)         # 84 observations
obs2 <- as.data.frame(obs)
obs2$count <- rowSums(obs)
obs2$tri.num <- myTriangular(obs2$count)
sum(obs2$tri.num) # 152
588/84            # 7 elephants?
588/28            # 21 interactions between 7 elephants?

df <- data.frame(node_1=numeric(), node_2=numeric(), social_event=numeric(), obs_id=numeric())
for (obs_id in 1:nrow(obs)) {
  for (i in which(obs[obs_id, ] == 1)) {
    for (j in 1:ncol(obs)) {
      if (i != j) {
        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        df[nrow(df) + 1, ] <- list(node_1=node_1, node_2=node_2, social_event=(obs[obs_id, i] == obs[obs_id, j]), obs_id=obs_id)
}}}}
str(df)
# node_1 = smaller value of the two individuals
# node_2 =  larger value of the two individuals
# social_event = 1 or 0, did it happen
# obs_id = event number (row number)

# DAN SAID THIS WOULD BE FOR CONVERTING GBI FORMAT TO TARGET FORMAT, BUT THE GBI OUTPUT I'VE GOT IS BY INDIVIDUALS IN BOTH DIRECTIONS, NOT BY EVENT NUMBER. HAVE I A) MISUNDERSTOOD "GBI", B) CREATED MY GBI WRONG, OR C) GOT THE WRONG CODE FROM DAN?

length(gbi[,1])           # 4262
length(eles_asnipe$Date)  # 4263
# gbi rows are events, but the first one is a pair of elephants so is removed -- rownames are based on ID column.
gbi2 <- gbi
eles_asnipe$event <- as.factor(paste(eles_asnipe$Date,eles_asnipe$Time,eles_asnipe$Location,sep = '_'))

#event <- data.frame(event = 1:4263,
#                    code = paste(eles_asnipe$Date,eles_asnipe$Time,eles_asnipe$Location,sep = '_'))
#which(event$code != eles_asnipe$event)

encounters <- data.frame(event = unique(eles_asnipe$event), event_id = 1:574)
asnipe <- left_join(x = eles_asnipe, y = encounters, by = 'event')
asnipe$event_id <- as.character(asnipe$event_id)
rownames(gbi2) <- asnipe$event_id[2:4263]

df <- data.frame(node_1=numeric(), node_2=numeric(), social_event=numeric(), obs_id=numeric())
for (obs_id in 1:nrow(gbi2)) {
  for (i in which(gbi2[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi2)) {
      if (i != j) {
        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        df[nrow(df) + 1, ] <- list(node_1=node_1, node_2=node_2, social_event=(gbi2[obs_id, i] == gbi[obs_id, j]), obs_id=obs_id)
      }
    }
  }
}
df

write_csv2(df, 'data_processed/mpnp_gbi_to_aggregated_dyad_test.csv') # only able to save the first 1048576

df1 <- df[c(1:1048576),]
2096563-1048576 # 1047987 -- will JUST fit into 2 Excel sheets
df2 <- df[c(1048577:2096563),]

write_csv2(df1, 'data_processed/mpnp_gbi_to_aggregated_dyad_test_firsthalf.csv')
write_csv2(df2, 'data_processed/mpnp_gbi_to_aggregated_dyad_test_secondhalf.csv')









target1 <- data.frame(node_1 = c('M1', rep('F9', 8), rep('U8',2)),
                      node_2 = c('M11','U8','F8','F13','M26','M87','F19','M19','M24','F8','F13'),
                      type_1 = c('AM', rep('AF',8), rep('CU',2)),
                      type_2 = c('AM','CU','AF','PF','JM','PM','AF','PM','CM','AF','PF'),
                      event = c(1,rep(2,10)),
                      location = c('A',rep('B',10)))

##### starting over again -- 4th January 2022 #####
## sightings
sightings <- read_delim(file = 'data_processed/motnp_encounters_correctid.csv', delim = ',')
sightings <- s[1:22]
str(sightings)

## elephant encounters
eles <- read_delim(file = 'data_processed/motnp_eles_long_correctid.csv', delim = ',')
eles$location <- paste(eles$gps_s, eles$gps_e, sep = '_')
colnames(eles)
eles <- eles[c(1,14,2,3,15,12,6,8,11)]
str(eles)

## nodes
nodes <- read_delim(file = 'data_processed/motnp_elenodes.csv', delim = ',')
colnames(nodes)
str(nodes)

## links
links <- read_delim('data_processed/motnp_elelinks.csv', delim = ',')

##### Dan code to convert gbi to target format #####
---
  title: "Convert group-by-individual data to long format"
output:
  rmarkdown::github_document
---
  
  ## Simulate some data
  
  ```{r}
obs <- t(sapply(1:20, function(x) rbinom(8, 1, 0.5)))
obs
```

## Convert to long format

```{r}
df <- data.frame(node_1=numeric(), node_2=numeric(), social_event=numeric(), obs_id=numeric())
for (obs_id in 1:nrow(obs)) {
  for (i in which(obs[obs_id, ] == 1)) {
    for (j in 1:ncol(obs)) {
      if (i != j) {
        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        df[nrow(df) + 1, ] <- list(node_1=node_1, node_2=node_2, social_event=(obs[obs_id, i] == obs[obs_id, j]), obs_id=obs_id)
      }
    }
  }
}
df
```










