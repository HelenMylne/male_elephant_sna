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


library(tidyverse)
df1 <- read_delim('data_processed/mpnp_gbi_to_aggregated_dyad_test_firsthalf.csv', delim = ';')
df2 <- read_delim('data_processed/mpnp_gbi_to_aggregated_dyad_test_secondhalf.csv', delim = ';')

df <- rbind(df1,df2)

write_csv(df, 'data_processed/test.csv')
readr::write_csv(df, 'data_processed/test.csv') # may have tried to use the wrong package before when it wouldn't go in


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

##### Dan Rmd code to convert gbi to target format #####
#---
#title: "Convert group-by-individual data to long format"
#output:
#rmarkdown::github_document
#---

## Simulate some data
#```{r}
#obs <- t(sapply(1:20, function(x) rbinom(8, 1, 0.5)))
#obs
#```

## Convert to long format
#```{r}
#df <- data.frame(node_1=numeric(), node_2=numeric(), social_event=numeric(), obs_id=numeric())
#for (obs_id in 1:nrow(obs)) {
#  for (i in which(obs[obs_id, ] == 1)) {
#    for (j in 1:ncol(obs)) {
#      if (i != j) {
#        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
#        if (i < j) {
#          node_1 <- i
#          node_2 <- j
#        } else {
#          node_1 <- j
#          node_2 <- i
#        }
#        df[nrow(df) + 1, ] <- list(node_1=node_1, node_2=node_2, social_event=(obs[obs_id, i] == obs[obs_id, j]), obs_id=obs_id)
#      }
#    }
#  }
#}
#df
#```

##### make short version of data to test Dan's code #####
eles_short <- eles[eles$encounter == 2 | eles$encounter == 6 | eles$encounter == 1,]
length(unique(eles_short$elephant)) # 22 rows, 18 elephants -- F13, F8, M87, M26 seen twice

# create gbi
eles_asnipe <- eles_short[,c(3,4,2,5)]
eles_asnipe
eles_asnipe$Date <- as.integer(eles_asnipe$date)
eles_asnipe$Date <- 1+eles_asnipe$Date - min(eles_asnipe$Date)
eles_asnipe$Time <- ifelse(eles_asnipe$time > 1, NA, eles_asnipe$time) # time = proportion of day so anything >1 has to be wrong. 1 visit is recorded at time 1599 -- cannot be corrected, so convert to missing value.
eles_asnipe$Time <- eles_asnipe$time*(24*60*60)   # convert time values to seconds through day
eles_asnipe$Time <- eles_asnipe$Time + (eles_asnipe$Date-1)*(24*60*60)
which(is.na(eles_asnipe$Time))

eles_asnipe <- eles_asnipe[,c(5,6,3,4)]
colnames(eles_asnipe) <- c('Date','Time','ID','Location')
eles_asnipe$ID <- as.character(eles_asnipe$ID)
str(eles_asnipe)

### make data frame that matches required structure for asnipe::get_associations_points_tw(): "Input data must be of the following form: Date is an integer for day (usually starting at 1 on the first day). Time are the number of seconds elapsed from the start (continuous across all dates). ID is a unique character string for each individual. Location is a unique character string for each location."
a <- get_associations_points_tw(eles_asnipe,            # data to input
                                time_window = 60,       # >1 minute apart = different groups
                                which_days = NULL,      # all days
                                which_locations = NULL) # all locations
gbi <- a[[1]]  # group-by-individual = first part of a -- gbi is empty
test <- a[[2]] # empty
test <- a[[3]] # empty
test <- a[[4]] # empty

## THIS ISN'T WORKING -- TRY WITH TIME PER DAY RATHER THAN CONTINUOUS:
eles_asnipe$Time <- eles_short$time
eles_asnipe$Time <- eles_asnipe$Time*(24*60*60)   # convert time values to seconds through day

a <- get_associations_points_tw(eles_asnipe,            # data to input
                                time_window = 60,       # >1 minute apart = different groups
                                which_days = NULL,      # all days
                                which_locations = NULL) # all locations
a # still empty

##### alternative way to generate gbi #####
#get_gbi generates a group by individual matrix. The function accepts a data.table with individual identifiers and a group column. The group by individual matrix can then be used to build a network using asnipe::get_network.
DT <- as.data.table(eles_asnipe)
DT$group <- as.integer(as.factor(DT$Time))
DT <- DT[,c(3,5)]
DT

gbi <- spatsoc::get_gbi(DT = DT, group = 'group', id = 'ID')

## 
obs <- t(sapply(1:20, function(x) rbinom(8, 1, 0.5)))
obs

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
# node_1 = smaller value of the two individuals
# node_2 =  larger value of the two individuals
# social_event = 1 or 0, did it happen
# obs_id = event number (row number)

gbi.df <- data.frame(node_1=numeric(), node_2=numeric(), social_event=numeric(), obs_id=numeric())
for (obs_id in 1:nrow(gbi)) {
  for (i in which(gbi[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi)) {
      if (i != j) {
        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi.df[nrow(gbi.df) + 1, ] <- list(node_1=node_1, node_2=node_2, social_event=(gbi[obs_id, i] == gbi[obs_id, j]), obs_id=obs_id)
      }
    }
  }
}
gbi.df
## seems to have worked, but need to check if numbers in gbi.df match to elephant IDs
str(gbi.df) ; str(gbi.test)
gbi.test <- data.frame(id_1 = colnames(gbi), node_1 = as.numeric(1:18))
gbi.check <- left_join(x = gbi.df, y = gbi.test, by = 'node_1')
gbi.test <- data.frame(id_2 = colnames(gbi), node_2 = as.numeric(1:18))
gbi.check <- left_join(x = gbi.check, y = gbi.test, by = 'node_2')
## check by eye -- should be:
# 1 for M1 and M11 and no others for that pairing -- yes
# 1 for all pairs in obs_id 2 and no others -- yes
# 1 for all pairs in obs_id 3 and no others -- yes

## all duplicated, and obs_id ≠ encounter number (encounters in test = 1, 2 and 6, some numbers missing in raw data, to avoid confusion will need to renumber one or the other)
gbi.check$unique <- paste(gbi.check$node_1, gbi.check$node_2, gbi.check$social_event, gbi.check$obs_id, sep = '_')
length(unique(gbi.check$unique)) # 282 (total length = 374)

gbi.distinct <- dplyr::distinct(gbi.check)

encounter.numbers <- data.frame(encounter = unique(eles_short$encounter), obs_id = 1:3)
gbi.check <- left_join(x = gbi.distinct, y = encounter.numbers, by = 'obs_id')  # perfect!

### TARGET FORM:
##   node_1 node_2 type_1   type_2   event location dyad_id location_id
## 1 Rey    Leia   Lifeform Lifeform     1 E              1           5
## 2 Rey    Leia   Lifeform Lifeform     1 C              1           3
## 3 Rey    Leia   Lifeform Lifeform     1 D              1           4
## 4 Rey    Leia   Lifeform Lifeform     1 D              1           4
## 5 Rey    Leia   Lifeform Lifeform     1 C              1           3
## 6 Rey    Leia   Lifeform Lifeform     1 F              1           6

### CURRENT FORM:
head(gbi.check)
##   node_1 node_2 social_event obs_id id_1 id_2         unique encounter
## 1      1      6            0      1  F13   M1 1_6_0_1_F13_M1         1
## 2      2      6            0      1  F19   M1 2_6_0_1_F19_M1         1
## 3      3      6            0      1  F21   M1 3_6_0_1_F21_M1         1
## 4      4      6            0      1   F8   M1  4_6_0_1_F8_M1         1
## 5      5      6            0      1   F9   M1  5_6_0_1_F9_M1         1
## 6      6      7            1      1   M1  M11 6_7_1_1_M1_M11         1

gbi <- data.frame(node_1 = gbi.check$id_1,
                  node_2 = gbi.check$id_2,
                  node_1_id = gbi.check$node_1,
                  node_2_id = gbi.check$node_2,
                  event = gbi.check$social_event,
                  encounter = gbi.check$encounter,
                  obs_id = gbi.check$obs_id)
                  #type_1 = NA,
                  #type_2 = NA,
                  #dyad_id = NA,
                  #location_id = NA)

head(eles, 10)  # individual eles, info on encounter
head(nodes, 10) # age_class = A/P/J/C, age_category = 50+/35-40/25-35/..., sex = M/F/U, dem_class = AF/AM/PF/PM/...

colnames(eles_short)  # "encounter","elephant","date","time","location","type","total_elephants","total_id_hkm","perc_id_hkm"
gbi.all <- left_join(x = gbi, y = eles_short, by = 'encounter')
gbi.all <- gbi.all[,c(1:7,9:15)]
names <- colnames(gbi.all) ; names[11] <- 'herd_type' ; colnames(gbi.all) <- names

nodes$node_1 <- nodes$id
unique(gbi$node_1)
nodes_short <- nodes[nodes$id == 'M1' | nodes$id == 'M3' | nodes$id == 'M5' | nodes$id == 'M8' | nodes$id == 'M9' |
                       nodes$id == 'M11' | nodes$id == 'M15' | nodes$id == 'M19' | nodes$id == 'M24' | nodes$id == 'M26' |
                       nodes$id == 'M52' | nodes$id == 'M87' | nodes$id ==  'F8' | nodes$id == 'F9' | nodes$id == 'F13' |
                       nodes$id == 'F19' | nodes$id == 'F21' | nodes$id == 'U8', ]
nodes_short$node_1 <- nodes_short$id
colnames(nodes_short) # "id_no","herd","name","age_class","age_category","family","days","comments","sex","num","id","count","dem_class","node_1"
gbi.all <- left_join(x = gbi.all, y = nodes_short, by = 'node_1')
gbi.all <- gbi.all[,c(1:14,17:19,23,26:27)]
names <- colnames(gbi.all)
names[c(15:20)] <- c('name_1','age_class_1','age_category_1','sex_1','count_1','dem_class_1')
colnames(gbi.all) <- names

nodes_short$node_2 <- nodes_short$id
nodes_short <- nodes_short[,c(1:13,15)]
gbi.all <- left_join(x = gbi.all, y = nodes_short, by = 'node_2')
colnames(gbi.all)
gbi.all <- gbi.all[,c(1:20,23:25,29,32:33)]
names <- colnames(gbi.all)
names[c(21:26)] <- c('name_2','age_class_2','age_category_2','sex_2','count_2','dem_class_2')
colnames(gbi.all) <- names

colnames(gbi.all)
gbi.all <- gbi.all[,c(1:15,21,16,22,17,23,18,24,19,25,20,26)]

##### remove duplicates, convert to final target format #####
### TARGET FORM:
##   node_1 node_2 type_1   type_2   event location dyad_id location_id
## 1 Rey    Leia   Lifeform Lifeform     1 E              1           5
## 2 Rey    Leia   Lifeform Lifeform     1 C              1           3
## 3 Rey    Leia   Lifeform Lifeform     1 D              1           4
## 4 Rey    Leia   Lifeform Lifeform     1 D              1           4
## 5 Rey    Leia   Lifeform Lifeform     1 C              1           3
## 6 Rey    Leia   Lifeform Lifeform     1 F              1           6

gbi.all$dyad <- paste(gbi.all$node_1, gbi.all$node_2, sep = '_')
gbi.all$dyad_id <- as.integer(as.factor(gbi.all$dyad))
gbi.all$location_id <- as.integer(as.factor(gbi.all$location))

gbi.short <- distinct(gbi.all)
gbi.dyads <- gbi.short[gbi.short$event == 1,]

df_agg <- gbi.short %>%
  group_by(node_1, node_2) %>%
  summarise(event_count=sum(event), dyad_id=cur_group_id()) %>%
  mutate(node_1_id_int=as.integer(as.factor(node_1)), node_2_id_int=as.integer(as.factor(node_2)))
head(df_agg)
##   node_1 node_2 event_count dyad_id node_1_id_int node_2_id_int
##   <chr>  <chr>        <dbl>   <int>         <int>         <int>
## 1 F13    F19              1       1             1             1
## 2 F13    F21              1       2             1             2
## 3 F13    F8               2       3             1             3
## 4 F13    F9               1       4             1             4
## 5 F13    M1               0       5             1             5
## 6 F13    M11              0       6             1             6
tail(df_agg)
##   node_1 node_2 event_count dyad_id node_1_id_int node_2_id_int
##   <chr>  <chr>        <dbl>   <int>         <int>         <int>
## 1 M8     M87              1     148             1             1
## 2 M8     M9               1     149             1             2
## 3 M8     U8               0     150             1             3
## 4 M87    M9               1     151             1             1
## 5 M87    U8               1     152             1             2
## 6 M9     U8               0     153             1             1
unique(df_agg$node_1_id_int)
unique(df_agg$node_2_id_int)
# yes the counts are correct, have to check dyad_id is correct, node_1_id_int all = 1, node_2_id_int 1-17 but missing 18. Too many are 1:1.
unique(df_agg$dyad_id)
df_agg$dyad <- paste(df_agg$node_1, df_agg$node_2, sep = '_')
which(df_agg$dyad == 'F13_F19')     # 1
which(gbi.short$dyad == 'F13_F19')  # 34 151
df_agg$dyad_id[1] == gbi.short$dyad_id[34]
df_agg$dyad_id[1] == gbi.short$dyad_id[151]
which(df_agg$dyad == 'M1_M11')      # 76
which(gbi.short$dyad == 'M1_M11')   # 6
df_agg$dyad_id[76] == gbi.short$dyad_id[6]
which(df_agg$dyad == 'M15_M24')     # 100
which(gbi.short$dyad == 'M15_M24')  # 112 204
df_agg$dyad_id[100] == gbi.short$dyad_id[112]
df_agg$dyad_id[100] == gbi.short$dyad_id[204]
# yes, dyad_id is correct

df_agg <- df_agg[,c(1:4)]
df_agg_corrected <- left_join(x = df_agg, y = gbi.short, by = 'dyad_id')
df_agg_corrected <- df_agg_corrected[,c(1:4,7,8)]
colnames(df_agg_corrected) <- c('node_1','node_2','event_count','dyad_id','node_1_id','node_2_id')

## I THINK this is now correct!

model_data <- list(
  N = nrow(gbi.short),                    # Number of observations
  M = nrow(df_agg),                       # Number of dyads
  L = length(unique(gbi.short$location)), # Number of locations
  dyad_ids = gbi.short$dyad_id,           # Vector of dyad IDs corresponding to each observation
  location_ids = gbi.short$location_id,   # Vector of location IDs corresponding to each observation
  event = gbi.short$event                 # Vector of binary values (0/1, presence/absence) corresponding to each observation
)

rm(a, df, df_agg, DT, eles_asnipe, encounter.numbers, gbi.check, gbi.df, gbi.distinct, gbi.dyads, gbi.test, obs, i,j,names,node_1,node_2,obs_id)
##### put into model to check output for 3 encounters only -- NOTE: I HAVEN'T WORKED THORUGH ANY DAGS OR ANYTHING HERE, JUST CAUSAL-SALAD-ED TO SEE IF I CAN REMEMBER HOW TO FIT A MODEL! #####
#model <- stan_model("../models/binary_model.stan")
#fit <- sampling(model, model_data, cores=4)
library(rstan)
library(rethinking)
model <- ulam(alist(
  event ~ dbinom(N, p),
  logit(p) <- beta_dyad[dyad_ids] + beta_loc[location_ids],
  beta_dyad[dyad_ids] ~ dnorm(0,2),
  beta_loc[location_ids] ~ dnorm(0,2)),
data = model_data, cores = 4, chains = 1, iter = 4000)
precis(model, 2)
traceplot(model) # that looks like a whole load of hairy caterpillars to me!
stancode(model)
# maybe I'm being overoptimistic, but this does actually look how I would expect it to! Not especially useful I don't think, but a start!

## try again, this time incorporating some information about the dyads other than just ID
head(df_agg_corrected)
df <- left_join(x = df_agg_corrected, y = gbi.short, by = 'dyad_id')
df <- distinct(df)
df$node_1.x == df$node_1.y
df$node_2.x == df$node_2.y
df <- df[,c(1:6,11:34)]
names <- colnames(df)
names[c(1,2,5,6)] <- c('node_1','node_2','node_1_id','node_2_id')
colnames(df) <- names
df$dem_type <- paste(df$dem_class_1, df$dem_class_2, sep = '_')
df$dem_type_id <- as.integer(as.factor(df$dem_type))
df$sex_type <- paste(df$sex_1, df$sex_2, sep = '_')
df$sex_type_id <- as.integer(as.factor(df$sex_type))
df$herd_type_id <- as.integer(as.factor(df$herd_type))

model_data_2 <- list(
  N = nrow(df),                    # Number of observations
  M = nrow(df),                    # Number of dyads
  L = length(unique(df$location)), # Number of locations
  dyad_ids = df$dyad_id,           # Vector of dyad IDs corresponding to each observation
  location_ids = df$location_id,   # Vector of location IDs corresponding to each observation
  event = df$event_count,          # Vector of binary values (0/1, presence/absence) corresponding to each observation
  dem_type = df$dem_type_id,       # add further information about sightings to model data
  sex_type = df$sex_type_id,
  herd_type = df$herd_type_id,
  total_count = df$total_elephants
)

model2 <- ulam(alist(
  event ~ dbinom(N, p),
  logit(p) <- b_dem[dem_type] + b_sex[sex_type] + b_herd[herd_type] + b_size*total_count,
  b_dem[dem_type] ~ dnorm(0,2),
  b_sex[sex_type] ~ dnorm(0,2),
  b_herd[herd_type] ~ dnorm(0,2),
  b_size ~ dnorm(0,2)),
  data = model_data_2, cores = 4, chains = 1, iter = 4000)
precis(model2,2)
traceplot(model2) # looking good! :D -- now do it again, but properly and actually follow all of the steps including producing DAGs, standardizing the variables, simulating data, performing prior predictive checks, checking the model on simulated data, THEN running it like this!
stancode(model2)

dev.off()
##### repeat with complete data, not just 3 encounters #####
## sightings
sightings <- read_delim(file = 'data_processed/motnp_encounters_correctid.csv', delim = ',')
sightings <- sightings[1:22]
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

## create gbi
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
DT <- as.data.table(eles_asnipe)
DT$d_pad <- str_pad(DT$Date, 3, pad = '0')
DT$encounter <- paste(DT$d_pad, DT$Time, DT$Location, sep = '_')
DT$group <- as.integer(as.factor(DT$encounter))
max(DT$group) # 574
DT <- DT[,c(3,7)]
DT

################ NOTE HERE -- NEED TO CHECK ELES_LONG DATAFRAME -- TOO FEW ENCOUNTERS OVERALL ################
gbi <- spatsoc::get_gbi(DT = DT, group = 'group', id = 'ID') ## THIS STEP WORKS ASSUMING THAT DT IS CORRECT. POTENTIAL ISSUE IS THAT THIS APPEARS TO INDICATE 574 SIGHTINGS, BUT SIGHITNGS SHEET CONTAINS 1313 ROWS OF DATA. IS THIS AN ERROR IN ELES_LONG??
length(unique(eles$encounter)) # 574 unique values
################ NOTE HERE -- NEED TO CHECK ELES_LONG DATAFRAME -- TOO FEW ENCOUNTERS OVERALL ################

myTriangular(472) * 2 # 111156 * 2 rows of data will be produced... I think...
gbi.df <- data.frame(node_1=numeric(), node_2=numeric(), social_event=numeric(), obs_id=numeric())
for (obs_id in 1:nrow(gbi)) {
  for (i in which(gbi[obs_id, ] == 1)) {
    for (j in 1:ncol(gbi)) {
      if (i != j) {
        # Hacky bit to make sure node_1 < node_2, not necessary but makes things a bit easier.
        if (i < j) {
          node_1 <- i
          node_2 <- j
        } else {
          node_1 <- j
          node_2 <- i
        }
        gbi.df[nrow(gbi.df) + 1, ] <- list(node_1=node_1, node_2=node_2, social_event=(gbi[obs_id, i] == gbi[obs_id, j]), obs_id=obs_id)
      }
    }
  }
}
gbi.df
write_delim(gbi.df, 'data_processed/gbi.df_first_run_2005989obs_4263individualencounters_22.01.06.csv', delim = ',')
#test <- read_delim('data_processed/gbi.df_first_run_2005989obs_4263individualencounters_22.01.06.csv', delim = ',')

## seems to have worked, but need to check if numbers in gbi.df match to elephant IDs
str(gbi.df)
gbi.test <- data.frame(id_1 = colnames(gbi), node_1 = as.numeric(1:472))
gbi.check <- left_join(x = gbi.df, y = gbi.test, by = 'node_1')
gbi.test <- data.frame(id_2 = colnames(gbi), node_2 = as.numeric(1:472))
gbi.check <- left_join(x = gbi.check, y = gbi.test, by = 'node_2')

# check
unique(eles_short$elephant)  # "M1","F9","F8","M11","U8","F13","F21","M15","M26","M87","F19","M3","M19","M5","M24","M8","M9","M52"
gbi.check.e1 <- gbi.check[gbi.check$obs_id == 1,]
gbi.check.e1[which(gbi.check.e1$social_event == 1),]             # M1 and M11 -- this is correct
gbi.check.e2 <- gbi.check[gbi.check$obs_id == 2,]
unique(gbi.check.e2$id_1[which(gbi.check.e2$social_event == 1)]) # "F13" "F19" "F8"  "F9"  "M19" "M24" "M26" "M87" - correct
unique(gbi.check.e2$id_2[which(gbi.check.e2$social_event == 1)]) # "F19" "F8"  "F9"  "M19" "M24" "M26" "M87" "U8"  - correct
gbi.check.e6 <- gbi.check[gbi.check$obs_id == 6,]    # 15073-15543
gbi.check.e6[which(gbi.check.e6$social_event == 1),] # all are 0 -- seems to be showing everyone against F40 and then F40 against everyone else, but that's it? -- because I've had to cut out a load of sightings at some point. Encounter numbers go 1,2,3,5,6,8, so the 6th registered encounter is actually encounter 8 which is 10 elephants of which only F40 was identified.
gbi.check.e6 <- gbi.check[gbi.check$obs_id == 5,]
unique(gbi.check.e6$id_1[which(gbi.check.e6$social_event == 1)]) # "F13","F21","F8","M15","M26","M3","M5","M52","M8","M87" - correct
unique(gbi.check.e6$id_2[which(gbi.check.e6$social_event == 1)]) # "F21","F8","M15","M26","M3","M5","M52","M8","M87","M9"  - correct

eles$obs_id <- as.integer(as.factor(eles$encounter))
eles <- eles[,c(1,14,2:13)]
gbi.encounter <- left_join(x = gbi.check, y = eles, by = 'obs_id')
length(unique(gbi.encounter$encounter)) # 574 -- correct from input data, but input should have been 1313

gbi.encounter$unique <- paste(gbi.encounter$node_1, gbi.encounter$node_2, gbi.encounter$social_event, gbi.encounter$obs_id, sep = '_')
length(unique(gbi.encounter$unique))           # 1963129 (total length = 42401775)
gbi.distinct <- dplyr::distinct(gbi.encounter) # 40708781 obs -- this is way too long, this shouldn't be happening -- too many to load, can't inspect visually
colnames(gbi.distinct) # "node_1","node_2","social_event","obs_id","id_1","id_2","encounter","date","time","gps_s","gps_e","total_elephants","total_id_dy","total_id_hkm","total_id_diff","perc_id_dy","perc_id_hkm","type","elephant","unique"
head(gbi.distinct, 10)
gbi.distinct <- gbi.distinct[,c(1:12,14,17,18)]

names <- colnames(gbi.distinct) ; names[c(7,13:15)] <- c('encounter_id','total_id','perc_id','herd_type') ; colnames(gbi.distinct) <- names
gbi.distinct$location <- paste(gbi.distinct$gps_s, gbi.distinct$gps_e, sep = '_')

## convert to model data format
md1 <- gbi.distinct[,c(1:9,16,12:15)]
md1$dyad <- paste(md1$id_1, md1$id_2, sep = '_')
md1$dyad_id <- as.integer(as.factor(md1$dyad))
md1$location_id <- as.integer(as.factor(md1$location))

md <- distinct(md1)
# md_dyads <- md[md$social_event == 1,] # not sure whether this is necessary
head(md)

df_agg <- md %>%
  group_by(id_1, id_2) %>%
  summarise(event_count=sum(social_event), dyad_id=cur_group_id()) %>%
  mutate(node_1_id=as.integer(as.factor(id_1)), node_2_id=as.integer(as.factor(id_2)))
length(df_agg$id_1) == myTriangular(472) # correct number of dyads!
head(df_agg)
##   id_1  id_2  event_count dyad_id node_1_id node_2_id
##   <chr> <chr>       <dbl>   <int>     <int>     <int>
## 1 F1    F10             0       1         1         1
## 2 F1    F100            0       2         1         2
## 3 F1    F101            0       3         1         3
## 4 F1    F102            0       4         1         4
## 5 F1    F103            0       5         1         5
## 6 F1    F104            0       6         1         6
tail(df_agg)
##  id_1  id_2  event_count dyad_id node_1_id node_2_id
##   <chr> <chr>       <dbl>   <int>     <int>     <int>
## 1 U67   U7              0  111151         1         1
## 2 U67   U8              0  111152         1         2
## 3 U67   U9              0  111153         1         3
## 4 U7    U8              1  111154         1         1
## 5 U7    U9              1  111155         1         2
## 6 U8    U9              2  111156         1         1

unique(df_agg$node_1_id)     # STILL GOT ALL ARE 1
unique(df_agg$node_2_id)     # 1-471 BUT head() AND tail() SHOW THAT node_2_id STILL HAS VALUES 1-3 FOR U7-U9 AND F10,F100 AND F101

length(unique(df_agg$dyad_id)) # 111156 -- correct
df_agg$dyad <- paste(df_agg$id_1, df_agg$id_2, sep = '_')
which(df_agg$dyad == 'F1_F10')        # 1
which(md1$dyad == 'F1_F10')           # 546 ENTRIES
df_agg$dyad_id[1] == md1$dyad_id[40181]
df_agg$dyad_id[1] == md1$dyad_id[14930761]
df_agg$dyad_id[1] == md1$dyad_id[39958855]
which(df_agg$dyad == 'M1_M11')        # 60765
which(md1$dyad == 'M1_M11')           # 500 ENTRIES
df_agg$dyad_id[60765] == md1$dyad_id[331]
df_agg$dyad_id[60765] == md1$dyad_id[15966960]
df_agg$dyad_id[60765] == md1$dyad_id[40513969]
which(df_agg$dyad == 'M15_M24')       # 76800
which(md1$dyad == 'M15_M24')          # 481 ENTRIES
df_agg$dyad_id[76800] == md1$dyad_id[24833]
df_agg$dyad_id[76800] == md1$dyad_id[11469918]
df_agg$dyad_id[76800] == md1$dyad_id[40194294]
# yes, dyad_id is correct

df_agg <- df_agg[,c(1:4)]
## tested way to correct node_1 and node_2 not working due to exceeded vector memory
#df_agg_corrected <- left_join(x = df_agg, y = md1, by = 'dyad_id')
#df_agg_corrected <- df_agg_corrected[,c(1:4,7,8)]
#colnames(df_agg_corrected) <- c('node_1','node_2','event_count','dyad_id','node_1_id','node_2_id')
## alternative way
df_agg$node_1 <- as.integer(as.factor(df_agg$id_1))
df_agg$node_2 <- as.integer(as.factor(df_agg$id_2))+1 # add 1 so starts at 2 and F1 is "1"
head(df_agg,10) ; tail(df_agg,10)
## I THINK this is now correct!

## add data about nodes
colnames(nodes)
node_data <- nodes[is.na(nodes$comments),c(1,3:5,9,11:13)]
node_data$id_1 <- node_data$id
node_data$id_2 <- node_data$id

unique(gbi.distinct$id_1) # 155 females, 250 males, 65 unknowns
unique(node_data$id_1)    # 151 females, 245 males, 66 unknowns

colnames(node_data) ; colnames(df_agg)
dyads <- left_join(x = df_agg, y = node_data, by = 'id_1')
names <- colnames(dyads)
names[c(2,7:15)] <- c('id_2','id_pad_1','name_1','age_class_1','age_category_1','sex_1','id1_deletecolumn','count_1','dem_class_1','deletecolumn1')
colnames(dyads) <- names
head(dyads)

dyads <- left_join(x = dyads, y = node_data, by = 'id_2')
names <- colnames(dyads)
names[c(1,16:24)] <- c('id_1','id_pad_2','name_2','age_class_2','age_category_2','sex_2','id2_deletecolumn','count_2','dem_class_2','deletecolumn2')
colnames(dyads) <- names

str(dyads)
dyads <- dyads[,c(4,1,2,5,6,3,7,16,8,17,9,18,10,19,11,20,13,22,14,23)]
str(dyads)

## create model data
model_data <- list(
  N = nrow(gbi.short),                    # Number of observations
  M = nrow(df_agg),                       # Number of dyads
  L = length(unique(gbi.short$location)), # Number of locations
  dyad_ids = gbi.short$dyad_id,           # Vector of dyad IDs corresponding to each observation
  location_ids = gbi.short$location_id,   # Vector of location IDs corresponding to each observation
  event = gbi.short$event                 # Vector of binary values (0/1, presence/absence) corresponding to each observation
)


