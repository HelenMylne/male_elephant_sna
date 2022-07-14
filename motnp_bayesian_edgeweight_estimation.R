## First attempt to roughly build a Bayesian model of the ALERT data
##### set up and data import ####
library(tidyverse)
library(rstan)
library(dplyr)
library(igraph)
## sightings
s <- read_delim(file = 'data_processed/motnp_encounters_correctid.csv', delim = ',')
s <- s[1:22]
str(s)

## elephant encounters
e <- read_delim(file = 'data_processed/motnp_eles_long_correctid.csv', delim = ',')
colnames(e)
e <- e[c(1:12,14)]
str(e)

## nodes
n <- read_delim(file = 'data_processed/motnp_elenodes.csv', delim = ',')
colnames(n)
str(n)

## links
el <- read_delim('data_processed/motnp_elelinks.csv', delim = ',')

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

##### example from Dan #####
# First of all we’ll load in Rstan for model fitting in Stan, dplyr for handling the data, and igraph for network plotting and computing network centrality. We also load in two custom R files: “simulations.R” to generate synthetic data for this example; and “sampler.R” to allow fitting models with uncertainty over network features respectively.
source("../scripts/simulations.R")
source("../scripts/sampler.R")

# Now we will simulate data using the simulate_binary() function. The rows of the resulting dataframe describe observations at the dyadic level between nodes. In this dataframe, event denotes whether or not an undirected social event was observed in this observation period. The exact definition of observation period will depend on the study, but is commonly a sampling period where at least one of the members of the dyad was observed. This can also be a sampling period where both members of the dyad were observed, and the distinction will affect the interpretation of social preferences. See the paper for further discussion on this. location denotes the location at which the observation took place, which may be relevant if location is likely to impact the visibility of social events.
set.seed(123)
data <- simulate_binary()
df <- data$df
head(df)

# Computationally it’s easier to work with dyad IDs rather than pairs of nodes in the statistical model, so we’ll map the pairs of nodes to dyad IDs before we put the data into the model. The same is true for the location factor, so we will also map the locations to location IDs. We can add these columns to the dataframe using the following code:
df <- df %>%
  group_by(node_1, node_2) %>%
  mutate(dyad_id=cur_group_id()) %>%
  mutate(location_id=as.integer(location))
head(df)

df_agg <- df %>%
  group_by(node_1, node_2) %>%
  summarise(event_count=sum(event), dyad_id=cur_group_id()) %>%
  mutate(node_1_id=as.integer(node_1), node_2_id=as.integer(node_2))
head(df_agg)


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

#df <- gather(data = test2, key = node_1, value = node_2, encounter)
#head(df)

test2$location <- paste(test2$gps_s, test2$gps_e, sep = '_')
test2 <- test2[c(1:3,14,6:13)]
head(test2)

sri <- read_delim('data_processed/motnp_dyads_noduplicates.csv', delim = ',')
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

eles <- sort(unique(events$elephant))
eles

m <- matrix(data = NA, nrow = 2, ncol = 11)
colnames(m) <- eles
pivot_wider(data = events, names_from = elephant, values_from = event)
#   location           M1    F9   M11    U8    F8   F13   M26   M87   F19   M19   M24
# 1 1752912_2549597     1    NA     1    NA    NA    NA    NA    NA    NA    NA    NA
# 2 1752031_2547860    NA     2    NA     2     2     2     2     2     2     2     2

# find values node_1[i] and node_2[i] in eles (j and k) --> look through columns m[,j] and m[,k] --> 
## ifelse(!is.na(m[l,j]) & !is.na(m[l,l]), put value l in dyad$event1, NA),

#test_short <- pivot_wider(data = events, names_from = event, values_from = elephant)
test_short <- pivot_wider(data = events, names_from = elephant, values_from = event)
eles <- colnames(test_short)

# dyad$event_count -- I can get a column for total events that the dyad were observed in together
for(i in 1:length(dyads$node_1)) {
  j <- match(dyads$node_1[i], eles)
  k <- match(dyads$node_2[i], eles)
  dyads$event_count[i] <- length(which(test_short[,j]==test_short[,k]))
}

# dyad$event_id -- haven't yet worked out how to generate a column for specific events where each pair were observed.
for(i in 1:length(dyads$node_1)) {
  j <- match(dyads$node_1[i], eles)
  k <- match(dyads$node_2[i], eles)
  for(l in 1:length(test_short$location)) {
   dyads$event[i] <- ifelse(test_short[l,j] == test_short[l,k], test_short[l,j], NA) 
  }
}
# almost worked -- M1 and M11 have no event, but others all do -- I THINK it's filled it with 1, and then gone to the next line and found an NA in both so filled with that instead?
dyads[21,]
j <- match(dyads$node_1[21], eles)
k <- match(dyads$node_2[21], eles)
for(l in 1:length(test_short$location)) {
  dyads$event[21] <- ifelse(test_short[l,j] == test_short[l,k], test_short[l,j], NA) 
}

head(dyads)
head(events)
events <- events %>% rename(node_1 = elephant)
inner_join(dyads, events, "node_1")

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

######## ANALYSABLE DATA ########
e$location <- paste(e$gps_s, e$gps_e, sep = '_')
e <- e[c(1:3,14,6:13)]

dyads_all <- data.frame(node_1 = el$id1, node_2 = el$id2,
                        type_1 = el$dem_class1, type_2 = el$dem_class2)
events_all <- data.frame(event = e$encounter, location = e$location, node_1 = e$elephant)

test_all <- full_join(dyads_all, events_all, by = 'node_1')
## not yet sure that this is working correctly - can't currently create the dataset containing only events where all males were identified to check it against, but just see if you can work it from here
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




