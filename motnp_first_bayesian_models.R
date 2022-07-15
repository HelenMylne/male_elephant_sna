## Bayesian analysis of ALERT data
#### Information ####
# Script to analyse association data from Mosi-Oa-Tunya National Park, Zambia.
# Input data produced by processing of raw files ("Raw_ALERT_ElephantDatabase_Youldon210811.xlsx") using script "22.01.06_ALERT_bayesian.R".

# Data collected: 19th May 2016-21st December 2017
# Collected by: Mr David Youldon, Mr Dabwiso Sakala, Miss Helen Mylne and other volunteers/interns/facilitated research students working with ALERT during this time
# Data supplied by: ALERT and Mr David Youldon (11th August 2021) and via Volunteer Encounter (Bex Saunders, 19th October 2021)

#### Set up ####
library(tidyverse)
library(dplyr)
library(rstan)
library(rethinking)
library(igraph)
library(dagitty)

################ Plan ################
### Process for a Bayesian model
# 1) Define DAG -- develop ideas for possible causal models, identify sources of confounding
# 2) Create data list -- standardise variables to make prior definition easier
# 3) Create simulated data set to run model through
# 4) Define likelihood distributions, model, and parameters of interest to estimate (if model is generative, convert to statistical)
# 5) Define priors and hyperpriors
# 6) Plot prior impacts to determine if model is reasonable
# 7) Run model on simulated data to confirm inference methods and understanding of processes. Use a single chain for debugging. --> if transitions are divergent then reparameterise
# 8) Run model on real standardised data using 4 chains
# 9) Examine posterior summary and traceplots
# 10) Sample posterior:
# a) extract samples
# b) calculate parameters of interest (e.g. reverse link functions, identify differences between treatments etc.)
# 11) Predict from posterior:
# a) define sequence of x and fix all other predictors
# b) predict many outcome values from x-sequence
# c) Calculate mean and HPDI/PI of outcomes per x
# d) Plot raw data, overlay posterior predictions
# 12) Compare to other models with WAIC/PSIS. Focus on dWAIC and dSE rather than actual values. Penalty terms indicate number of effective parameters. NOTE: model comparison not selection.
# 13) Convert estimates from parameter scale to outcome scale
# 14) Return to DAGs --> decide which is most likely

################ 1) Draw DAGS ################
### reminding myself how dagitty works
test <- dagitty::dagitty("dag{
                         x [exposure];
                         y [outcome];
                         z [unobserved];
                         y <- x;
                         y <- z}")
coordinates(test) <- list(x = c(X = 0, Y = 1, Z = 2),          # x working -- positions at 0,1,2
                          y = c(X = 0, Y = 1, Z = 0))          # y not working -- positions at 0,1,2
drawdag(test)                                                  # produce DAG diagram

### Binomial model -- 10th January 2022
binom <- dagitty::dagitty("dag{
                         age_dyad [exposure];
                         sex_dyad [exposure];
                         weight [outcome];
                         relatedness_dyad [unobserved];
                         identifiable [unobserved];
                         familiarity <- age_dyad -> weight <- sex_dyad -> familiarity;
                         age_dyad -> others <- sex_dyad <- identifiable;
                         familiarity <- others <- familiarity;
                         weight <- relatedness_dyad -> familiarity;
                         others -> weight <- familiarity
                         }")
coordinates(binom) <- list(x = c(age_dyad = 0, sex_dyad = 2, weight = 1, relatedness = 0.5, others = 1.5, familiarity = 1),
                           y = c(age_dyad = 1, sex_dyad = 1, weight = 1, relatedness = 1.5, others = 1.5, familiarity = 2, identifiable = 0))
drawdag(binom)

binom <- dagitty::dagitty("dag{
                         A [exposure];
                         S [exposure];
                         EW [outcome];
                         R [unobserved];
                         I [unobserved];
                         F <- A -> EW <- S -> F;
                         A -> OD <- S <- I;
                         F <- OD <- F;
                         EW <- R -> F;
                         OD -> EW <- F
                         }")
coordinates(binom) <- list(x = c(A = 0, S = 2, EW = 1, R = 0, OD = 2, `F` = 1),
                           y = c(A = 0, S = 0, EW = 0, R = 1, OD = 1, `F` = 1))
drawdag(binom, radius = 6, cex = 1.6)
drawdag(binom, radius = 6, cex = 1.6, col_labels = c('blue','blue','red','black','black','black')) # colours don't work - all blue. wanted to have outcome in red, exposures in blue, and others in black
drawdag(binom, radius = 6, cex = 1.6, #col_labels = list(A = 'blue', S = 'blue', EW = 'red', R = 'black', OD = 'black', `F` = 'black')) # nope
        col_arrow = c('red','blue','blue','black','black','black','black','black','black','red','blue','blue'))

### Bernoulli model -- 10th January 2022
bern <- dagitty::dagitty("dag{
                         EdgeWeight [outcome];
                         Relatedness [unobserved];
                         Personality [unobserved];
                         SexualState [unobserved];
                         CropRaider [unobserved];
                         ResourceAvailability [unobserved];
                         CropRaider <- OtherDyads <- Relatedness -> Familiarity <- Age_dyad -> EdgeWeight <- Sex_dyad -> Familiarity <- OtherDyads <- Familiarity;
                         Methods -> EdgeWeight <- SexualState <- Age_dyad -> OtherDyads <- GroupSize;
                         EdgeWeight <- Relatedness -> Familiarity -> EdgeWeight <- OtherDyads;
                         EdgeWeight <- CropRaider <- ResourceAvailability <- Relatedness <- Location -> Methods -> NoSightings -> EdgeWeight;
                         GroupSize <- ResourceAvailability <- Season -> GroupSize <- Relatedness;
                         Personality -> EdgeWeight <- NoSightings <- Methods <- Location -> RainFall <- Season;
                         RainFall -> ResourceAvailability
                         }")
coordinates(bern) <- list(x = c(Age_dyad = 4,
                                CropRaider = 2,
                                EdgeWeight = 3,
                                Familiarity = 3,
                                GroupSize = 1,
                                Location = 0,
                                Methods = 1,
                                NoSightings = 2,
                                OtherDyads = 2,
                                Personality = 4,
                                RainFall = 0,
                                Relatedness = 2,
                                ResourceAvailability = 1,
                                Season = 0,
                                Sex_dyad = 3,
                                SexualState = 4),
                          y = c(Age_dyad = 3,
                                CropRaider = 1.5,
                                EdgeWeight = 2,
                                Familiarity = 4,
                                GroupSize = 3,
                                Location = 0,
                                Methods = 0,
                                NoSightings = 0,
                                OtherDyads = 3,
                                Personality = 1,
                                RainFall = 1,
                                Relatedness = 4,
                                ResourceAvailability = 1.5,
                                Season = 2,
                                Sex_dyad = 0,
                                SexualState = 2))
drawdag(bern)

brn2 <- dagitty::dagitty("dag{
                         EW [outcome];
                         R [unobserved];
                         P [unobserved];
                         SS [unobserved];
                         CR [unobserved];
                         RA [unobserved];
                         CR <- OD <- R -> F <- A -> EW <- S -> F <- OD <- F;
                         M -> EW <- SS <- A -> OD <- GS;
                         EW <- R -> F -> EW <- OD;
                         EW <- CR <- RA <- R <- L -> M -> NS -> EW;
                         GS <- RA <- WD -> GS <- R;
                         P -> EW <- NS <- M <- L -> RF <- WD;
                         RF -> RA
                         }")
coordinates(brn2) <- list(x = c(A = 4,
                                CR = 2,
                                EW = 3,
                                `F` = 3,
                                GS = 0,
                                L = 1,
                                M = 2,
                                NS = 3,
                                OD = 2,
                                P = 4,
                                RF = 0,
                                R = 2,
                                RA = 1,
                                WD = 0,
                                S = 4,
                                SS = 4),
                          y = c(A = 2.75,
                                CR = 2,
                                EW = 2,
                                `F` = 3.25,
                                GS = 3,
                                L = 1.5,
                                M = 1.5,
                                NS = 1.5,
                                OD = 2.75,
                                P = 1.5,
                                RF = 1.5,
                                R = 3.25,
                                RA = 2.5,
                                WD = 2,
                                S = 3.25,
                                SS = 2))
col_vect <- c('red','blue','blue','blue','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','red','blue','black','black','black','black','black','black')
drawdag(brn2, col_arrow = col_vect, radius = 4.5)

### Binomial model -- 12th January 2022
binom <- dagitty::dagitty("dag{
                         age_dyad [exposure];
                         sex_dyad [exposure];
                         weight [outcome];
                         relatedness_dyad [unobserved];
                         age_dyad -> weight <- sex_dyad;
                         weight <- relatedness_dyad;
                         }")
dagitty::coordinates(binom) <- list(x = c(age_dyad = 0, sex_dyad = 1, weight = 2, relatedness = 3),
                                    y = c(age_dyad = 0, sex_dyad = 2, weight = 1, relatedness = 0))
drawdag(binom)

binom <- dagitty::dagitty("dag{
                         A [exposure];
                         S [exposure];
                         W [outcome];
                         R [unobserved];
                         A -> W <- S;
                         W <- R;
                         }")
dagitty::coordinates(binom) <- list(x = c(A = 0.5, S = 1, W = 1.0, R = 1.5),
                                    y = c(A = 0.5, S = 2, W = 1.2, R = 0.5))
drawdag(binom, radius = 6, cex = 1.6)

### Bernoulli model -- 12th January 2022
bern <- dagitty::dagitty("dag{
                         EdgeWeight [outcome];
                         Relatedness [unobserved];
                         SexualState [unobserved];
                         ResourceAvailability [unobserved];
                         CropRaider <- Relatedness -> EdgeWeight <- Personality; 
                         Age_dyad -> EdgeWeight <- Sex_dyad;
                         EdgeWeight <- SexualState <- Age_dyad;
                         EdgeWeight <- CropRaider <- ResourceAvailability <- RainFall <- Location -> GroupSize;
                         GroupSize <- ResourceAvailability <- Season -> GroupSize <- Relatedness;
                         Location -> RainFall <- Season
                         }")
coordinates(bern) <- list(x = c(Age_dyad = 4,
                                CropRaider = 2,
                                EdgeWeight = 3,
                                GroupSize = 1,
                                Location = 0,
                                RainFall = 0,
                                Relatedness = 2,
                                ResourceAvailability = 1,
                                Season = 0,
                                Sex_dyad = 3,
                                SexualState = 4),
                          y = c(Age_dyad = 3,
                                CropRaider = 1.5,
                                EdgeWeight = 2,
                                GroupSize = 3,
                                Location = 0,
                                RainFall = 1,
                                Relatedness = 4,
                                ResourceAvailability = 3,
                                Season = 1.5,
                                Sex_dyad = 0,
                                SexualState = 2))
drawdag(bern)

bern <- dagitty::dagitty("dag{
                         W [outcome];
                         R [unobserved];
                         SS [unobserved];
                         RA [unobserved];
                         CR [unobserved]
                         S -> CR <- R -> W; 
                         A -> W <- S;
                         W <- SS <- A;
                         W <- CR <- RA <- L;
                         GS <- RA <- WD -> GS <- R
                         }")
coordinates(bern) <- list(x = c(A = 4,
                                CR = 2,
                                W = 3,
                                GS = 1,
                                L = 0,
                                R = 2,
                                RA = 1,
                                WD = 0,
                                S = 3,
                                SS = 4),
                          y = c(A = 0,
                                CR = 0,
                                W = 0,
                                GS = 1,
                                L = 0,
                                R = 1,
                                RA = 0,
                                WD = 1,
                                S = 1,
                                SS = 1))
drawdag(bern, radius = 4.5)
drawdag(bern, cex = 1.6, radius = 6)

paths(bern, from = exposures(A), to = outcomes(W), Z = list(), limit = 100, directed = FALSE) # don't know why not working


################ 2) Create data lists ################
### import data for aggregated model (binomial)
counts_df <- read_delim('data_processed/motnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',')

# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
str(counts_df)

# create variable for age difference
unique(counts_df$age_category_1) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
counts_df$age_category_1 <- ifelse(counts_df$age_category_1 == '1-2','0-3',counts_df$age_category_1)
counts_df$age_cat_id_1 <- ifelse(counts_df$age_category_1 == '0-3', 1,
                                 ifelse(counts_df$age_category_1 == '3-4', 1,
                                        ifelse(counts_df$age_category_1 == '4-5', 1,
                                               ifelse(counts_df$age_category_1 == '5-6', 2,
                                                      ifelse(counts_df$age_category_1 == '6-7', 2,
                                                             ifelse(counts_df$age_category_1 == '7-8', 2,
                                                                    ifelse(counts_df$age_category_1 == '8-9', 2,
                                                                           ifelse(counts_df$age_category_1 == '9-10', 2,
                                                                                  ifelse(counts_df$age_category_1 == '10-15', 3,
                                                                                         ifelse(counts_df$age_category_1 == '15-19', 4,
                                                                                                ifelse(counts_df$age_category_1 == '20-25', 5,
                                                                                                       ifelse(counts_df$age_category_1 == '20-35', 5,
                                                                                                              ifelse(counts_df$age_category_1 == '25-40', 6,
                                                                                                                     ifelse(counts_df$age_category_1 == '35-50', 6,
                                                                                                                            ifelse(counts_df$age_category_1 == '40+', 7,
                                                                                                                                   ifelse(counts_df$age_category_1 == '50+', 7, counts_df$age_category_1))))))))))))))))
counts_df[is.na(counts_df$age_cat_id_1),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_cat_id_1[which(is.na(counts_df$age_cat_id_1))] <- 1

unique(counts_df$age_category_1[counts_df$age_class_1 == 'Calf'])     # shouldn't include any ages over 4-5
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Juvenile']) # shouldn't include any ages under 5-6

counts_df$age_class_1 <- ifelse(counts_df$age_class_1 == 'Calf' & counts_df$age_category_1 == '6-7', 'Juvenile',
                                ifelse(counts_df$age_class_1 == 'Calf' & counts_df$age_category_1 == '7-8', 'Juvenile',
                                       ifelse(counts_df$age_class_1 == 'Juvenile' & counts_df$age_category_1 == '3-4', 'Calf',
                                              ifelse(counts_df$age_class_1 == 'Juvenile' & counts_df$age_category_1 == '4-5', 'Calf',
                                                     counts_df$age_class_1))))

unique(counts_df$age_category_2) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
counts_df$age_category_2 <- ifelse(counts_df$age_category_2 == '1-2','0-3',counts_df$age_category_2)
counts_df$age_cat_id_2 <- ifelse(counts_df$age_category_2 == '0-3', 1,
                                 ifelse(counts_df$age_category_2 == '3-4', 1,
                                        ifelse(counts_df$age_category_2 == '4-5', 1,
                                               ifelse(counts_df$age_category_2 == '5-6', 2,
                                                      ifelse(counts_df$age_category_2 == '6-7', 2,
                                                             ifelse(counts_df$age_category_2 == '7-8', 2,
                                                                    ifelse(counts_df$age_category_2 == '8-9', 2,
                                                                           ifelse(counts_df$age_category_2 == '9-10', 2,
                                                                                  ifelse(counts_df$age_category_2 == '10-15', 3,
                                                                                         ifelse(counts_df$age_category_2 == '15-19', 4,
                                                                                                ifelse(counts_df$age_category_2 == '20-25', 5,
                                                                                                       ifelse(counts_df$age_category_2 == '20-35', 5,
                                                                                                              ifelse(counts_df$age_category_2 == '25-40', 6,
                                                                                                                     ifelse(counts_df$age_category_2 == '35-50', 6,
                                                                                                                            ifelse(counts_df$age_category_2 == '40+', 7,
                                                                                                                                   ifelse(counts_df$age_category_2 == '50+', 7, counts_df$age_category_2))))))))))))))))
counts_df[is.na(counts_df$age_cat_id_2),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
counts_df$age_cat_id_2[which(is.na(counts_df$age_cat_id_2))] <- 1

str(counts_df)
counts_df$age_cat_id_1 <- as.numeric(counts_df$age_cat_id_1)
counts_df$age_cat_id_2 <- as.numeric(counts_df$age_cat_id_2)

counts_df$age_diff <- abs(counts_df$age_cat_id_1 - counts_df$age_cat_id_2)
summary(counts_df$age_diff) # 0-6

# standardise variables
counts_df$sexf_1 <- as.integer(as.factor(counts_df$sex_1))
counts_df$sexf_2 <- as.integer(as.factor(counts_df$sex_2))
counts_df$dem_class_1 <- as.factor(counts_df$dem_class_1)
counts_df$dem_class_2 <- as.factor(counts_df$dem_class_2)
counts_df$dem_type <- as.factor(paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'))
counts_df$dem_diff <- ifelse(counts_df$dem_class_1 == counts_df$dem_class_2, 1, 2)
counts_df$age_type <- as.factor(paste(counts_df$age_class_1, counts_df$age_class_2, sep = '_'))
counts_df$sex_type <- as.factor(paste(counts_df$sex_1, counts_df$sex_2, sep = '_'))
counts_df$sex_diff <- ifelse(counts_df$sex_1 == counts_df$sex_2, 1, 2)
str(counts_df)

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$all_events

summary(counts_df$apart) # 0-82
counts_df$dyad[which(counts_df$apart == 82)]       # M14_M56
counts_df$all_events[which(counts_df$apart == 82)] # never together
counts_df$count_1[which(counts_df$apart == 82)]    # 36 sightings
counts_df$count_2[which(counts_df$apart == 82)]    # 46 sightings
36+46 # 82 -- correct
never_together <- counts_df[counts_df$all_events == 0,]
unique(never_together$apart == never_together$count_1 + never_together$count_2) # all correct

### create data list for aggregated model (binomial)
counts_ls <- list(
  #N = ...,                           # Number of observations -- unnecessary for aggregated data?
  #L = ...,                           # Number of locations -- unnecessary for aggregated data?
  #location_ids = ...,                # Vector of location IDs corresponding to each observation -- unnecessary for aggregated data?
  #herd_type = df$herd_type_id,       # Males, females or mixed group -- unnecessary for aggregated data?
  #total_count = df$total_elephants,  # Total group size in which elephant was observed -- unnecessary for aggregated data?
  #M = nrow(counts_df),               # Number of dyads -- unnecessary for aggregated data?
  dyad_ids  = counts_df$dyad_id,      # Vector of dyad IDs corresponding to each observation
  all_dyad  = counts_df$count_dyad,   # D = total number of times one or other of the dyad was observed
  event_all = counts_df$all_events,   # Count number of sightings
  events_mo = counts_df$mo_events,    # Count number of Male Only sightings
  events_mx = counts_df$mx_events,    # Count number of Mixed sightings
  events_bh = counts_df$bh_events,    # Count number of Breeding Herd sightings
  dem_dyad  = counts_df$dem_type,     # Dyad demographic pairing
  dem_diff  = counts_df$dem_diff,     # Dyad difference in demography (same or different)
  age_dyad  = counts_df$age_type,     # Dyad age pairing
  age_dyad  = counts_df$age_diff,     # Dyad difference in age (subtracted categories)
  sex_dyad  = counts_df$sex_type,     # Dyad sex pairing
  sex_dyad  = counts_df$sex_diff      # Dyad difference in sex (same or different)
)

################ 4) Define likelihood distributions, model and parameters to estimate ################
# X_ij ∼ Binomial(D_ij, p_ij)
# logit(p_ij) = w_ij + ...

# event ~ Binomial(D, p) --> parameter to determine is p = probability of individuals interacting
# logit(p) = w           --> parameter defining p is w = edge weight between individuals. other additional terms might be home range overlap, but in this case all elephants have an overlap in MOTNP. can't include dyad-level effects in binomial model because that would be one parameter for every data point? 

# if there are no additional parameters to include, can skip out the binomial bit and just put in the weight parameter as y to obtain an output. Prior is then "x", and can be beta rather than binomial:
# weight_ij <- beta(together, apart)

################ 3) Create simulated data set ################
### Population
# 120 individuals: 50 male (10 each of ages 3-7), 50 female (10 each of ages 3-7), 20 unknown (10 each of ages 1-2)
# males no preference for any other individual -- edge weight 0.2
# females in 8 cliques -- edge weight 0.8, otherwise 0.2
# unknowns attached to a single female -- edge weight 0.95

# assign females to family units
cliques_F <- sample(1:8,50,replace=T)
table(cliques_F)

# create data frame of all individuals
population <- data.frame(id = rep(NA,120),
                         num = c(1:50,1:50,1:20),
                         sex = c(rep('M',50), rep('F',50), rep('U',20)),
                         age = c(rep(c(3:7), each = 10), rep(c(3:7), each = 10), rep(c(1:2), each = 10)),
                         clique = c(rep(NA, 50), cliques_F, rep(NA, 20)))
population$id <- paste(population$sex, population$num, sep = '')

# randomly assign calves to females
mothers <- data.frame(id = population$id[101:120],
                      family = sample(population$id[61:100], replace = F, size = 20))
mothers
population <- left_join(x = population, y = mothers, by = 'id')

# need a clique for U elephants and a family for F elephants
females <- population[population$sex == 'F', c(1,5,6)]
unknown <- population[population$sex == 'U', c(1,5,6)]
females$mother <- females$id ; unknown$mother <- unknown$family   # create variable to join by
families <- left_join(x = unknown, y = females, by = 'mother')    # join to give unknowns a clique
unknown <- families[,c(1,4,6)]

families <- left_join(x = females, y = unknown, by = 'mother')    # join to give mothers their calf
females <- families[,c(1,5,2)]

colnames(unknown) <- c('id','mother','clique')
colnames(females) <- c('id','mother','clique')
families <- rbind(females, unknown)                               # put mother and calf data into a single data frame

# join to combine male and female/unknown data together
population <- left_join(x = population, y = families, by = 'id')
population <- population[,c(1,3,4,7,8)]
colnames(population)[c(4,5)] <- c('family','clique')

### Create dyadic data frame
dyads <- data.frame(id_1 = rep(population$id, each = 120), id_2 = rep(population$id, 120))               # create new data frame of all dyads
dyads$same <- ifelse(dyads$id_1 == dyads$id_2, 'yes', 'no') ; dyads <- dyads[dyads$same == 'no', c(1,2)] # remove self-dyads (e.g. M1-M1)
dyads$node_1 <- as.integer(as.factor(dyads$id_1)) ; dyads$node_2 <- as.integer(as.factor(dyads$id_2))    # create factor of nodes
dyads$dyad <- ifelse(dyads$node_1 < dyads$node_2,                   # create variable of dyads
                     paste(dyads$id_1, dyads$id_2, sep = '_'),
                     paste(dyads$id_2, dyads$id_1, sep = '_'))

dyads <- data.frame(dyad = unique(dyads$dyad))                                  # create new data frame of only unique dyads
dyads <- separate(dyads, col = dyad, into = c('id_1','id_2'), remove = FALSE)   # nodes columns
dyads$node_1 <- as.integer(as.factor(dyads$id_1))                               # nodes columns
dyads$node_2 <- as.integer(as.factor(dyads$id_2))                               # nodes columns
dyads$dyad_id <- as.integer(as.factor(dyads$dyad))                              # dyad index variable

population$id_1 <- population$id ; population$id_2 <- population$id  # variable to join information on each individual
dyads <- left_join(x = dyads, y = population, by = 'id_1')           # add information about node 1
dyads <- dyads[,c(1:6,8:11)]
colnames(dyads)[c(3,7:10)] <- c('id_2','sex_1','age_1','family_1','clique_1')

dyads <- left_join(x = dyads, y = population, by = 'id_2')           # add information about node 2
dyads <- dyads[,c(1:10,12:15)]
colnames(dyads)[c(2,11:14)] <- c('id_1','sex_2','age_2','family_2','clique_2')

dyads <- dyads[,c(1:7,11,8,12,9,13,10,14)]

### Assign edge weight to dyad pairs
for(i in 1:nrow(dyads)){
  if(is.na(dyads$clique_1[i]) | is.na(dyads$clique_2[i])) {
    dyads$edge[i] <- rethinking::rbeta2(1,0.20,50)                     # not in a family: edge weight around 0.2
  } else 
    ifelse(dyads$clique_1[i] == dyads$clique_2[i],
           ifelse(dyads$family_1[i] == dyads$id_2[i],
                  dyads$edge[i] <- rethinking::rbeta2(1,0.95,50),      # mother-calf: edge weight around 0.95
                  dyads$edge[i] <- rethinking::rbeta2(1,0.8,50)),      # same family, not mother-calf: edge weight around 0.8
           dyads$edge[i] <- rethinking::rbeta2(1,0.20,50))             # different family: edge weight around 0.2
}
unique(dyads$edge)
summary(dyads$edge)

rm(families, females, mothers, unknown, cliques_F)

### Sample observations
N <- 20
dyads$event_count <- rbinom(nrow(dyads), N, prob = dyads$edge)

### add additional columns (dyadic sex, age and dem_class)
head(dyads)
# assign each node to age category (adult, pubescent, juvenile or calf) from age class
dyads$age_cat_1 <- ifelse(dyads$age_1 == 1, 'C', ifelse(dyads$age_1 == 2, 'J', ifelse(dyads$age_1 < 5, 'P', 'A')))
dyads$age_cat_2 <- ifelse(dyads$age_2 == 1, 'C', ifelse(dyads$age_2 == 2, 'J', ifelse(dyads$age_2 < 5, 'P', 'A')))

# combine age categories into single variable for dyad
dyads$age_dyad <- ifelse(dyads$age_cat_1 == 'A',
                         ifelse(dyads$age_cat_2 == 'A', 'A_A',
                                ifelse(dyads$age_cat_2 == 'P', 'A_P',
                                       ifelse(dyads$age_cat_2 == 'J', 'A_J', 'A_C'))),
                         ifelse(dyads$age_cat_1 == 'P',
                                ifelse(dyads$age_cat_2 == 'A', 'A_P',
                                       ifelse(dyads$age_cat_2 == 'P', 'P_P',
                                              ifelse(dyads$age_cat_2 == 'J', 'P_J', 'P_C'))),
                                ifelse(dyads$age_cat_1 == 'J',
                                       ifelse(dyads$age_cat_2 == 'A', 'A_J',
                                              ifelse(dyads$age_cat_2 == 'P', 'P_J',
                                                     ifelse(dyads$age_cat_2 == 'J', 'J_J', 'J_C'))),
                                       ifelse(dyads$age_cat_2 == 'A', 'A_C',
                                              ifelse(dyads$age_cat_2 == 'P', 'P_C',
                                                     ifelse(dyads$age_cat_2 == 'J', 'J_C', 'C_C'))))))
unique(dyads$age_dyad) # "P_P" "P_A" "A_P" "P_C" "P_J" "A_A" "A_C" "A_J" "C_C" "C_J" "J_C" "J_J"
dyads$age_diff <- abs(dyads$age_1 - dyads$age_2)

# create composite measure of demography per node (e.g. AM = adult male, CU = calf of unknown sex)
dyads$dem_class_1 <- paste(dyads$age_cat_1, dyads$sex_1, sep = '')
dyads$dem_class_2 <- paste(dyads$age_cat_2, dyads$sex_2, sep = '')

# combine demographies into single variable for dyad
dyads$demf_1 <- as.integer(as.factor(dyads$dem_class_1)) ; dyads$demf_2 <- as.integer(as.factor(dyads$dem_class_2))
dyads$dem_dyad <- ifelse(dyads$demf_1 < dyads$demf_2,
                         paste(dyads$dem_class_1, dyads$dem_class_2, sep = '_'),
                         paste(dyads$dem_class_2, dyads$dem_class_1, sep = '_'))
dyads$dem_dyad <- ifelse(dyads$dem_dyad == 'CU_JU', 'JU_CU',        # standardise dem_dyad so always in same order
                         ifelse(dyads$dem_dyad == 'CU_PF', 'PF_CU',
                                ifelse(dyads$dem_dyad == 'CU_PM','PM_CU',
                                       ifelse(dyads$dem_dyad == 'JU_PF','PF_JU',
                                              ifelse(dyads$dem_dyad == 'JU_PM','PM_JU',dyads$dem_dyad)))))
dyads$dem_dyad_f <- as.integer(as.factor(dyads$dem_dyad))           # make index variable

dyads$dem_diff <- ifelse(dyads$dem_class_1 == dyads$dem_class_2, 0, 1) # binary: same or different
dyads$sex_dyad <- paste(dyads$sex_1, dyads$sex_2, sep = '_')           # composite sex variable for dyad
dyads$sex_diff <- ifelse(dyads$sex_1 == dyads$sex_2, 0, 1)             # binary: same or different
dyads <- dyads[,c(1:22,25:29)] # remove demf variables -- not needed

# together vs apart per dyad
summary(dyads$event_count)             # together, all out of 20
dyads$apart <- N - dyads$event_count   # apart = 20-together

### convert to data list -- unnecessary now? weight <- beta(together,apart)
colnames(dyads)
simdat_ls <- list(
  M = nrow(dyads),               # Number of dyads
  dyad_id  = dyads$dyad_id,      # Vector of dyad IDs
  event    = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart#,        # Number of sightings of each dyad apart
  #dem_dyad = dyads$dem_dyad,     # Dyad demographic pairing
  #dem_diff = dyads$dem_diff,     # Dyad difference in demography (same or different)
  #age_dyad = dyads$age_dyad,     # Dyad age pairing
  #age_diff = dyads$age_diff,     # Dyad difference in age (subtracted categories)
  #sex_dyad = dyads$sex_dyad,     # Dyad sex pairing
  #sex_diff = dyads$sex_diff      # Dyad difference in sex (same or different)
)

################ 5-6) Define priors and hyperpriors, run prior predictive checks ################
### calculate weight from beta distribution -- no uncertainty
# weight_ij = beta(together, apart)
dyads$weight <- rbeta(n = nrow(dyads), shape1 = dyads$event_count, shape2 = dyads$apart)

plot(event_count ~ edge, data = dyads,
     las = 1, ylab = 'sightings together', xlab = 'assigned edge weight',
     pch = 16, col = col.alpha(rangi2, 0.2))   # this is good!

plot(weight ~ event_count, data = dyads,
     las = 1, xlab = 'sightings together', ylab = 'calculated edge weight',
     pch = 16, col = col.alpha(rangi2, 0.2))   # yes...

plot(weight ~ edge, data = dyads, xlim = c(0,1), ylim = c(0,1),
     las = 1, xlab = 'assigned edge weight', ylab = 'calculated edge weight',
     pch = 16, col = col.alpha(rangi2, 0.2))   # yay!
abline(h = c(0,0.2,0.4,0.6,0.8,1), v = c(0,0.2,0.4,0.6,0.8,1), col = 'gray80', lty = 3)

# to add uncertainty, include it directly in the model, or use beta() instead of rbeta() so it's defining the distribution, not drawing from it
### calculate weight from beta distribution -- no uncertainty
# weight_ij = beta(together, apart)
dyads$weight_dist <- beta(a = dyads$event_count, b = dyads$apart)
summary(dyads$weight_dist)
dyads$weight_dist[1] # this still just gives a number, but running it repeatedly always gives the same number... way too small and no variation

dens(dyads$weight)
dens(dyads$weight_dist)

cliques <- dyads[which(dyads$clique_1 == dyads$clique_2), c(1:5,11:16,28:30)]
dens(cliques$weight_dist)   # nope
dens(cliques$weight)        # this looks good, it's just not right.... because it's individual values not distributions

### retain rbeta not just beta, put into a proper model instead of just a calculation
# my method would be to draw many samples from beta using rbeta, and then calculate means and CIs from that (probably only really works for the graphing, not for the actual calculations)
for(i in 1:nrow(dyads)){
  dyads$weight_dist[i] <- rbeta(n = 10, shape1 = dyads$event_count[i], shape2 = dyads$apart[i])
}
dyads$weight_dist[1]     # still only reports a single number -- does it just save the last draw?
dens(dyads$weight_dist)  # much better

for(i in 1:nrow(dyads)){
  dyads$weight_mean[i] <- mean(dyads$weight_dist[i])
}

plot(dyads$weight_dist ~ dyads$weight_mean) # only used the single value reported in weight_dist, not the whole distribution

for(i in 1:nrow(dyads)){
  dyads$weight_dist[i] <- list(rbeta(n = 10, shape1 = dyads$event_count[i], shape2 = dyads$apart[i]))
}
dyads$weight_dist[1]       # saves all
for(i in 1:nrow(dyads)){
  dyads$weight_mean[i] <- mean(as.numeric(dyads$weight_dist[i])) # doesn't work
}
summary(dyads$weight_mean)


for(i in 1:nrow(dyads)){
  j <- rbeta(n = 10, shape1 = dyads$event_count[i], shape2 = dyads$apart[i])
  dyads$weight_mean[i] <- mean(j)
  dyads$weight_upr[i]  <- quantile(j,0.945)
  dyads$weight_lwr[i]  <- quantile(j,0.055)
}
summary(dyads$weight_mean)
summary(dyads$weight_upr)
summary(dyads$weight_lwr)

plot(dyads$weight_mean ~ dyads$edge, las = 1, pch = 19, col = col.alpha(rangi2, 0.25),
     xlab = 'assigned edge weight', ylab = 'calculated edge weight', xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, lty = 2)
for(i in 1:nrow(dyads)){  # plot CI lines, but for large lump at the bottom only plot every 10th one so they are vaguely visible
  if(dyads$edge[i] < 0.5 & i%%10 == 0){
    lines(x = c(dyads$edge[i],dyads$edge[i]), y = c(dyads$weight_upr[i],dyads$weight_lwr[i]))  
  } else {
    if(dyads$edge[i] > 0.5) {
      lines(x = c(dyads$edge[i],dyads$edge[i]), y = c(dyads$weight_upr[i],dyads$weight_lwr[i]))
    }
  }
}

draws <- 1000
for(i in 1:nrow(dyads)){
  j <- rbeta(n = draws, shape1 = dyads$event_count[i], shape2 = dyads$apart[i])
  dyads$weight_mean[i] <- mean(j)
  dyads$weight_upr[i]  <- quantile(j,0.945)
  dyads$weight_lwr[i]  <- quantile(j,0.055)
}
summary(dyads$weight_mean)
summary(dyads$weight_upr)
summary(dyads$weight_lwr)

plot(dyads$weight_mean ~ dyads$edge, las = 1, pch = 19, col = col.alpha(rangi2, 0.25),
     xlab = 'assigned edge weight', ylab = 'calculated edge weight', xlim = c(0,1), ylim = c(0,1))
abline(a = 0, b = 1, lty = 2)
for(i in 1:nrow(dyads)){  # plot CI lines, but for large lump at the bottom only plot every 10th one so they are vaguely visible
  if(dyads$edge[i] < 0.5 & i%%10 == 0){
    lines(x = c(dyads$edge[i],dyads$edge[i]), y = c(dyads$weight_upr[i],dyads$weight_lwr[i]))  
  } else {
    if(dyads$edge[i] > 0.5) {
      lines(x = c(dyads$edge[i],dyads$edge[i]), y = c(dyads$weight_upr[i],dyads$weight_lwr[i]))
    }
  }
}

# Dan recommended using beta(), so try putting that into model
?beta()
# weight = beta(together, apart) BUT with additional parameters to avoid a flat prior weight = beta(a+together,b+apart)
# weight ~ beta(a+together, b+apart)
# (a,b) ~ rnorm(5,2) -- 5 was the value Dan suggested but I think he may have pulled it from the air with no particular reason/thought beyond it being positive to make a dome shaped beta prior, rather than a flat line (a=b=0). If anything, I would think that the extremes of the beta distribution should be more likely than the centre, in which case want shape values to be <1, but not sure that would work as entire shape value can't be negative?

simdat_ls <- list(together = dyads$event_count ,
                  apart = dyads$apart,
                  n_dyads = nrow(dyads))

ew1 <- ulam(alist(
  weight ~ dbeta(a+event,b+apart),            # nope - can't handle the idea of a real number + integer....
  a ~ dnorm(5,2),
  b ~ dnorm(5,2)
), data = simdat_ls, cores = 4, chains = 1)

ew1 <- ulam(alist(
  weight ~ dbeta(event,apart)                 # try without additional parameters to start with -- can't handle 0 values?
), data = simdat_ls, cores = 4, chains = 1)

ew1 <- ulam(alist(
  weight ~ dbeta(1+event,1+apart)             # nope - can't do this again
), data = simdat_ls, cores = 4, chains = 1)

simdat_ls$event_1 <- simdat_ls$event+1
simdat_ls$apart_1 <- simdat_ls$apart+1
ew1 <- ulam(alist(
  weight ~ dbeta(event_1,apart_1)                 # this works, but I don't know how just adding 1 to both sides affects the shapes themselves -- would think that it will alter the shapes quite substantially
), data = simdat_ls, cores = 4, chains = 1)
precis(ew1)     # pretty poor, SD = 0, mean = 5.5% = 94.5%, n_eff v low
traceplot(ew1)  # actually doesn't look so bad... 

ew1 <- ulam(alist(
  weight ~ dbeta(a+event_1,b+apart_1),    # still no, even once already had some added
  a ~ dnorm(5,2),
  b ~ dnorm(5,2)
), data = simdat_ls, cores = 4, chains = 1)

ew1 <- ulam(alist(
  weight ~ dbeta(a+event,b+apart),    # doesn't complain about this line
  a ~ dpois(5),                       # parameter can't be integer array
  b ~ dpois(5)                        # parameter can't be integer array
), data = simdat_ls, cores = 4, chains = 1)

simdat_ls$event_01 <- simdat_ls$event+0.1
simdat_ls$apart_01 <- simdat_ls$apart+0.1

ew1 <- ulam(alist(
  weight ~ dbeta(event_01,apart_01)                 # much smaller change forced on model
), data = simdat_ls, cores = 4, chains = 1)
precis(ew1)     # pretty poor, SD = 0, mean = 5.5% = 94.5%, n_eff v low
traceplot(ew1)  # actually doesn't look so bad... 

### 
clique1 <- dyads[dyads$clique_1 == 1 | dyads$clique_2 == 1,]
clique1 <- clique1[!is.na(clique1$clique_1) & !is.na(clique1$clique_2),]
plot(clique1$weight ~ clique1$edge)
simdat_cliques <- list()

ew1 <- ulam(alist(
  weight ~ dbeta(event_01,apart_01)
), data = simdat_ls, cores = 4, chains = 1)
precis(ew1)     # pretty poor, SD = 0, mean = 5.5% = 94.5%, n_eff v low

###
simdat_ls <- list(together = as.numeric(dyads$event_count),
                  apart = as.numeric(dyads$apart))
  
ew1 <- ulam(alist(
  weight ~ dbeta((1+a+together),(1+b+apart)),
  a ~ dnorm(5,2),
  b ~ dnorm(5,2)
), data = simdat_ls, cores = 4, chains = 4) 

ew1 <- ulam(alist(
  weight ~ dbeta((1+together),(1+apart))
), data = simdat_ls, cores = 4, chains = 4) 

simdat_ls$add_1 <- rep(1, nrow(dyads))
ew1 <- ulam(alist(
  weight ~ dbeta((add_1+together),(add_1+apart))
), data = simdat_ls, cores = 4, chains = 4)


################ 8) Run model on real standardised data -- Binomial model to calculate SRI with uncertainty ################
### create data list -- can contain no NA values in any column, even if column is not specified in model
counts_ls <- list(
  #total   = counts_df$count_dyad,     # total number of times one or other of the dyad was observed
  together = counts_df$all_events+0.1, # count number of sightings seen together, adjusted so no 0
  apart    = counts_df$apart+0.1)      # count number of sightings seen apart, adjusted so no 0
unique(is.na(counts_ls$together))      # FALSE -- no NA values in data
unique(is.na(counts_ls$apart))         # FALSE -- no NA values in data
unique(counts_ls$together)             # 32 levels
unique(counts_ls$apart)                # 81 levels

### run model
sri1 <- ulam(alist(
  weight ~ dbeta(together, apart)
), data = counts_ls, chains = 1)
traceplot(sri1)
precis(sri1)

sri1 <- ulam(alist(
  weight ~ dbeta(a + together, b + apart),
  a ~ dnorm(5,2),
  b ~ dnorm(5,2)
), data = counts_ls, chains = 1)
traceplot(sri1)
precis(sri1)

stancode(sri1)


##
################ 9-10) Examine and sample posterior ################
################ 11&13) Predict from posterior, convert to outcome scale ################
################ Jordan Hart, BioRxiv paper example ################
# This example covers fitting an edge weight model to count data (where the count of social events per observation is recorded) with an observation-level location effect, basic model checking and diagnostics, visualising networks with uncertainty, calculating probability distributions over network centrality, and propagating network uncertainty into subsequent analyses.

# First of all we'll load in Rstan for model fitting in Stan, dplyr for handling the data, and igraph for network plotting and computing network centrality. We also load in two custom R files: "simulations.R" to generate synthetic data for this example; and "sampler.R" to allow fitting models with uncertainty over network features respectively.

# simulations.R
metropolis <- function(target, initial, iterations=10000, warmup=2000, thin=4, chain_id=1, refresh=2000) {
  k <- length(initial)
  chain <- matrix(0, iterations + warmup, k)
  ll_obs <- NULL
  pred <- NULL
  chain[1, ] <- initial
  acceptances <- c(1)
  prop_dispersion <- 2.38
  for (i in 2:(iterations + warmup)) {
    current <- chain[i - 1, ]
    # Single Component Adaptive Metropolis (SCAM)
    if (i <= 100) {
      candidate <- rnorm(k, mean=current, sd=5)
    } else {
      prop_var <- prop_dispersion * (current_var + 0.05)
      candidate <- rnorm(k, mean=current, sd=sqrt(prop_var))
    }
    current_lk <- target(current, i)
    candidate_lk <- target(candidate, i)
    A <- exp(candidate_lk - current_lk)
    if (runif(1) < A) { # Accept
      chain[i, ] <- candidate
      acceptances[length(acceptances) + 1] <- 1
    } else {
      acceptances[length(acceptances) + 1] <- 0
      chain[i, ] <- current
    }
    # Update current mean and current var for next iteration.
    if (i == 100) {
      current_mean <- colMeans(chain[1:i, ])
      current_var <- apply(chain[1:i, ], 2, var)
    } else if (i >= 100) {
      prev_mean <- current_mean
      prev_var <- current_var
      current_mean <- (i - 1)/i * prev_mean + (1/i) * chain[i, ]
      current_var <- (i - 2)/(i - 1) * prev_var + prev_mean^2 + 1/(i - 1) * (chain[i, ])^2 - i/(i - 1) * current_mean^2
    }
    # Adjust acceptance rate for batch.
    if (i %% 100 == 0) {
      acc_batch <- mean(acceptances[(i - 99):i])
      if (acc_batch < 0.23) {
        prop_dispersion <- prop_dispersion/exp(sqrt(1/(i/100)))
      } else {
        prop_dispersion <- prop_dispersion * exp(sqrt(1/(i/100)))
      }
    }
    # Print out progress
    if (refresh != 0 && i %% refresh == 0) {
      if (i < warmup) {
        cat(paste0("Chain: ", chain_id, " | Iteration: ", i, "/", warmup + iterations, " (Warmup)\n"))
      } else {
        cat(paste0("Chain: ", chain_id, " | Iteration: ", i, "/", warmup + iterations, " (Sampling)\n"))
      }
    }
  }
  # close(pb)
  cat(paste0("Acceptance Rate: ", mean(acceptances), "\n"))
  return(chain[seq(warmup, iterations + warmup - 1, thin), ])
}

# sampler.R
#simulate_binary <- function() {
## code to prepare `example.4.2` dataset goes here
node_names <- c("Rey", "Leia", "Obi-Wan", "Luke", "C-3PO", "BB-8", "R2-D2", "D-O")
node_types <- c("Lifeform", "Lifeform", "Lifeform", "Lifeform", "Droid", "Droid", "Droid", "Droid")
location_names <- c("A", "B", "C", "D", "E", "F")
n <- 8
node_types_binary <- 1 * (node_types == "Lifeform")
node_types_binary %*% t(node_types_binary)
p <- matrix(runif(n^2, max=0.25), n, n) + 0.75 * node_types_binary %*% t(node_types_binary)
p <- p * upper.tri(p)
d <- matrix(sample(10:50, size=n^2, replace=TRUE), n, n)
d <- d * upper.tri(d)
beta_loc <- rnorm(6, 0, 1)
df <- data.frame(matrix(nrow=0, ncol=6))
colnames(df) <- c("node_1", "node_2", "type_1", "type_2", "event", "location")
for (i in 1:n) {
  for (j in 1:n) {
    if (i < j) {
      for (k in 1:d[i, j]) {
        location_id <- sample(1:6, size=1)
        # At least one of them was visible, did they associate?
        logit_p <- qlogis(p[i, j])
        logit_pn <- logit_p + beta_loc[location_id]
        df[nrow(df) + 1, ] <- c(node_names[i], node_names[j], node_types[i], node_types[j], rbinom(1, 1, plogis(logit_pn)), location_names[location_id])
      }
    }
  }
}

df$node_1 <- factor(df$node_1, levels=node_names)
df$node_2 <- factor(df$node_2, levels=node_names)
df$type_1 <- factor(df$type_1, levels=c("Lifeform", "Droid"))
df$type_2 <- factor(df$type_2, levels=c("Lifeform", "Droid"))
df$location <- factor(df$location, levels=location_names)
df$event <- as.integer(df$event)
list(df=df, p=p)

simulate_count <- function() {
  ## code to prepare `example.4.2` dataset goes here
  node_names <- c("Rey", "Leia", "Obi-Wan", "Luke", "C-3PO", "BB-8", "R2-D2", "D-O")
  node_types <- c("Lifeform", "Lifeform", "Lifeform", "Lifeform", "Droid", "Droid", "Droid", "Droid")
  location_names <- c("A", "B", "C", "D", "E", "F")
  n <- 8
  node_types_binary <- 1 * (node_types == "Lifeform")
  node_types_binary %*% t(node_types_binary)
  p <- matrix(runif(n^2, max=5), n, n) + 5 * node_types_binary %*% t(node_types_binary)
  p <- p * upper.tri(p)
  d <- matrix(sample(1:50, size=n^2, replace=TRUE), n, n)
  d <- d * upper.tri(d)
  beta_loc <- rnorm(6, 0, 1)
  df <- data.frame(matrix(nrow=0, ncol=6))
  colnames(df) <- c("node_1", "node_2", "type_1", "type_2", "event_count", "location")
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        for (k in 1:d[i, j]) {
          location_id <- sample(1:6, size=1)
          # At least one of them was visible, did they associate?
          log_p <- log(p[i, j])
          log_pn <- log_p + beta_loc[location_id]
          df[nrow(df) + 1, ] <- c(node_names[i], node_names[j], node_types[i], node_types[j], rpois(1, exp(log_pn)), location_names[location_id])
        }
      }
    }
  }
  
  df$node_1 <- factor(df$node_1, levels=node_names)
  df$node_2 <- factor(df$node_2, levels=node_names)
  df$type_1 <- factor(df$type_1, levels=c("Lifeform", "Droid"))
  df$type_2 <- factor(df$type_2, levels=c("Lifeform", "Droid"))
  df$event_count <- as.integer(df$event_count)
  df$location <- factor(df$location, levels=location_names)
  list(df=df, p=p)
}

simulate_duration <- function() {
  # Define node names and node types
  node_names <- c("Rey", "Leia", "Obi-Wan", "Luke", "C-3PO", "BB-8", "R2-D2", "D-O")
  node_types <- c("Lifeform", "Lifeform", "Lifeform", "Lifeform", "Droid", "Droid", "Droid", "Droid")
  location_names <- c("A", "B", "C", "D", "E", "F")
  # Duration of each sampling/observation period
  obs_time <- 600
  # Create underlying edge weights, rho.
  n <- 8
  logit_p <- matrix(rnorm(n^2, -4, 1), n, n)
  node_types_binary <- 1 * (node_types == "Lifeform")
  logit_p <- logit_p + 3.0 * (node_types_binary %*% t(node_types_binary))
  logit_p <- logit_p * upper.tri(logit_p)
  # Create right-skewed distribution of mean event times, where max_obs_time is the maximum observation time (and therefore maximum event time).
  # mu <- matrix(rbeta(n^2, 20, 100), n, n) * obs_time
  # mu <- mu * upper.tri(mu)
  loc <- rnorm(6)
  lmbd <- matrix(0.001 * rgamma(n^2, 2, 1), n, n)
  lmbd <- lmbd * upper.tri(lmbd)
  # How to set lmbd and mu when there are nuisance effects? Do the maths!
  # Okay, done the maths, now check and implement it!
  d <- matrix(sample(seq(40, 50), n^2, replace=TRUE), n, n) * obs_time
  d <- d * upper.tri(d)
  obs <- data.frame(node_1=character(), node_2=character(), duration=numeric(), event=numeric(), location=character())
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        num_events <- rpois(1, lmbd[i, j] * d[i, j])
        for (k in 1:num_events) {
          location_id <- sample(1:6, 1)
          # logit_pn <- rho[i, j] + loc[location_id]
          # pn <- plogis(logit_pn)
          mu <- (1/lmbd[i, j]) * plogis(logit_p[i, j] - loc[location_id])
          duration <- round(min(rexp(1, 1/mu), obs_time)) # Above-truncated by maximum observation time
          obs[nrow(obs) + 1, ] <- list(node_names[i], node_names[j], duration, 1, location_names[location_id])
        }
        if (num_events == 0) {
          obs[nrow(obs) + 1, ] <- list(node_names[i], node_names[j], 0, 0)
        }
      }
    }
  }
  obs$node_1 <- factor(obs$node_1, levels=node_names)
  obs$node_2 <- factor(obs$node_2, levels=node_names)
  
  obs_agg <- obs %>%
    group_by(node_1, node_2) %>%
    summarise(total_event_time=sum(duration), num_events=sum(event))
  obs_agg$total_obs_time <- t(d)[lower.tri(d)]
  obs_agg$node_1_type <- factor(node_types[obs_agg$node_1], levels=c("Lifeform", "Droid"))
  obs_agg$node_2_type <- factor(node_types[obs_agg$node_2], levels=c("Lifeform", "Droid"))
  obs$location <- factor(obs$location, levels=location_names)
  list(df_obs=obs, df_obs_agg=obs_agg, mu=mu, lmbd=lmbd)
}

#### Simulating data ####
# Now we will simulate data using the `simulate_binary()` function. The rows of the resulting dataframe describe observations at the dyadic level between nodes. In this dataframe, `event` denotes whether or not an undirected social event was observed in this observation period. The exact definition of observation period will depend on the study, but is commonly a sampling period where at least one of the members of the dyad was observed. This can also be a sampling period where both members of the dyad were observed, and the distinction will affect the interpretation of edge weights. See the paper for further discussion on this. `location` denotes the location at which the observation took place, which may be relevant if location is likely to impact the visibility of social events.
set.seed(1)
data <- simulate_count()
df <- data$df
head(df)

### Preparing the data
# Computationally it's easier to work with dyad IDs rather than pairs of nodes in the statistical model, so we'll map the pairs of nodes to dyad IDs before we put the data into the model. The same is true for the location factor, so we will also map the locations to location IDs. We can add these columns to the dataframe using the following code:
df <- df %>%
  group_by(node_1, node_2) %>%
  mutate(dyad_id=cur_group_id()) %>%
  mutate(location_id=as.integer(location))
head(df)

# It will also be useful later to aggregate the dataframe at the dyad level, assign dyad IDs corresponding to each dyad, and calculate total event counts for each dyad. We can do this using:
df_agg <- df %>%
  group_by(node_1, node_2) %>%
  summarise(event_count_total=sum(event_count), dyad_id=cur_group_id()) %>%
  mutate(node_1_id=as.integer(node_1), node_2_id=as.integer(node_2))
head(df_agg)

### Put it into a list object. The data required by the statistical model is defined in `binary_model.stan`.
model_data <- list(
  N=nrow(df), # Number of observations
  M=nrow(df_agg), # Number of dyads
  L=6, # Number of locations
  dyad_ids=df$dyad_id, # Vector of dyad IDs corresponding to each observation
  location_ids=df$location_id, # Vector of location IDs corresponding to each observation
  event_count=df$event_count # Vector of event counts corresponding to each observation
)

#### Fitting the model ####
# To fit the model, we first must compile it and load it into memory using the function `stan_model()` and providing the filepath to the model. The working directory will need to be set to the directory of the model for this to work properly.
model <- rstan::stan_model('data {
  int<lower=0> N;               // Number of data points
  int<lower=0> M;               // Number of dyads
  int<lower=0> L;               // Number of locations
  int<lower=0> dyad_ids[N];     // Dyad ID corresponding to each data point
  int<lower=0> event_count[N];  // Outcome corresponding to each data point (presence/absence)
  int<lower=0> location_ids[N]; // Location ID corresponding to each data point
}
parameters {
  vector[M] log_p;              // Logit edge weights for each dyad.
  vector[L] beta_loc;           // Parameters for location effects.
  real<lower=0> loc_sigma;      // Hyperparameter for location effect adaptive prior standard deviation.
}
transformed parameters {
  vector[N] log_pn = log_p[dyad_ids] + beta_loc[location_ids]; // Logit probability of a social event for each observation.
}
model {
  # // Main model
  event_count ~ poisson(exp(log_pn));
  # // Adaptive prior over location effects
  beta_loc ~ normal(0, loc_sigma);
  # // Priors
  log_p ~ normal(0, 2.5);
  loc_sigma ~ normal(0, 1);
}
generated quantities {
  int event_pred[N] = poisson_rng(exp(log_pn));
}')
#####
# Compiling the model may take a minute or two, but once this is done, the model can be fit using `sampling()`. The argument `cores` sets the number of CPU cores to be used for fitting the model, if your computer has 4 or more cores, it's worth setting this to 4.
fit <- sampling(model, model_data, cores=4, iter=5000, refresh=500)

#### Model checking ####
# The R-hat values provided by Stan indicate how well the chains have converged, with values very close to 1.00 being ideal. Values diverging from 1.00 indicate that the posterior samples may be very unreliable, and shouldn't be trusted. The chains can be plotted using Rstan's `traceplot` function to verify this visually:
traceplot(fit)

# Good R-hat values don't necessarily indicate that the model is performing well, only that the parameter estimates appear to be robust. To check that the model is performing as it should, a predictive check can be used. A predictive check uses the fitted model to make predictions, and compares those predictions to the observed data. The predictions should indicate that the observed data are concordant with the predictions from the model. There are many ways to perform a predictive check, as data can be summarised in many different ways. For the purposes of this example, we'll use a simple density check where the probability distributions of the aggregated event counts are compared against the predictions from the model. Note that this isn't a guarantee that the model predictions are good, only that the predictions have the same event count distribution as the data. Ideally several predictive checks would be used to check the performance of the model.

# This check uses predictions generated by the Stan model as the quantity `event_pred`, with one set of predictions for each step in the MCMC chain. The predictive check will randomly sample 10 of these steps, compute the event counts for each dyad, and plot the densities against the density of the observed event counts from the data.

# Extract event predictions from the fitted model
event_pred <- rstan::extract(fit)$event_pred
num_iterations <- dim(event_pred)[1]

# Plot the density of the observed event counts
plot(density(df_agg$event_count_total), main="", xlab="Dyadic event counts", xlim=c(0, 600), ylim=c(0, 0.006), frame.plot=FALSE)
# Plot the densities of the predicted event counts, repeat for 10 samples
df_copy <- df
for (i in 1:50) {
  df_copy$event_count <- event_pred[sample(1:num_iterations, size=1), ]
  df_agg_copy <- df_copy %>% 
    group_by(node_1, node_2) %>%
    summarise(event_count_total=sum(event_count))
  lines(density(df_agg_copy$event_count_total), ylim=c(0, 0.007), col="#387780")
}
lines(density(df_agg$event_count_total), lwd=3, main="", xlab="Dyadic event counts", ylim=c(0, 0.007))
# axis(side = 1)

log_p_samples <- extract(fit)$log_p
p_quantiles <- apply(log_p_samples, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
plot(t(data$p)[lower.tri(data$p)], p_quantiles[2, ])
df_comparison <- data.frame(true=log(t(data$p)[lower.tri(data$p)]), est=p_quantiles[2, ], est_lower=p_quantiles[1, ], est_upper=p_quantiles[3, ])
ggplot(df_comparison, aes(x=true, y=est)) +
  geom_point(color="#387780") +
  geom_errorbar(aes(ymin=est_lower, ymax=est_upper)) +
  geom_abline() +
  labs(x="Underlying edge weight", y="Estimated edge weight") +
  coord_cartesian(ylim=c(-1, 3)) +
  theme_classic()
ggsave("true_vs_est.png", dpi=600)

# This plot shows that the observed data falls well within the predicted densities, and the predictions suggest the model has captured the main features of the data well. Now we can be reasonably confident that the model has fit correctly and describes the data well, so we can start to make inferences from the model.

### Extracting edge weights
# The main purpose of this part of the framework is to estimate edge weights of dyads. We can access these using the `logit_p` quantity. This will give a distribution of logit-scale edge weights for each dyad, akin to an edge list. A more useful format for network data is usually adjacency matrices, rather than edge lists, so instead we'll convert the distribution of edge lists to a distribution of adjacency matrices, and store the result in an 8 x 8 x 4000 tensor, as there are 8 nodes and 4000 samples from the posterior. 

log_p_samples <- extract(fit)$log_p

adj_tensor <- array(0, c(8, 8, num_iterations))
for (dyad_id in 1:model_data$M) {
  dyad_row <- df_agg[df_agg$dyad_id == dyad_id, ]
  adj_tensor[dyad_row$node_1_id, dyad_row$node_2_id, ] <- log_p_samples[, dyad_id]
}
adj_tensor[, , 1] # Print the first sample of the posterior distribution over adjacency matrices

# The adjacency matrix above corresponds to a single draw of the posterior adjacency matrices. You'll notice that many of the entries are negative, because the edge weights are on the logit scale. These can be transformed back to the [0, 1] range using the logistic function. If there are no additional effects (such as location in our case), the transformed edge weights will be probabilities and the median will be approximately the same as the simple ratio index for each dyad. However, when additional effects are included, the transformed values can no longer be interpreted as probabilities, though they may be useful for visualisation and analysis purposes. We can logistic transform an adjacency matrix using the logistic function (`plogis()` in base R). This will also map 0 values to 0.5, so it will be necessary to set those values back to zero again. This transformation can be achieved using the following code:

plogis(adj_tensor[, , 1]) * upper.tri(adj_tensor[, , 1])

# It will be necessary to use this transformation for the visualisations and analyses we have planned, so we'll apply the transformation to the entire tensor:

adj_tensor_transformed <- adj_tensor
for (i in 1:dim(adj_tensor)[3]) {
  adj_tensor_transformed[, , i] <- plogis(adj_tensor[, , i]) * upper.tri(adj_tensor[, , i])
}

### Visualising uncertainty
# The aim of our network visualisation is to plot a network where the certainty in edge weights (edge weights) can be seen. To do this we'll use a semi-transparent line around each edge with a width that corresponds to a standardised uncertainty measures. The uncertainty measure will simply be the normalised difference between the 97.5% and 2.5% credible interval estimate for each edge weight. We can calculate this from the transformed adjacency tensor object, generate two igraph objects for the main network and the uncertainty in edges, and plot them with the same coordinates.
minmax_norm <- function(x) {
  (x - min(x))/(max(x) - min(x))
}

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles <- apply(adj_tensor_transformed, c(1, 2), function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
adj_lower <- adj_quantiles[1, , ]
adj_mid <- adj_quantiles[2, , ]
adj_upper <- adj_quantiles[3, , ]

# Calculate standardised width/range of credible intervals.
adj_range <- ((adj_upper - adj_lower))
adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one form the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid * (adj_mid > 0.65), mode="undirected", weighted=TRUE)
g_range <- graph_from_adjacency_matrix(adj_range * (adj_mid > 0.65), mode="undirected", weighted=TRUE)

# Plot the median graph first and then the standardised width graph to show uncertainty over edges.
coords <- igraph::layout_nicely(g_mid)
plot(g_mid, edge.width=10 * minmax_norm(E(g_mid)$weight), edge.color="black",  layout=coords)
plot(g_mid, edge.width=60 * minmax_norm(E(g_range)$weight), edge.color=rgb(0, 0, 0, 0.25), vertex.color="#387780",
     # vertex.label=c("1", "2", "3", "4", "5", "6", "7", "D-O"),
     vertex.label.dist=0, vertex.label.color="white", vertex.label.cex=2.5, vertex.label.family="Helvetica", layout=coords, add=TRUE)


# This plot can be extended in multiple ways, for example by thresholding low edge weights to visualise the network more tidily, or by adding halos around nodes to show uncertainty around network centrality, and so on.

### Extracting network centralities
# Uncertainty around network metrics such as centrality can be calculated quite simply by drawing adjacency matrices from the posterior distribution over adjacency matrices, generating a network from them, and calculating the network metric of interest. It is important to sample over the adjacency matrices rather than by the edges on their own, as this maintains the joint distribution of edge weights and will generate more reliable and accurate estimates of network centrality.

centrality_matrix <- matrix(0, nrow=num_iterations, ncol=8)
for (i in 1:num_iterations) {
  g <- graph_from_adjacency_matrix(adj_tensor[, , i], mode="undirected", weighted=TRUE)
  centrality_matrix[i, ] <- eigen_centrality(g)$vector
}
colnames(centrality_matrix) <- c("Rey", "Leia", "Obi-Wan", "Luke", "C-3PO", "BB-8", "R2-D2", "D-O")
head(centrality_matrix) # Each column in this matrix corresponds to one of the nodes in the network, and each row is its centrality in one sample of the posterior distribution of the adjacency matrices. We can calculate the credible intervals using the `quantile` function as follows:

centrality_quantiles <- t(apply(centrality_matrix, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
centrality_quantiles

plot(centrality_quantiles[, 2], 1:nrow(centrality_quantiles), xlim=c(min(centrality_quantiles), max(centrality_quantiles)), xlab="Node strength", ylab="Node ID", col="#387780", )
arrows(x0=centrality_quantiles[, 3], y0=1:nrow(centrality_quantiles), x1=centrality_quantiles[, 1], y1=1:nrow(centrality_quantiles), lwd=3, code=0)

library(ggplot2)
library(tidyr)
df_wide <- data.frame(centrality_matrix)
colnames(df_wide) <- 1:8
df_long <- pivot_longer(df_wide, cols=1:8, names_to="node_id", values_to="Centrality")
ggplot(df_long, aes(x=Centrality)) +
  geom_density(fill="#387780", alpha=0.7, size=0.8) +
  facet_grid(rows=vars(as.factor(node_id)), scales="free") +
  labs(x="Eigenvector centrality") + 
  theme_void() + 
  theme(strip.text.y=element_text(size=30), axis.text.x = element_text(angle = 0, size=30, debug = FALSE), axis.title.x=element_text(size=30), plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

### Maintaining uncertainty in regression on centralities
# The key challenge to quantifying uncertainty in network analysis is to incorporate uncertainty due to sampling into downstream analyses, commonly regression. This can be achieved by modifying the likelihood function of a regression model to treat the network centralities with uncertainty. We have written a custom MCMC sampler function that samples from the joint distribution of network centralities calculated earlier and treats those samples as the data in the likelihood function. Likelihood functions for the sampler use the `index` variable to keep track of which data points are being compared internally in the sampler, to ensure that candidate steps in the MCMC are not accepted or rejected because they are being compared to different data points, rather than because of the parameter space. 

# Custom likelihood functions take a similar form to the `target +=` syntax in Stan, but for more specific resources the following document is a good start: https://www.ime.unicamp.br/~cnaber/optim_1.pdf. We will implement a linear regression to test if lifeforms are more central in the social network than droids. We have included a coefficient for both lifeform and droid, unlike standard frequentist models. This is because using a reference category (such as droid) would imply that there is less uncertainty around the centrality of droids than around lifeforms. It also allows for easy comparison between categories by calculating the difference in posteriors. 

loglik <- function(params, Y, X, index) {
  # Define parameters
  # intercept <- params[1]
  beta_lifeform <- params[1]
  beta_droid <- params[2]
  sigma <- exp(params[3]) # Exponential keeps underlying value unconstrained, which is much easier for the sampler.
  # Sample data according to index
  y <- Y[index %% dim(Y)[1] + 1, ]
  # Define model
  target <- 0
  target <- target + sum(dnorm(y, mean=beta_lifeform * X[, 1] + beta_droid * X[, 2], sd=sigma, log=TRUE)) # Main model
  # target <- target + dnorm(intercept, mean=0, sd=1, log=TRUE) # Prior on intercept
  target <- target + dnorm(beta_lifeform, mean=0, sd=1, log=TRUE) # Prior on lifeform coefficient
  target <- target + dnorm(beta_droid, mean=0, sd=1, log=TRUE) # Prior on droid coefficient
  target <- target + dexp(sigma, 1, log=TRUE) # Prior on sigma
  return(target)
}

# Now we will prepare data for fitting the model. The predictor matrix is simply a matrix with 2 columns and 8 rows, corresponding to whether each of the 8 nodes is a lifeform (column 1) or a droid (column 2).

predictor_matrix <- matrix(0, nrow=8, ncol=2)
colnames(predictor_matrix) <- c("lifeform", "droid")
predictor_matrix[1:4, 1] <- 1
predictor_matrix[5:8, 2] <- 1
predictor_matrix

# Since network strength is strictly positive, a Gaussian error is not a reasonable model for the data. The Gaussian family model is much easier to implement as well as interpret than many other models, so we will standardise the centralities by taking z-scores.

centrality_matrix_std <- (centrality_matrix - apply(centrality_matrix, 1, mean))/apply(centrality_matrix, 1, sd)
centrality_matrix_std[is.nan(centrality_matrix_std)] <-0
head(centrality_matrix_std)

# Now we're in a position to fit the model. To do this, we define the target function, which is simply a function that maps candidate parameters and a network centrality index to the log-likelihood of that function for the given sample of the centrality posterior. This means the target function can be written as a function of the data `centrality_matrix_std` and `predictor_matrix`.

target <- function(params, index) loglik(params, centrality_matrix_std, predictor_matrix, index)

# The function `metropolis` from `sampler.R` can now be used to fit the model using the provided target function, an initial set of parameters, and some additional MCMC options.

chain <- metropolis(target, c(0, 0, 0), iterations=200000, thin=100, refresh=10000)
colnames(chain) <- c("beta_lifeform", "beta_droid", "sigma")
head(chain)

### Checking the regression
# The resulting chain of MCMC samples forms the posterior distribution of parameter estimates for the regression model. But before we look at these too closely, we should check that the chains have converged:
par(mfrow=c(3, 1))
for (i in 1:3) {
  plot(chain[, i], type="l")
}
# These chains appear to be quite healthy. Ideally we would run multiple additional chains starting at different points to check that they converge and mix properly. For the sake of this example we won't go into that here.
# Again, the performance of the sampler doesn't necessarily guarantee the performance of the model, so we'll use predictive checks to test the performance of the model. In this case, the data are not fixed, and there are multiple possible values they can take. Therefore we'll plot the distribution of centrality values on different draws of the adjacency matrices as well as the distribution of predicted centrality values on different draws.
plot(density(centrality_matrix_std[1, ]), ylim=c(0, 0.7), main="", xlab="Standardised node strength", cex.lab=2, frame.plot=FALSE)
sample_ids <- sample(1:1000, size=100)
preds <- sapply(sample_ids, function(i) rnorm(8, mean=chain[i, "beta_lifeform"] * predictor_matrix[, 1] + chain[i, "beta_droid"] * predictor_matrix[, 2], sd=exp(chain[i, "sigma"])))
for (i in 1:length(sample_ids)) {
  pred <- preds[, i]
  lines(density(pred), col=rgb(56/255, 119/255, 128/255, 0.5))
}
for (i in 1:length(sample_ids)) {
  lines(density(centrality_matrix_std[sample_ids[i], ]), col=rgb(0, 0, 0, 0.25))
}
# The model appears to fit reasonably well, and the observed data are completely consistent with the predictions of the model, so we can go ahead with the analysis.

### Interpreting the regression
# The regression coefficients and parameters can be summarised by calculating their percentile credible interval similar to before:

coefficient_quantiles <- t(apply(chain, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
coefficient_quantiles

# A frequentist analysis (and some Bayesian ones too) would have only one category, lifeform or droid, and the other category would be the implicit reference category, absorbed by the intercept. In this type of analysis, the coefficients for the two categories correspond to the average difference between the centrality of nodes in that category compared to the population average (the intercept). Therefore, to look for a difference between the two categories, we can simply calculate the difference in the posterior distributions of those two categories:

beta_difference <- chain[, "beta_lifeform"] - chain[, "beta_droid"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

# The mass of probability is with there being a positive difference of around 1.57 standard deviations between the centralities of lifeforms compared to droids. Many of the benefits of Bayesian analysis only apply when significance testing is avoided. Though it is reasonably common for a result such as the one above not overlapping zero to be interpreted as being "significant", using such a decision rule leaves Bayesian analysis open to the same flaws as frequentist analyses often have. For this reason we caution strongly against using such a rule.

d <- density(beta_difference)
plot(d, lwd=2, xlab="Posterior difference in centrality between sexes", main="", cex.lab=2, frame.plot=FALSE, xlim=c(0, 3))
polygon(d, col=rgb(56/255, 119/255, 128/255, 0.75))

