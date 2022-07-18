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
library(cmdstanr)
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

boxplot(dyads$edge ~ dyads$pair_type)                                  # plotting some really low values for in-group dyads
points(y = dyads[dyads$pair_type == 'group',]$edge, x = rep(2, length(which(dyads$pair_type == 'group'))),
       col = col.alpha('black',0.2), pch = 19)

for(i in 1:nrow(dyads)){
  if(dyads$pair_type[i] == 'family')  {dyads$edge[i] <- rethinking::rbeta2(1,0.95,50)   # mother-calf: edge weight around 0.95
  } else 
    if(dyads$pair_type[i] == 'group') {dyads$edge[i] <- rethinking::rbeta2(1,0.8,50)    # same family, not mother-calf: edge ~ 0.8
    } else 
      dyads$edge[i] <- rethinking::rbeta2(1,0.20,50)                                        # different family: edge weight around 0.2
}
summary(dyads$edge)

boxplot(dyads$edge ~ dyads$pair_type)                                          # much better
points(y = dyads$edge, x = rep(2, length(which(dyads$pair_type == 'group'))),
       col = col.alpha('black',0.2), pch = 19)

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
  n_dyads = nrow(dyads),         # Number of dyads
  #dyad_id  = dyads$dyad_id,      # Vector of dyad IDs
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

#####
simdat_ls <- list(
  n_dyads  = nrow(dyads),        # Number of dyads
  together = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart         # Number of sightings of each dyad apart
)

### Dan model stancode provided 23rd January
# "simpleBetaNet_DWF_22.01.23.stan"
#'data{
#  int n_dyads;
#  vector[n_dyads] apart;
#  vector[n_dyads] together;
#}
#parameters {
#  vector<lower=0,upper=1>[n_dyads] weight;            // Outcome variable
#}
#model {
#  for (n in 1:n_dyads){
#    weight[n] ~ beta( 1 + together[n], 1 +  apart[n] );
#  }
#}'

# Compile Stan model
set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')
mod <- cmdstan_model("simpleBetaNet_DWF_22.01.23.stan")

# Fit model
edge_weight <- mod$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

# check model
edge_weight # summary, precis, traceplot do not work
edge_weight$output_files()
output1 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_DWF_22.01.23-202201251132-1-11bff9.csv")
output2 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_DWF_22.01.23-202201251132-2-11bff9.csv")
output3 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_DWF_22.01.23-202201251132-3-11bff9.csv")
output4 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_DWF_22.01.23-202201251132-4-11bff9.csv")

draws1 <- as.data.frame(output1$post_warmup_draws)
draws2 <- as.data.frame(output1$post_warmup_draws)
draws3 <- as.data.frame(output1$post_warmup_draws)
draws4 <- as.data.frame(output1$post_warmup_draws)

colnames(draws1[1:10])
colnames(draws2[1:10])
colnames(draws3[1:10])
colnames(draws4[1:10])

draws <- rbind(draws1, draws2, draws3, draws4)

# try again in rethinking because that I can work out how to check the output -- nope, can't make this work
ew2 <- ulam(alist(
  weight ~ beta(together[n_dyads], apart[n_dyads])
), data = simdat_ls, chains = 4, cores = 4)
precis(ew2) # much smaller n_eff
simdat_ls2 <- list(n_dyads  = nrow(dyads),
                   dyad_id  = as.integer(as.factor(dyads$dyad)),
                   together = dyads$event_count,
                   apart    = dyads$apart)
ew2 <- ulam(alist(
  weight ~ dbeta(together[dyad_id], apart[dyad_id])
), data = simdat_ls2, chains = 4, cores = 4)

# give up and try again to work out the Stan model output
plot(draws$`1.weight[1]`, type = 'l',     # looks like a decent traceplot to me
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws$`1.weight[2]`, col = 'tan')
lines(draws$`1.weight[3]`, col = 'orange')
lines(draws$`1.weight[4]`, col = 'green')
lines(draws$`1.weight[5]`, col = 'chocolate')
lines(draws$`1.weight[6]`, col = 'blue')
lines(draws$`1.weight[7]`, col = 'red')
lines(draws$`1.weight[8]`, col = 'seagreen')
lines(draws$`1.weight[9]`, col = 'purple')
lines(draws$`1.weight[10]`,col = 'magenta')
# all look pretty good, all overlap nicely
lines(draws[which(dyads$family_1 == dyads$id_2)[1]+1], col = 'black')      # +1 because first column is not linked to a dyad
lines(draws[which(dyads$family_1 == dyads$id_2)[2]+1], col = 'tan')
lines(draws[which(dyads$family_1 == dyads$id_2)[3]+1], col = 'orange')
lines(draws[which(dyads$family_1 == dyads$id_2)[4]+1], col = 'green')
lines(draws[which(dyads$family_1 == dyads$id_2)[5]+1], col = 'chocolate')
lines(draws[which(dyads$family_1 == dyads$id_2)[6]+1], col = 'blue')
lines(draws[which(dyads$family_1 == dyads$id_2)[7]+1], col = 'red')
lines(draws[which(dyads$family_1 == dyads$id_2)[8]+1], col = 'seagreen')
lines(draws[which(dyads$family_1 == dyads$id_2)[9]+1], col = 'purple')
lines(draws[which(dyads$family_1 == dyads$id_2)[10]+1], col = 'magenta')
# all much higher than non-family, all overlap nicely
which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[1]+1], col = 'black')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[2]+1], col = 'tan')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[3]+1], col = 'orange')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[4]+1], col = 'green')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[5]+1], col = 'chocolate')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[6]+1], col = 'blue')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[7]+1], col = 'red')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[8]+1], col = 'seagreen')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[9]+1], col = 'purple')
lines(draws[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[10]+1], col = 'magenta')
# in between family and non-group, all overlap nicely

which(dyads$family_1 == dyads$id_2) # 5410 5541 5648 5699 5810... 6561 6607 6630 6673 6772
which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2) # 4886 4940 4943 4956 4962... 7117 7125 7126 7128 7132
dens(draws$`1.weight[1]`, xlim = c(0,1), ylim = c(0,10), las = 1, lty = 3)
dens(add = T, lty = 3, draws$`1.weight[2]`, col = 'tan')
dens(add = T, lty = 3, draws$`1.weight[3]`, col = 'orange')
dens(add = T, lty = 3, draws$`1.weight[4]`, col = 'green')
dens(add = T, lty = 3, draws$`1.weight[5]`, col = 'chocolate')
dens(add = T, lty = 3, draws$`1.weight[6]`, col = 'blue')
dens(add = T, lty = 3, draws$`1.weight[7]`, col = 'red')
dens(add = T, lty = 3, draws$`1.weight[8]`, col = 'seagreen')
dens(add = T, lty = 3, draws$`1.weight[9]`, col = 'purple')
dens(add = T, lty = 3, draws$`1.weight[10]`,col = 'magenta')
dens(add = T, lty = 2, draws$`1.weight[5410]`, col = 'black')
dens(add = T, lty = 2, draws$`1.weight[5541]`, col = 'tan')
dens(add = T, lty = 2, draws$`1.weight[5648]`, col = 'orange')
dens(add = T, lty = 2, draws$`1.weight[5699]`, col = 'green')
dens(add = T, lty = 2, draws$`1.weight[5810]`, col = 'chocolate')
dens(add = T, lty = 2, draws$`1.weight[6561]`, col = 'blue')
dens(add = T, lty = 2, draws$`1.weight[6607]`, col = 'red')
dens(add = T, lty = 2, draws$`1.weight[6630]`, col = 'seagreen')
dens(add = T, lty = 2, draws$`1.weight[6673]`, col = 'purple')
dens(add = T, lty = 2, draws$`1.weight[6772]`, col = 'magenta')
dens(add = T, lty = 1, draws$`1.weight[4886]`, col = 'black')
dens(add = T, lty = 1, draws$`1.weight[4940]`, col = 'tan')
dens(add = T, lty = 1, draws$`1.weight[4943]`, col = 'orange')
dens(add = T, lty = 1, draws$`1.weight[4956]`, col = 'green')
dens(add = T, lty = 1, draws$`1.weight[4962]`, col = 'chocolate')
dens(add = T, lty = 1, draws$`1.weight[7117]`, col = 'blue')
dens(add = T, lty = 1, draws$`1.weight[7125]`, col = 'red')
dens(add = T, lty = 1, draws$`1.weight[7126]`, col = 'seagreen')
dens(add = T, lty = 1, draws$`1.weight[7128]`, col = 'purple')
dens(add = T, lty = 1, draws$`1.weight[7132]`, col = 'magenta')


mean(draws[1:1000,2])
mean(draws[1001:2000,2])
mean(draws[2001:3000,2])
mean(draws[3001:4000,2])
mean(draws[1:1000,3])
mean(draws[1001:2000,3])
mean(draws[2001:3000,3])
mean(draws[3001:4000,3])
# all the chains produce exactly the same mean??

means <- data.frame(dyad = dyads$dyad, mean = colMeans(draws[2:7141])#, mean_ch1 = mean(draws[1:1000,2]), mean_ch2 = mean(draws[1001:2000,2:7141]), mean_ch3 = mean(draws[2001:3000,2:7141]), mean_ch4 = mean(draws[3001:4000,2:7141])
                    )
draws_pi <- apply(draws[2:7141], 2, rethinking::PI)
means$pi_lwr <- draws_pi[1,]
means$pi_upr <- draws_pi[2,]

plot_data <- left_join(x = means, y = dyads, by = 'dyad')
plot_data$age_diff_jitter <- jitter(plot_data$age_diff,1)

plot(mean ~ age_diff_jitter, data = plot_data,
     xlab = 'age category difference', ylab = 'edge weight',
     ylim = c(0,1), pch = 19, col = col.alpha(rangi2, 0.2), las = 1)
for(i in 1:nrow(plot_data)){
  lines(x = c(plot_data$age_diff_jitter[i], plot_data$age_diff_jitter[i]),
        y = c(plot_data$pi_lwr[i]  , plot_data$pi_upr[i]  ),
        lty = ifelse(plot_data$sex_diff[i] == 0, 1, 2),
        col = ifelse(plot_data$sex_dyad[i] == 'M_M' | plot_data$sex_dyad[i] == 'F_M' | 
                       plot_data$sex_dyad[i] == 'M_U', 'red', 'green'))
}

#cliques  <- plot_data[plot_data$clique_1 == plot_data$clique_2 & plot_data$id_1 != plot_data$family_2,]
#cliques <- cliques[!is.na(cliques$dyad),]
#families <- plot_data[plot_data$id_1 == plot_data$family_2,] ; families <- families[!is.na(families$dyad),]
#others <- anti_join(plot_data, cliques)

plot_data$pair_type <- ifelse(plot_data$family_1 == plot_data$id_2 | plot_data$family_2 == plot_data$id_1,
                              'family',
                              ifelse(plot_data$clique_1 == plot_data$clique_2,
                                     'clique', 'no_group'))
dyads$pair_type <- ifelse(dyads$family_1 == dyads$id_2 | dyads$family_2 == dyads$id_1,
                          'family',
                          ifelse(dyads$clique_1 == dyads$clique_2,
                                 'clique', 'no_group'))
cliques  <- dyads[dyads$pair_type == 'clique',] ; cliques <- cliques[!is.na(cliques$dyad),]
families <- dyads[dyads$pair_type == 'family',] ; families <- families[!is.na(families$dyad),]
others <- dyads[dyads$pair_type == 'no_group',] ; others <- others[!is.na(others$dyad),]

table(plot_data$pair_type)
boxplot(plot_data$mean ~ plot_data$pair_type, # perfect
        las = 1, xlab = 'dyad relationship', ylab = 'mean edge weight', ylim = c(0,1),
        col = c('turquoise1','orchid','firebrick1')) 
abline(h = mean(cliques$edge),  lty = 2, col = 'blue')     # plot true mean for the group
abline(h = mean(families$edge), lty = 2, col = 'purple3')  # plot true mean for the group
abline(h = mean(others$edge),   lty = 2, col = 'darkred')  # plot true mean for the group
points(x = rnorm(nrow(cliques),  1,0.15), y = cliques$edge,  pch = 4, cex = 0.8, col = 'blue')     # true values
points(x = rnorm(nrow(families), 2,0.15), y = families$edge, pch = 4, cex = 0.8, col = 'purple3')  # true values
points(x = rnorm(nrow(others),   3,0.15), y = others$edge,   pch = 4, cex = 0.8, col = 'darkred')  # true values

### add a and b priors
# Compile Stan model
mod2 <- cmdstan_model("simpleBetaNet_HKM_22.01.25.stan")

simdat_ls <- list(
  n_dyads  = nrow(dyads),        # Number of dyads
  together = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart#,        # Number of sightings of each dyad apart
  #a        = rpois(nrow(dyads), 10),
  #b        = rpois(nrow(dyads), 10)
)

rpois(nrow(dyads), 5) %>% table() # gives a fair number of 0s
rpois(nrow(dyads), 10) %>% table() # gives a fair number of 0s

# Fit model
edge_weight_ab <- mod2$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

a <- rnorm(1000, 10,1) ; min(a) # never reaches a negative number, so how is model failing due to negative shape parameters in beta?

edge_weight_ab$output_files()
read_cmdstan_csv('/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_HKM_22.01.25-202201261044-1-36e648.csv')

# check model
edge_weight_ab # summary, precis, traceplot do not work
edge_weight_ab$output_files()
output_ab <- read_cmdstan_csv('/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_HKM_22.01.25-202201261044-1-36e648.csv')
draws_ab <- as.data.frame(output_ab$post_warmup_draws)

################ 8) Run model on real standardised data -- Binomial model to calculate SRI with uncertainty ################
### create data list -- can contain no NA values in any column, even if column is not specified in model
counts_ls <- list(
  n_dyads  = nrow(counts_df),          # total number of times one or other of the dyad was observed
  together = counts_df$all_events,     # count number of sightings seen together
  apart    = counts_df$apart)          # count number of sightings seen apart
unique(is.na(counts_ls$together))      # FALSE -- no NA values in data
unique(is.na(counts_ls$apart))         # FALSE -- no NA values in data
unique(counts_ls$together)             # 32 levels
unique(counts_ls$apart)                # 81 levels
unique(counts_ls$n_dyads)

### run model
mod

# Fit model (slow)
weight_motnp <- mod$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

# check model
weight_motnp
weight_motnp$summary()
weight_motnp$output_files()
output1 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_DWF_22.01.23-202201251305-1-668382.csv")
output2 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_DWF_22.01.23-202201251305-2-668382.csv")
output3 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_DWF_22.01.23-202201251305-3-668382.csv")
output4 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpCmjRdp/simpleBetaNet_DWF_22.01.23-202201251305-4-668382.csv")

draws1 <- as.data.frame(output1$post_warmup_draws)
draws2 <- as.data.frame(output1$post_warmup_draws)
draws3 <- as.data.frame(output1$post_warmup_draws)
draws4 <- as.data.frame(output1$post_warmup_draws)

colnames(draws1[1:10])
colnames(draws2[1:10])
colnames(draws3[1:10])
colnames(draws4[1:10])

draws <- rbind(draws1, draws2, draws3, draws4)
colnames(draws)[2:106954] <- counts_df$dyad

# build traceplots
plot_cols <- sample(x = 1:106954, size = 30, replace = F)
plot(draws[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws[,plot_cols[2]], col = 'tan')
lines(draws[,plot_cols[3]], col = 'orange')
lines(draws[,plot_cols[4]], col = 'green')
lines(draws[,plot_cols[5]], col = 'chocolate')
lines(draws[,plot_cols[6]], col = 'blue')
lines(draws[,plot_cols[7]], col = 'red')
lines(draws[,plot_cols[8]], col = 'seagreen')
lines(draws[,plot_cols[9]], col = 'purple')
lines(draws[,plot_cols[10]],col = 'magenta')
lines(draws[,plot_cols[11]],col = 'black')
lines(draws[,plot_cols[12]], col = 'tan')
lines(draws[,plot_cols[13]], col = 'orange')
lines(draws[,plot_cols[14]], col = 'green')
lines(draws[,plot_cols[15]], col = 'chocolate')
lines(draws[,plot_cols[16]], col = 'blue')
lines(draws[,plot_cols[17]], col = 'red')
lines(draws[,plot_cols[18]], col = 'seagreen')
lines(draws[,plot_cols[19]], col = 'purple')
lines(draws[,plot_cols[20]],col = 'magenta')
lines(draws[,plot_cols[21]],col = 'black')
lines(draws[,plot_cols[22]], col = 'tan')
lines(draws[,plot_cols[23]], col = 'orange')
lines(draws[,plot_cols[24]], col = 'green')
lines(draws[,plot_cols[25]], col = 'chocolate')
lines(draws[,plot_cols[26]], col = 'blue')
lines(draws[,plot_cols[27]], col = 'red')
lines(draws[,plot_cols[28]], col = 'seagreen')
lines(draws[,plot_cols[29]], col = 'purple')
lines(draws[,plot_cols[30]],col = 'magenta')

# plot for a single individual -- Sierra = F52, herd members = F60+U21+F98, calf = U17
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,4000))
lines(draws$F52_M40,col = 'black')      # non-herd member, adult male
lines(draws$F52_M15, col = 'tan')       # non-herd member, pubescent male
lines(draws$F52_M203, col = 'orange')   # non-herd member, adult male
lines(draws$F52_M26, col = 'green')     # non-herd member, calf
lines(draws$F52_F8, col = 'chocolate')  # non-herd member, adult female
lines(draws$F52_U9, col = 'blue')       # non-herd member, calf
lines(draws$F52_F98, col = 'red')       # herd member most frequently absent from sightings
lines(draws$F52_U17, col = 'purple')    # calf
lines(draws$F52_U21, col = 'seagreen')  # sister
lines(draws$F52_F60,col = 'magenta')    # sister's calf
# IT WORKED!!!!

# plot each chain individually to check mixing
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1$`1.weight[1]`, col = 'purple')
lines(draws2$`1.weight[1]`, col = 'green')
lines(draws3$`1.weight[1]`, col = 'blue')
lines(draws4$`1.weight[1]`, col = 'magenta')
lines(draws1$`1.weight[2]`, col = 'purple')
lines(draws2$`1.weight[2]`, col = 'green')
lines(draws3$`1.weight[2]`, col = 'blue')
lines(draws4$`1.weight[2]`, col = 'magenta')
lines(draws1$`1.weight[3]`, col = 'purple')
lines(draws2$`1.weight[3]`, col = 'green')
lines(draws3$`1.weight[3]`, col = 'blue')
lines(draws4$`1.weight[3]`, col = 'magenta')
which(counts_df$dyad == 'F52_U17') # 41192
lines(draws1$`1.weight[41192]`, col = 'purple')
lines(draws2$`1.weight[41192]`, col = 'green')
lines(draws3$`1.weight[41192]`, col = 'blue')
lines(draws4$`1.weight[41192]`, col = 'magenta')
which(counts_df$dyad == 'F52_F60') # 40896
lines(draws1$`1.weight[40896]`, col = 'purple')
lines(draws2$`1.weight[40896]`, col = 'green')
lines(draws3$`1.weight[40896]`, col = 'blue')
lines(draws4$`1.weight[40896]`, col = 'magenta')
which(counts_df$dyad == 'F52_U21') # 41197
lines(draws1$`1.weight[41197]`, col = 'purple')
lines(draws2$`1.weight[41197]`, col = 'green')
lines(draws3$`1.weight[41197]`, col = 'blue')
lines(draws4$`1.weight[41197]`, col = 'magenta')
which(counts_df$dyad == 'F52_F98') # 40937
lines(draws1$`1.weight[40937]`, col = 'purple')
lines(draws2$`1.weight[40937]`, col = 'green')
lines(draws3$`1.weight[40937]`, col = 'blue')
lines(draws4$`1.weight[40937]`, col = 'magenta')
# OK so every chain has identical outputs, with a perfect match of the draws made -- this doesn't seem right...

# density plots
dens(draws[,plot_cols[1]], ylim = c(0,50), xlim = c(0,1), las = 1, lty = 3)
dens(add = T, draws[,plot_cols[2]], col = 'tan')
dens(add = T, draws[,plot_cols[3]], col = 'orange')
dens(add = T, draws[,plot_cols[4]], col = 'green')
dens(add = T, draws[,plot_cols[5]], col = 'chocolate')
dens(add = T, draws[,plot_cols[6]], col = 'blue')
dens(add = T, draws[,plot_cols[7]], col = 'red')
dens(add = T, draws[,plot_cols[8]], col = 'seagreen')
dens(add = T, draws[,plot_cols[9]], col = 'purple')
dens(add = T, draws[,plot_cols[10]],col = 'magenta')
dens(add = T, draws[,plot_cols[11]],col = 'black')
dens(add = T, draws[,plot_cols[12]], col = 'tan')
dens(add = T, draws[,plot_cols[13]], col = 'orange')
dens(add = T, draws[,plot_cols[14]], col = 'green')
dens(add = T, draws[,plot_cols[15]], col = 'chocolate')
dens(add = T, draws[,plot_cols[16]], col = 'blue')
dens(add = T, draws[,plot_cols[17]], col = 'red')
dens(add = T, draws[,plot_cols[18]], col = 'seagreen')
dens(add = T, draws[,plot_cols[19]], col = 'purple')
dens(add = T, draws[,plot_cols[20]],col = 'magenta')
dens(add = T, draws[,plot_cols[21]],col = 'black')
dens(add = T, draws[,plot_cols[22]], col = 'tan')
dens(add = T, draws[,plot_cols[23]], col = 'orange')
dens(add = T, draws[,plot_cols[24]], col = 'green')
dens(add = T, draws[,plot_cols[25]], col = 'chocolate')
dens(add = T, draws[,plot_cols[26]], col = 'blue')
dens(add = T, draws[,plot_cols[27]], col = 'red')
dens(add = T, draws[,plot_cols[28]], col = 'seagreen')
dens(add = T, draws[,plot_cols[29]], col = 'purple')
dens(add = T, draws[,plot_cols[30]],col = 'magenta')

dens(draws$F52_M40,
     ylim = c(0,75), xlim = c(0,1), las = 1)    # non-herd member, adult male
dens(add = T, draws$F52_M15, col = 'tan')       # non-herd member, pubescent male
dens(add = T, draws$F52_M203, col = 'orange')   # non-herd member, adult male
dens(add = T, draws$F52_M26, col = 'green')     # non-herd member, calf
dens(add = T, draws$F52_F8, col = 'chocolate')  # non-herd member, adult female
dens(add = T, draws$F52_U9, col = 'blue')       # non-herd member, calf
dens(add = T, draws$F52_F98, col = 'red')       # herd member most frequently absent from sightings
dens(add = T, draws$F52_U17, col = 'purple')    # calf
dens(add = T, draws$F52_U21, col = 'seagreen')  # sister
dens(add = T, draws$F52_F60,col = 'magenta')    # sister's calf

means <- data.frame(dyad = counts_df$dyad, mean = colMeans(draws[2:106954]))

draws_pi <- apply(draws[2:106954], 2, rethinking::PI)
means$pi_lwr <- draws_pi[1,]
means$pi_upr <- draws_pi[2,]

plot_data <- left_join(x = means, y = counts_df, by = 'dyad')
plot_data$age_diff_jitter <- jitter(plot_data$age_diff,1)

plot(mean ~ age_diff_jitter, data = plot_data,
     xlab = 'age category difference', ylab = 'edge weight',
     ylim = c(0,1), pch = 19, col = col.alpha(rangi2, 0.2), las = 1)
for(i in 1:nrow(plot_data)){
  lines(x = c(plot_data$age_diff_jitter[i], plot_data$age_diff_jitter[i]),
        y = c(plot_data$pi_lwr[i]  , plot_data$pi_upr[i]  ),
        lty = ifelse(plot_data$sex_diff[i] == 0, 1, 2),
        col = ifelse(plot_data$sex_dyad[i] == 'M_M' | plot_data$sex_dyad[i] == 'F_M' | 
                       plot_data$sex_dyad[i] == 'M_U', 'red', 'green'))
}

par(mar = c(4,6,1,2))
boxplot(plot_data$mean ~ plot_data$dem_type, # perfect
        las = 1, xlab = 'mean edge weight', ylab = '', ylim = c(0,1),
        horizontal = T, outpch = 4, outcex = 0.4, cex.axis = 0.5,
        col = ifelse(plot_data$dem_class_1 == 'AM' | plot_data$dem_class_2 == 'AM', 'blue',
                     ifelse(plot_data$dem_class_1 == 'PM' | plot_data$dem_class_2 == 'PM', 'purple',
                            'magenta')))

unique(plot_data$dem_type)
unique(counts_df$dem_type)

plot_data$dem_class_1_cat <- as.integer(as.factor(plot_data$dem_class_1))
plot_data$dem_class_2_cat <- as.integer(as.factor(plot_data$dem_class_2))
plot_data$dem_type_short <- ifelse(plot_data$dem_class_1_cat <= plot_data$dem_class_2_cat,
                                   paste(plot_data$dem_class_1, plot_data$dem_class_2, sep = '_'),
                                   paste(plot_data$dem_class_2, plot_data$dem_class_1, sep = '_'))
sort(unique(plot_data$dem_type_short))
boxplot(plot_data$mean ~ plot_data$dem_type_short,
        las = 1, xlab = 'mean edge weight', ylab = '', ylim = c(0,1),
        horizontal = T, outpch = 4, outcex = 0.4, cex.axis = 0.5)
plot_data$dem_type_vshort <- ifelse(plot_data$age_type == 'Calf_Calf', 'CC',
                                    ifelse(plot_data$age_type == 'Juvenile_Calf', 'JC',
                                           ifelse(plot_data$age_type == 'Juvenile_Juvenile', 'JJ',
                                                  plot_data$dem_type_short)))
unique(plot_data$dem_type_vshort)
unique(plot_data$age_type)

errors <- plot_data[plot_data$dem_type_vshort == 'CU_CU',]
plot_data$age_class_1 <- ifelse(plot_data$id_1 == 'U8', 'Calf', plot_data$age_class_1)
plot_data$age_class_2 <- ifelse(plot_data$id_2 == 'U8', 'Calf', plot_data$age_class_2)
plot_data$age_type[106953] <- 'Calf_Calf'
plot_data$dem_type_vshort <- ifelse(plot_data$age_type == 'Calf_Calf', 'CC',
                                    ifelse(plot_data$age_type == 'Juvenile_Calf' | plot_data$age_type == 'Calf_Juvenile', 'JC',
                                           ifelse(plot_data$age_type == 'Juvenile_Juvenile', 'JJ',
                                                  plot_data$dem_type_short)))
plot_data$dem_type_vshort <- ifelse(plot_data$age_class_1 == 'Calf',
                                    paste(plot_data$dem_class_2, '_C', sep = ''),
                                    ifelse(plot_data$age_class_2 == 'Calf',
                                           paste(plot_data$dem_class_2, '_C', sep = ''),
                                           plot_data$dem_type_vshort))
plot_data$dem_type_vshort <- ifelse(plot_data$age_class_1 == 'Juvenile',
                                    paste(plot_data$dem_class_2, '_J', sep = ''),
                                    ifelse(plot_data$age_class_2 == 'Juvenile',
                                           paste(plot_data$dem_class_2, '_J', sep = ''),
                                           plot_data$dem_type_vshort))
sort(unique(plot_data$dem_type_vshort))
errors <- plot_data[plot_data$dem_type_vshort == 'CF_J' | plot_data$dem_type_vshort == 'CF_C' | 
                      plot_data$dem_type_vshort == 'CM_C' | plot_data$dem_type_vshort == 'CM_J' | 
                      plot_data$dem_type_vshort == 'CU_C' | plot_data$dem_type_vshort == 'CU_J' |
                      plot_data$dem_type_vshort == 'JF_J' | plot_data$dem_type_vshort == 'JM_J' | 
                      plot_data$dem_type_vshort == 'JU_J' ,]

plot_data$dem_type_vshort2 <- paste(plot_data$dem_type_vshort,' ')
boxplot(plot_data$mean ~ plot_data$dem_type_vshort2,
        las = 1, xlab = 'mean edge weight', ylab = '', ylim = c(0,1),
        horizontal = T, outpch = 4, outcex = 0.4, cex.axis = 0.8)


types <- data.frame(all = sort(unique(plot_data$dem_type_short)),
                    type1 = NA,
                    type2 = NA)
types <- separate(types, col = all, into = c('type1', 'type2'), remove = F)
types <- separate(types, type1, into = c('class1','sex1'), sep = 1, remove = F)
types <- separate(types, type2, into = c('class2','sex2'), sep = 1, remove = F)
head(types)
types$join <- types$all
types$join <- case_when(types$class1 == 'C' & types$class2 == 'C' ~ 'CC',
                        types$class1 == 'J' & types$class2 == 'C' ~ 'JC',
                        types$class1 == 'C' & types$class2 == 'J' ~ 'JC',
                        types$class1 == 'J' & types$class2 == 'J' ~ 'JJ',
                        types$class1 == 'C' & types$class2 == 'P' ~ paste(types$type2,'C',sep='_'),
                        types$class1 == 'P' & types$class2 == 'C' ~ paste(types$type1,'C',sep='_'),
                        types$class1 == 'C' & types$class2 == 'A' ~ paste(types$type2,'C',sep='_'),
                        types$class1 == 'A' & types$class2 == 'C' ~ paste(types$type1,'C',sep='_'),
                        types$class1 == 'J' & types$class2 == 'P' ~ paste(types$type2,'J',sep='_'),
                        types$class1 == 'P' & types$class2 == 'J' ~ paste(types$type1,'J',sep='_'),
                        types$class1 == 'J' & types$class2 == 'A' ~ paste(types$type2,'J',sep='_'),
                        types$class1 == 'A' & types$class2 == 'J' ~ paste(types$type1,'J',sep='_'),
                        TRUE ~ types$all)
unique(types$join)
types <- types[,c(1,8)]
colnames(types)[1] <- 'dem_type_short'

plot_data2 <- left_join(plot_data, types, by = 'dem_type_short')
plot_data2$join
boxplot(plot_data2$mean ~ plot_data2$join,
        las = 1, xlab = 'mean edge weight', ylab = '', ylim = c(0,1),
        horizontal = T, outpch = 4, outcex = 0.4, cex.axis = 0.8,
        col = c('magenta','purple','grey','grey','magenta','purple','grey','blue','purple','purple','purple','blue','purple','grey',
                'grey','grey','grey','grey','magenta','purple','grey','purple','purple','blue','purple','grey','grey'))



draws$`1.lp__`
mean(draws$`1.lp__`)
plot(draws$`1.lp__`, type = 'l')
summary(draws$`1.lp__`)

### Summary plots to send to Dan and Colin -- ggplot it!
random_cols <- sample(1:106953, size = 10, replace = F)
f52_cols <- c(which(counts_df$dyad == 'F52_M40'), which(counts_df$dyad == 'F52_M15'), which(counts_df$dyad == 'F52_M203'), which(counts_df$dyad == 'F52_M26'), which(counts_df$dyad == 'F52_F100'), which(counts_df$dyad == 'F52_F8'), which(counts_df$dyad == 'F52_U9'), which(counts_df$dyad == 'F52_F98'), which(counts_df$dyad == 'F52_U17'), which(counts_df$dyad == 'F52_U21'), which(counts_df$dyad == 'F52_F60'))

summary_data <- draws1[2:106954]
colnames(summary_data) <- counts_df$dyad
summary_data <- summary_data[,c(random_cols, f52_cols)]
summary_data$index <- 1:1000
names <- colnames(summary_data)[1:20]
summary_data <- pivot_longer(summary_data, col = names)
summary_data <- summary_data %>% 
  separate(name, into = c('id_1','id_2'), remove = FALSE)
summary_data$sierra <- ifelse(summary_data$id_1 == 'F52', 'F52','Random')
head(summary_data)

cols <- c('darkred','orange','yellow','green','darkgreen','darkred','orange','yellow','green','darkgreen','seagreen','blue','turquoise','purple','orchid','seagreen','blue','turquoise','purple','orchid')
ggplot(data = summary_data, aes(x = index, y = value, col = name))+
  geom_line()+
  facet_wrap(~sierra)+
  scale_colour_manual(values = cols)+
  scale_x_continuous(name = 'draw number')+
  scale_y_continuous(name = 'edge weight')+
  theme_light()+
  labs(colour = 'dyad')

# plot for a single individual -- Sierra = F52, herd members = F60+U21+F98, calf = U17
dens(draws$F52_M40, ylim = c(0,70), las = 1, ylab = 'density', xlim = c(0,1), xlab = 'edge weight', col = 'orchid')
dens(add = T, draws$F52_M15, col = 'orange')       # non-herd member, pubescent male
dens(add = T, draws$F52_M203, col = 'green')   # non-herd member, adult male
dens(add = T, draws$F52_M26, col = 'darkgreen')     # non-herd member, calf
dens(add = T, draws$F52_F8, col = 'seagreen')  # non-herd member, adult female
dens(add = T, draws$F52_U9, col = 'blue')       # non-herd member, calf
dens(add = T, draws$F52_F98, col = 'yellow')       # herd member most frequently absent from sightings
dens(add = T, draws$F52_U17, col = 'turquoise')    # calf
dens(add = T, draws$F52_U21, col = 'darkred')  # sister
dens(add = T, draws$F52_F60,col = 'purple')    # sister's calf


