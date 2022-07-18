#### Bayesian analysis of ALERT data ####
# Script to process association data from Mosi-Oa-Tunya National Park, Zambia.
# Data collected: 19th May 2016-21st December 2017
# Collected by: Mr David Youldon, Mr Dabwiso Sakala, Miss Helen Mylne and other volunteers/interns/facilitated research students working with ALERT during this time
# Data supplied by: ALERT and Mr David Youldon (11th August 2021) and via Volunteer Encounter (Bex Saunders, 19th October 2021)
# Data input: raw data provided by ALERT processed using script 22.01.13_ALERT_bayesian.R

#### Set up ####
# load packages
library(tidyverse)
library(dplyr)
library(rstan)
library(rethinking)
library(igraph)
library(dagitty)
library(cmdstanr)

# information
sessionInfo()
R.Version()
rstan::stan_version()
#packageVersion('')
#citation('')

# set stan path
set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')

# set seed
set.seed(12345)

################ 1) Draw DAGS ################
# plot with full names
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

# plot with letters
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

# clear environment and reset plot window
rm(binom)
dev.off()
################ 2) Create data lists ################
### import data for aggregated model (binomial)
counts_df <- read_delim('data_processed/motnp_bayesian_allpairwiseevents_splitbygrouptype_22.01.13.csv', delim = ',')

# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
counts_df$sex_1 <- as.character(counts_df$sex_1)
str(counts_df)  # sex_1 still comes up as logical?

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

unique(counts_df$age_category_1[counts_df$age_class_1 == 'Calf'])      # shouldn't include any ages over 4-5
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(counts_df$age_category_1[counts_df$age_class_1 == 'Adult'])     # shouldn't include any ages under 20-25

counts_df$age_class_1 <- ifelse(counts_df$age_cat_id_1 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_1 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_1 > 4, 'Adult','Pubescent')))

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

counts_df$age_class_2 <- ifelse(counts_df$age_cat_id_2 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_2 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_2 > 4, 'Adult','Pubescent')))

unique(counts_df$age_category_2[counts_df$age_class_2 == 'Calf'])      # shouldn't include any ages over 4-5
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(counts_df$age_category_2[counts_df$age_class_2 == 'Adult'])     # shouldn't include any ages under 20-25

### correct dem_class with corrected age classes
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))

### correct dem_type with new dem_class
counts_df$age_class_id_1 <- ifelse(counts_df$age_class_1 == 'Adult',4,
                                   ifelse(counts_df$age_class_1 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_1 == 'Juvenile',2,1)))
counts_df$age_class_id_2 <- ifelse(counts_df$age_class_2 == 'Adult',4,
                                   ifelse(counts_df$age_class_2 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_2 == 'Juvenile',2,1)))
counts_df$dem_type <- ifelse(counts_df$age_class_id_1 > counts_df$age_class_id_2,
                             paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'),
                             ifelse(counts_df$age_class_id_1 < counts_df$age_class_id_2,
                                    paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_'),
                                    paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_')))
sort(unique(counts_df$dem_type))

### add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$all_events

### create data list -- can contain no NA values in any column, even if column is not specified in model
counts_ls <- list(
  n_dyads  = nrow(counts_df),          # total number of times one or other of the dyad was observed
  together = counts_df$all_events,     # count number of sightings seen together
  apart    = counts_df$apart)          # count number of sightings seen apart

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
dyads$node_1 <- as.integer(as.factor(dyads$id_1)) ; dyads$node_2 <- as.integer(as.factor(dyads$id_2))    # create factor of  <- 
dyads$dyad <- ifelse(dyads$node_1 < dyads$node_2,                               # create variable of dyads
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

### Assign type of dyad
dyads$pair_type <- ifelse(is.na(dyads$clique_1) | is.na(dyads$clique_2), 'no_group',
                          ifelse(dyads$clique_1 == dyads$clique_2,
                                 ifelse(is.na(dyads$family_1) | is.na(dyads$family_2), 'group',
                                        ifelse(dyads$family_1 == dyads$id_2 | dyads$family_2 == dyads$id_1,
                                               'family', 'group')),
                                 'no_group'))

### Assign edge weight to dyad pairs
for(i in 1:nrow(dyads)){
  if(dyads$pair_type[i] == 'family')  {dyads$edge[i] <- rethinking::rbeta2(1,0.95,50)   # mother-calf: edge weight around 0.95
  } else 
    if(dyads$pair_type[i] == 'group') {dyads$edge[i] <- rethinking::rbeta2(1,0.8,50)    # same family, not mother-calf: edge ~ 0.8
    } else 
      dyads$edge[i] <- rethinking::rbeta2(1,0.20,50)                                        # different family: edge weight around 0.2
}
summary(dyads$edge)

boxplot(dyads$edge ~ dyads$pair_type)                                                   # check simulations are around what I expected
points(y = dyads$edge, x = rep(2, length(which(dyads$pair_type == 'group'))),           # had a few problems with this one, but looks good now
       col = col.alpha('black',0.2), pch = 19)

### clean up environment
rm(families, females, mothers, unknown, cliques_F)

### Sample observations
N <- 20
dyads$event_count <- rbinom(nrow(dyads), N, prob = dyads$edge)

### plot simulated observations against assigned edge weight to check model input is correlates to true value
plot(event_count ~ edge, data = dyads,
     las = 1, ylab = 'sightings together', xlab = 'assigned edge weight',
     pch = 16, col = col.alpha(rangi2, 0.2))

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
unique(dyads$age_dyad) # "P_P" "A_P" "P_C" "P_J" "A_A" "A_C" "A_J" "C_C" "J_C" "J_J"
dyads$age_diff <- abs(dyads$age_1 - dyads$age_2)

# create composite measure of demography per node (e.g. AM = adult male, CU = calf of unknown sex)
dyads$dem_class_1 <- paste(dyads$age_cat_1, dyads$sex_1, sep = '')
dyads$dem_class_2 <- paste(dyads$age_cat_2, dyads$sex_2, sep = '')

# combine demographies into single variable for dyad
dyads$demf_1 <- as.integer(as.factor(dyads$dem_class_1)) ; dyads$demf_2 <- as.integer(as.factor(dyads$dem_class_2))
dyads$dem_dyad <- ifelse(dyads$demf_1 < dyads$demf_2,
                         paste(dyads$dem_class_1, dyads$dem_class_2, sep = '_'),
                         paste(dyads$dem_class_2, dyads$dem_class_1, sep = '_'))
dyads$dem_dyad <- ifelse(dyads$dem_dyad == 'CU_JU', 'JU_CU',           # standardise dem_dyad so always in same order
                         ifelse(dyads$dem_dyad == 'CU_PF', 'PF_CU',
                                ifelse(dyads$dem_dyad == 'CU_PM','PM_CU',
                                       ifelse(dyads$dem_dyad == 'JU_PF','PF_JU',
                                              ifelse(dyads$dem_dyad == 'JU_PM','PM_JU',dyads$dem_dyad)))))
dyads$dem_dyad_f <- as.integer(as.factor(dyads$dem_dyad))              # make index variable

dyads$dem_diff <- ifelse(dyads$dem_class_1 == dyads$dem_class_2, 0, 1) # binary: same or different
dyads$sex_dyad <- paste(dyads$sex_1, dyads$sex_2, sep = '_')           # composite sex variable for dyad
dyads$sex_diff <- ifelse(dyads$sex_1 == dyads$sex_2, 0, 1)             # binary: same or different
dyads <- dyads[,c(1:23,26:30)]                                         # remove demf variables -- not needed

# together vs apart per dyad
summary(dyads$event_count)             # together, all out of 20
dyads$apart <- N - dyads$event_count   # apart = 20-together

### create data list -- can contain no NA values in any column, even if column is not specified in model
simdat_ls <- list(
  n_dyads  = nrow(dyads),       # total number of times one or other of the dyad was observed
  together = dyads$event_count, # count number of sightings seen together
  apart    = dyads$apart)       # count number of sightings seen apart

# clear environment
rm(population, i, N)

################ 4) Define likelihood distributions, model and parameters to estimate ################
# Binomial model using a beta distribution where shape 1 = times together and shape 2 = times apart

### Compile Stan model
mod_1.1 <- cmdstan_model("models/simpleBetaNet_DWF_22.01.23.stan")
mod_1.1

################ 5) Run model on simulated data ################
# create data list
simdat_ls <- list(
  n_dyads  = nrow(dyads),        # Number of dyads
  together = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart         # Number of sightings of each dyad apart
)

### Fit model
edge_weight_1.1 <- mod_1.1$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
edge_weight_1.1$summary()
#   variable        mean      median       sd     mad          q5        q95     rhat ess_bulk ess_tail
#   <chr>          <dbl>       <dbl>    <dbl>   <dbl>       <dbl>      <dbl>    <dbl>    <dbl>    <dbl>
# 1 lp__         -83508.      -83508     61.5    60.8      -83608.    -83407.    1.00     1347.    1966.
# 2 weight[1]      0.183       0.173   0.0806  0.0791      0.0680       0.330    1.00     8950.    2451.
# 3 weight[2]      0.272       0.267   0.0914  0.0896       0.134       0.434    1.00     9029.    2347.
# 4 weight[3]      0.272       0.265   0.0928  0.0974       0.133       0.436    1.00     8872.    2522.
# 5 weight[4]      0.408       0.406   0.0987   0.102       0.251       0.572   0.999    10311.    2992.
# 6 weight[5]     0.0457      0.0322   0.0450  0.0325     0.00221       0.138    1.00     5098.    2204.
# 7 weight[6]      0.137       0.125   0.0727  0.0689      0.0393       0.276    1.00     8291.    2548.
# 8 weight[7]      0.179       0.169   0.0778  0.0760      0.0669       0.325    1.01     9580.    3062.
# 9 weight[8]      0.272       0.266   0.0932  0.0934       0.132       0.437    1.00     8875.    2823.
#10 weight[9]     0.0446      0.0314   0.0434  0.0326     0.00219       0.132    1.00     4795.    1798.
# … with 7,131 more rows
output1 <- read_cmdstan_csv(edge_weight_1.1$output_files()[1])
output2 <- read_cmdstan_csv(edge_weight_1.1$output_files()[2])
output3 <- read_cmdstan_csv(edge_weight_1.1$output_files()[3])
output4 <- read_cmdstan_csv(edge_weight_1.1$output_files()[4])

draws1_1.1 <- as.data.frame(output1$post_warmup_draws)
draws2_1.1 <- as.data.frame(output2$post_warmup_draws)
draws3_1.1 <- as.data.frame(output3$post_warmup_draws)
draws4_1.1 <- as.data.frame(output4$post_warmup_draws)
draws_simdat_1.1 <- rbind(draws1_1.1, draws2_1.1, draws3_1.1, draws4_1.1)
rm(output1, output2, output3, output4)

# check chain mixing -- well mixed
mean(draws1_1.1$`1.weight[1]`) ; mean(draws2_1.1$`1.weight[1]`) ; mean(draws3_1.1$`1.weight[1]`) ; mean(draws4_1.1$`1.weight[1]`)
mean(draws1_1.1$`1.weight[2]`) ; mean(draws2_1.1$`1.weight[2]`) ; mean(draws3_1.1$`1.weight[2]`) ; mean(draws4_1.1$`1.weight[2]`)
mean(draws1_1.1$`1.weight[3]`) ; mean(draws2_1.1$`1.weight[3]`) ; mean(draws3_1.1$`1.weight[3]`) ; mean(draws4_1.1$`1.weight[3]`)

# plot each chain individually to check mixing -- all mixed well
plot(draws1_1.1$`1.weight[1]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_1.1$`1.weight[1]`, col = 'green')
lines(draws3_1.1$`1.weight[1]`, col = 'blue')
lines(draws4_1.1$`1.weight[1]`, col = 'magenta')
plot(draws1_1.1$`1.weight[2]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_1.1$`1.weight[2]`, col = 'green')
lines(draws3_1.1$`1.weight[2]`, col = 'blue')
lines(draws4_1.1$`1.weight[2]`, col = 'magenta')
plot(draws1_1.1$`1.weight[3]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_1.1$`1.weight[3]`, col = 'green')
lines(draws3_1.1$`1.weight[3]`, col = 'blue')
lines(draws4_1.1$`1.weight[3]`, col = 'magenta')

# build traceplots -- highly variable, but generally fits the inputs
plot(draws1_1.1$`1.weight[1]`, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1_1.1$`1.weight[2]`, col = 'tan')
lines(draws1_1.1$`1.weight[3]`, col = 'orange')
lines(draws1_1.1$`1.weight[4]`, col = 'green')
lines(draws1_1.1$`1.weight[5]`, col = 'chocolate')
lines(draws1_1.1$`1.weight[6]`, col = 'blue')
lines(draws1_1.1$`1.weight[7]`, col = 'red')
lines(draws1_1.1$`1.weight[8]`, col = 'seagreen')
lines(draws1_1.1$`1.weight[9]`, col = 'purple')
lines(draws1_1.1$`1.weight[10]`,col = 'magenta')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[1]+1], col = 'black')       # traceplot of mother-calf relationships
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[2]+1], col = 'tan')         # +1 because first column is not linked to a dyad
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[3]+1], col = 'orange')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[4]+1], col = 'green')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[6]+1], col = 'blue')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[7]+1], col = 'red')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[9]+1], col = 'purple')
lines(draws1_1.1[which(dyads$family_1 == dyads$id_2)[10]+1], col = 'magenta')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[1]+1], col = 'black')      # traceplots of elephants in the same group but which are not mother-calf pairs
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[2]+1], col = 'tan')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[3]+1], col = 'orange')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[4]+1], col = 'green')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[6]+1], col = 'blue')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[7]+1], col = 'red')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[9]+1], col = 'purple')
lines(draws1_1.1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[10]+1], col = 'magenta')

# mean and credible interval for each dyad weight
means_1.1 <- data.frame(dyad = dyads$dyad,
                        mean = apply(draws_simdat_1.1[2:7141], 2, mean), 
                        median = apply(draws_simdat_1.1[2:7141], 2, median))
draws_1.1_pi <- apply(draws_simdat_1.1[2:7141], 2, rethinking::PI)
means_1.1$pi_lwr <- draws_1.1_pi[1,]
means_1.1$pi_upr <- draws_1.1_pi[2,]

# boxplot comparing types of dyads
plot_data_1.1 <- left_join(x = means_1.1, y = dyads, by = 'dyad')
ggplot(data = plot_data_1.1, aes(x = pair_type, y = mean, fill = pair_type))+
  geom_boxplot()+
  geom_jitter(aes(x = pair_type, y = edge, fill = pair_type),
              pch = 4, colour = col.alpha('red', 0.2))+
  scale_fill_viridis_d()+
  theme_classic()

# clear large objects from environment
rm(draws_1.1_pi, draws1_1.1, draws1_1.1_pi, draws2_1.1, draws3_1.1, draws4_1.1, means_1.1, plot_data_1.1, draws_simdat_1.1, dyads, simdat_ls)


################ 6) Run model on real standardised data ################
### run model
mod_1.1

### Fit model (slow) - just over an hour (1:06:30)
weight_motnp_1.1 <- mod_1.1$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_motnp_1.1
# variable       mean     median     sd    mad         q5        q95 rhat ess_bulk ess_tail
#lp__      -819804.08 -819806.00 248.32 249.82 -820201.05 -819385.90 1.00     1300     2048
#weight[1]       0.05       0.04   0.03   0.03       0.01       0.12 1.00     6812     2442
#weight[2]       0.05       0.04   0.03   0.03       0.01       0.12 1.00     6335     2736
#weight[3]       0.07       0.06   0.05   0.04       0.01       0.16 1.00     6288     2513
#weight[4]       0.05       0.04   0.03   0.03       0.01       0.12 1.00     7594     2252
#weight[5]       0.07       0.06   0.05   0.04       0.01       0.17 1.00     7396     2252
#weight[6]       0.05       0.04   0.03   0.03       0.01       0.12 1.00     7002     2412
#weight[7]       0.07       0.06   0.05   0.04       0.01       0.16 1.01     6784     2471
#weight[8]       0.05       0.04   0.03   0.03       0.01       0.11 1.00     7885     1953
#weight[9]       0.06       0.05   0.04   0.03       0.01       0.13 1.00     6779     2644
# showing 10 of 106954 rows (change via 'max_rows' argument or 'cmdstanr_max_rows' option)
#weight_motnp_1.1$summary()  # exceeds vector memory on my laptop
output1 <- read_cmdstan_csv(weight_motnp_1.1$output_files()[1])
output2 <- read_cmdstan_csv(weight_motnp_1.1$output_files()[2])
output3 <- read_cmdstan_csv(weight_motnp_1.1$output_files()[3])
output4 <- read_cmdstan_csv(weight_motnp_1.1$output_files()[4])

draws1_motnp1.1 <- as.data.frame(output1$post_warmup_draws)
draws2_motnp1.1 <- as.data.frame(output2$post_warmup_draws)
draws3_motnp1.1 <- as.data.frame(output3$post_warmup_draws)
draws4_motnp1.1 <- as.data.frame(output4$post_warmup_draws)
draws_motnp1.1 <- rbind(draws1_motnp1.1, draws2_motnp1.1, draws3_motnp1.1, draws4_motnp1.1)
rm(output1, output2, output3, output4)

colnames(draws_motnp1.1)[2:106954] <- counts_df$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:106954, size = 30, replace = F)

# tidy data -- vector memory usually exhausted and won't run
tidy_draws_1.1 <- pivot_longer(draws_motnp1.1[,2:106954], cols = everything(), names_to = 'dyad', values_to = 'draw')
tidy_draws_1.1$chain <- rep(1:4, each = 106953000)
tidy_draws_1.1$index <- rep(rep(1:1000, each = 106953),4)
head(tidy_draws_1.1, 10)
tail(tidy_draws_1.1, 10)

# subset for Sierra (F52) -- family = U17, F60, U21 and F98
tidy_sierra <- draws_motnp1.1[,c('F52_M40','F52_M26','F52_F8','F52_M15','F52_F98','F52_U17','F52_U21','F52_F60')]
tidy_sierra <- pivot_longer(tidy_sierra, cols = everything(), names_to = 'dyad', values_to = 'draw')
tidy_sierra <- separate(tidy_sierra, dyad, sep = '_', remove = F, into = c('F52','id'))
tidy_sierra$label <- ifelse(tidy_sierra$id == 'U17', 'U17 (calf)',
                         ifelse(tidy_sierra$id == 'U21', 'U21 (relative)',
                                ifelse(tidy_sierra$id == 'F60', 'F60 (relative)',
                                       ifelse(tidy_sierra$id == 'F98', 'F98 (relative)', 
                                              ifelse(tidy_sierra$id == 'F8', 'F8 (unrelated female)',
                                                     ifelse(tidy_sierra$id == 'M26', 'M26 (adult male)',
                                                            ifelse(tidy_sierra$id == 'M15', 'M15 (unrelated pubescent male)', 'M40 (adult male)')))))))
tidy_sierra$label <- factor(tidy_sierra$label, levels = c('U17 (calf)', 'F60 (relative)', 'U21 (relative)', 'F98 (relative)', 'M40 (adult male)','M26 (adult male)', 'F8 (unrelated female)','M15 (unrelated pubescent male)'))
tidy_sierra$chain <- rep(1:4, each = length(tidy_sierra$dyad)/4)
tidy_sierra$index <- rep(rep(1:1000, each = length(unique(tidy_sierra$dyad))),4)

### save data 
write_csv(draws_motnp1.1, 'data_processed/motnp_bayesian_edgedistributions_a1.b1_22.03.03.csv')
#write_csv(tidy_draws, 'data_processed/motnp_bayesian_edgedistributions_tidy_22.03.03.csv')

################ 7) Summarise and plot edge weights ################
# draws_motnp1.1 <- read_csv('data_processed/motnp_bayesian_edgedistributions_a1.b1_22.03.03.csv')
### build traceplots ####
# random dyads -- many are very wide, high uncertainty
plot(draws_motnp1.1[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_motnp1.1[,plot_cols[2]], col = 'tan')
lines(draws_motnp1.1[,plot_cols[3]], col = 'orange')
lines(draws_motnp1.1[,plot_cols[4]], col = 'green')
lines(draws_motnp1.1[,plot_cols[5]], col = 'chocolate')
lines(draws_motnp1.1[,plot_cols[6]], col = 'blue')
lines(draws_motnp1.1[,plot_cols[7]], col = 'red')
lines(draws_motnp1.1[,plot_cols[8]], col = 'seagreen')
lines(draws_motnp1.1[,plot_cols[9]], col = 'purple')
lines(draws_motnp1.1[,plot_cols[10]],col = 'magenta')
lines(draws_motnp1.1[,plot_cols[11]],col = 'black')
lines(draws_motnp1.1[,plot_cols[12]], col = 'tan')
lines(draws_motnp1.1[,plot_cols[13]], col = 'orange')
lines(draws_motnp1.1[,plot_cols[14]], col = 'green')
lines(draws_motnp1.1[,plot_cols[15]], col = 'chocolate')
lines(draws_motnp1.1[,plot_cols[16]], col = 'blue')
lines(draws_motnp1.1[,plot_cols[17]], col = 'red')
lines(draws_motnp1.1[,plot_cols[18]], col = 'seagreen')
lines(draws_motnp1.1[,plot_cols[19]], col = 'purple')
lines(draws_motnp1.1[,plot_cols[20]],col = 'magenta')
lines(draws_motnp1.1[,plot_cols[21]],col = 'black')
lines(draws_motnp1.1[,plot_cols[22]], col = 'tan')
lines(draws_motnp1.1[,plot_cols[23]], col = 'orange')
lines(draws_motnp1.1[,plot_cols[24]], col = 'green')
lines(draws_motnp1.1[,plot_cols[25]], col = 'chocolate')
lines(draws_motnp1.1[,plot_cols[26]], col = 'blue')
lines(draws_motnp1.1[,plot_cols[27]], col = 'red')
lines(draws_motnp1.1[,plot_cols[28]], col = 'seagreen')
lines(draws_motnp1.1[,plot_cols[29]], col = 'purple')
lines(draws_motnp1.1[,plot_cols[30]],col = 'magenta')

# Sierra = F52, herd members = F60+U21+F98, calf = U17. This looks really good, and not too wide (all seen a good number of times).
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,4000))
lines(draws_motnp1.1$F52_M40, col = 'black')      # non-herd member, adult male
lines(draws_motnp1.1$F52_M15, col = 'tan')        # non-herd member, pubescent male
lines(draws_motnp1.1$F52_M203,col = 'orange')     # non-herd member, adult male
lines(draws_motnp1.1$F52_M26, col = 'green')      # non-herd member, calf
lines(draws_motnp1.1$F52_F8,  col = 'chocolate')  # non-herd member, adult female
lines(draws_motnp1.1$F52_U9,  col = 'blue')       # non-herd member, calf
lines(draws_motnp1.1$F52_F98, col = 'red')        # herd member most frequently absent from sightings
lines(draws_motnp1.1$F52_U17, col = 'purple')     # calf
lines(draws_motnp1.1$F52_U21, col = 'seagreen')   # sister
lines(draws_motnp1.1$F52_F60, col = 'magenta')    # sister's calf

# ggplot version if tidying worked: plot for Sierra (F52)
(traceplot_sierra <- ggplot(tidy_sierra[tidy_sierra$chain == 1,], aes(x = index, y = draw, colour = label))+
    geom_line()+
    scale_color_viridis_d()+
    theme_classic()+
    scale_x_continuous('MCMC chain position',  expand = c(0,0))+
    scale_y_continuous('edge weight estimate', expand = c(0.02,0))+
    labs(colour = 'Partner (relationship)'))
ggsave(path = 'outputs/motnp_22.03.03/', filename = 'motnp_sierra_traceplot.pdf',
       plot = traceplot_sierra,
       width = 30, height = 24, units = 'cm')

(traceplot_sierra2 <- ggplot(tidy_sierra, aes(x = index, y = draw, colour = as.factor(chain)))+
    geom_line()+
    scale_color_viridis_d()+
    theme_light()+
    facet_wrap(.~label)+
    scale_x_continuous('MCMC chain position',  expand = c(0,0))+
    scale_y_continuous('edge weight estimate', expand = c(0.02,0))+
    theme(legend.position = c(0.85, 0.15), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))+
    labs(colour = 'Chain ID'))
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_sierra_traceplot.pdf',
       plot = traceplot_sierra,
       width = 30, height = 24, units = 'cm')

# check chain mixing
plot(draws1_motnp1.1$`1.weight[1]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_motnp1.1$`1.weight[1]`, col = 'green')
lines(draws3_motnp1.1$`1.weight[1]`, col = 'blue')
lines(draws4_motnp1.1$`1.weight[1]`, col = 'magenta')
plot(draws1_motnp1.1$`1.weight[2]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_motnp1.1$`1.weight[2]`, col = 'green')
lines(draws3_motnp1.1$`1.weight[2]`, col = 'blue')
lines(draws4_motnp1.1$`1.weight[2]`, col = 'magenta')
plot(draws1_motnp1.1$`1.weight[3]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_motnp1.1$`1.weight[3]`, col = 'green')
lines(draws3_motnp1.1$`1.weight[3]`, col = 'blue')
lines(draws4_motnp1.1$`1.weight[3]`, col = 'magenta')
which(counts_df$dyad == 'F52_U17') # 41192
plot(draws1_motnp1.1$`1.weight[41192]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_motnp1.1$`1.weight[41192]`, col = 'green')
lines(draws3_motnp1.1$`1.weight[41192]`, col = 'blue')
lines(draws4_motnp1.1$`1.weight[41192]`, col = 'magenta')
which(counts_df$dyad == 'F52_F60') # 40896
plot(draws1_motnp1.1$`1.weight[40896]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_motnp1.1$`1.weight[40896]`, col = 'green')
lines(draws3_motnp1.1$`1.weight[40896]`, col = 'blue')
lines(draws4_motnp1.1$`1.weight[40896]`, col = 'magenta')
which(counts_df$dyad == 'F52_U21') # 41197
plot(draws1_motnp1.1$`1.weight[41197]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_motnp1.1$`1.weight[41197]`, col = 'green')
lines(draws3_motnp1.1$`1.weight[41197]`, col = 'blue')
lines(draws4_motnp1.1$`1.weight[41197]`, col = 'magenta')
which(counts_df$dyad == 'F52_F98') # 40937
plot(draws1_motnp1.1$`1.weight[40937]`, col = 'purple', type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws2_motnp1.1$`1.weight[40937]`, col = 'green')
lines(draws3_motnp1.1$`1.weight[40937]`, col = 'blue')
lines(draws4_motnp1.1$`1.weight[40937]`, col = 'magenta')

### density plots ####
# random columns
dens(draws_motnp1.1[,plot_cols[1]], ylim = c(0,50), xlim = c(0,1), las = 1)
dens(add = T, draws_motnp1.1[,plot_cols[2]], col = 'tan')
dens(add = T, draws_motnp1.1[,plot_cols[3]], col = 'orange')
dens(add = T, draws_motnp1.1[,plot_cols[4]], col = 'green')
dens(add = T, draws_motnp1.1[,plot_cols[5]], col = 'chocolate')
dens(add = T, draws_motnp1.1[,plot_cols[6]], col = 'blue')
dens(add = T, draws_motnp1.1[,plot_cols[7]], col = 'red')
dens(add = T, draws_motnp1.1[,plot_cols[8]], col = 'seagreen')
dens(add = T, draws_motnp1.1[,plot_cols[9]], col = 'purple')
dens(add = T, draws_motnp1.1[,plot_cols[10]],col = 'magenta')
dens(add = T, draws_motnp1.1[,plot_cols[11]],col = 'black')
dens(add = T, draws_motnp1.1[,plot_cols[12]],col = 'tan')
dens(add = T, draws_motnp1.1[,plot_cols[13]],col = 'orange')
dens(add = T, draws_motnp1.1[,plot_cols[14]],col = 'green')
dens(add = T, draws_motnp1.1[,plot_cols[15]],col = 'chocolate')
dens(add = T, draws_motnp1.1[,plot_cols[16]],col = 'blue')
dens(add = T, draws_motnp1.1[,plot_cols[17]],col = 'red')
dens(add = T, draws_motnp1.1[,plot_cols[18]],col = 'seagreen')
dens(add = T, draws_motnp1.1[,plot_cols[19]],col = 'purple')
dens(add = T, draws_motnp1.1[,plot_cols[20]],col = 'magenta')
dens(add = T, draws_motnp1.1[,plot_cols[21]],col = 'black')
dens(add = T, draws_motnp1.1[,plot_cols[22]],col = 'tan')
dens(add = T, draws_motnp1.1[,plot_cols[23]],col = 'orange')
dens(add = T, draws_motnp1.1[,plot_cols[24]],col = 'green')
dens(add = T, draws_motnp1.1[,plot_cols[25]],col = 'chocolate')
dens(add = T, draws_motnp1.1[,plot_cols[26]],col = 'blue')
dens(add = T, draws_motnp1.1[,plot_cols[27]],col = 'red')
dens(add = T, draws_motnp1.1[,plot_cols[28]],col = 'seagreen')
dens(add = T, draws_motnp1.1[,plot_cols[29]],col = 'purple')
dens(add = T, draws_motnp1.1[,plot_cols[30]],col = 'magenta')
dens(add = T, draws_motnp1.1[,which(counts_df$dyad == 'F52_U17')+1],col = 'magenta')
dens(add = T, draws_motnp1.1[,which(counts_df$dyad == 'F52_F60')+1],col = 'purple')
dens(add = T, draws_motnp1.1[,which(counts_df$dyad == 'F52_U21')+1],col = 'blue')
dens(add = T, draws_motnp1.1[,which(counts_df$dyad == 'F52_F98')+1],col = 'red')

# density plots -- only run if tidying data worked
dens(tidy_sierra[tidy_sierra$id == 'M40',]$draw, xlab = 'edge weight estimation',       # all chains, plotted together
     ylim = c(0,60), xlim = c(0,1), las = 1)                                 # non-herd member, adult male
dens(tidy_sierra[tidy_sierra$id == 'M15',]$draw, add = T, col = 'tan')       # non-herd member, pubescent male
dens(tidy_sierra[tidy_sierra$id == 'M26',]$draw, add = T, col = 'green')     # non-herd member, calf
dens(tidy_sierra[tidy_sierra$id == 'F8',]$draw,  add = T, col = 'chocolate') # non-herd member, adult female
dens(tidy_sierra[tidy_sierra$id == 'F98',]$draw, add = T, col = 'red')       # herd member
dens(tidy_sierra[tidy_sierra$id == 'U17',]$draw, add = T, col = 'purple')    # calf
dens(tidy_sierra[tidy_sierra$id == 'U21',]$draw, add = T, col = 'seagreen')  # sister
dens(tidy_sierra[tidy_sierra$id == 'F60',]$draw, add = T, col = 'magenta')   # sister's calf
text('M40', x = 0.15, y = 05, col = 'black')
text('M15', x = 0.05, y = 55, col = 'tan')
text('M26', x = 0.05, y = 22, col = 'green')
text('F8' , x = 0.06, y = 18, col = 'chocolate')
text('F98', x = 0.60, y = 08, col = 'red')
text('U17', x = 0.95, y = 16, col = 'purple')
text('U21', x = 0.80, y = 10, col = 'seagreen')
text('F60', x = 0.86, y = 11, col = 'magenta')

(sierra_density <-
    ggplot(data = tidy_sierra[tidy_sierra$chain == 1], aes(x = draw, fill = ID))+       # only 1 chain
    geom_density()+  # why different values to the rethinking::dens() version?
    scale_fill_viridis_d(alpha = 0.4)+
    scale_x_continuous('edge weight estimation', expand = c(0,0), limits = c(0,1.1))+
    scale_y_continuous(expand = c(0,0), limits = c(0,55))+
    theme_classic()+
    theme(legend.position = c(0.5, 0.6)))
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_sierra_density.pdf',
       plot = sierra_density,
       width = 20, height = 24, units = 'cm')

(sierra_density <-
    ggplot(data = tidy_sierra, aes(x = draw, col = as.factor(chain)))+                 # all chains, plotted separately
    geom_density()+ 
    scale_fill_viridis_d(alpha = 0.2)+
    scale_x_continuous('edge weight estimation', expand = c(0,0), limits = c(0,1.1))+
    scale_y_continuous(expand = c(0,0), limits = c(0,55))+
    theme_light()+
    facet_wrap(. ~ ID)+
    theme(legend.position = c(0.85, 0.15), legend.title = element_text(size = 14))+
    labs(colour = 'Chain ID'))
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_sierra_density.pdf',
       plot = sierra_density,
       width = 20, height = 24, units = 'cm')

### summarise data ####
# summarise -- look for any anomalies in draw values or chain variation
summaries <- data.frame(dyad = colnames(draws_motnp1.1[2:106954]),
                        min = rep(NA, ncol(draws_motnp1.1)-1),
                        max = rep(NA, ncol(draws_motnp1.1)-1),
                        mean = rep(NA, ncol(draws_motnp1.1)-1),
                        median = rep(NA, ncol(draws_motnp1.1)-1),
                        sd = rep(NA, ncol(draws_motnp1.1)-1))
for(i in 72137:nrow(summaries)){
  summaries$min[i]    <- min(draws_motnp1.1[,i+1])
  summaries$max[i]    <- max(draws_motnp1.1[,i+1])
  summaries$mean[i]   <- mean(draws_motnp1.1[,i+1])
  summaries$median[i] <- median(draws_motnp1.1[,i+1])
  summaries$sd[i]     <- sd(draws_motnp1.1[,i+1])
}

summary(summaries)
summaries$dyad[which(summaries$min == min(summaries$min))] # F145_M68 = lowest value drawn
summaries$dyad[which(summaries$min == max(summaries$min))] # F17_U1 = highest lowest value drawn
summaries$dyad[which(summaries$max == min(summaries$max))] # M14_M56 = lowest highest value drawn
summaries$dyad[which(summaries$max == max(summaries$max))] # F44_U15 M131_M178 = highest value drawn (1)
summaries$dyad[which(summaries$sd == min(summaries$sd))]   # M14_M56
summaries$dyad[which(summaries$sd == max(summaries$sd))]   # F113_F117

plot(draws_motnp1.1$F145_M68, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')       # still shows a fair amount of variation
lines(draws_motnp1.1$F17_U1, col = 'red')                # still shows a fair amount of variation
lines(draws_motnp1.1$M14_M56, col = 'blue')              # very low variation, but enough to trust the chain
lines(draws_motnp1.1$F44_U15, col = 'green')             # still shows a fair amount of variation
lines(draws_motnp1.1$M131_M178, col = 'orange')          # looks pretty crazy! -- 2 young males, each only seen once
summaries$min[which(summaries$dyad == 'M131_M178')]
summaries$sd[which(summaries$dyad == 'M131_M178')]
summaries$mean[which(summaries$dyad == 'M131_M178')]
lines(draws_motnp1.1$F113_F117, col = 'purple')          # OK this one DOES look a bit crazy...

# organise dem_class
plot_data_motnp1.1 <- left_join(x = summaries, y = counts_df, by = 'dyad')
head(plot_data_motnp1.1)
plot_data_motnp1.1$dem_class_1_cat <- as.integer(as.factor(plot_data_motnp1.1$dem_class_1))
plot_data_motnp1.1$dem_class_2_cat <- as.integer(as.factor(plot_data_motnp1.1$dem_class_2))
plot_data_motnp1.1$dem_type_short <- ifelse(plot_data_motnp1.1$dem_class_1_cat <= plot_data_motnp1.1$dem_class_2_cat,
                                            paste(plot_data_motnp1.1$dem_class_1, plot_data_motnp1.1$dem_class_2, sep = '_'),
                                            paste(plot_data_motnp1.1$dem_class_2, plot_data_motnp1.1$dem_class_1, sep = '_'))
sort(unique(plot_data_motnp1.1$dem_type_short))
types <- data.frame(all = sort(unique(plot_data_motnp1.1$dem_type_short)), type1 = NA, type2 = NA)
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
plot_data_motnp1.1 <- left_join(plot_data_motnp1.1, types, by = 'dem_type_short')
plot_data_motnp1.1$join <- paste(plot_data_motnp1.1$join, ' ', sep = ' ')
which(is.na(plot_data_motnp1.1$join))

# boxplot edge weights by demographic type of dyad -- all types
colours <- c('magenta','purple','grey','grey','magenta','purple','grey','blue','purple','purple','purple','blue','purple','grey','grey','grey','grey','grey','magenta','purple','grey','purple','purple','blue','purple','grey','grey')
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_motnp1.1, aes(y = mean, x = join))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 fill = colours)+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown',
         y = 'mean edge weight')+
    coord_flip())
ggsave(path = 'outputs/motnp_22.03.03/', filename = 'motnp_all_edge_vs_demtype.pdf',
       plot = edge_vs_demtype_all,
       width = 20, height = 24, units = 'cm')  

# boxplot edge weights by demographic type of dyad -- only adults/pubescents
adults_motnp1.1 <- plot_data_motnp1.1[plot_data_motnp1.1$age_cat_id_1 > 2 & plot_data_motnp1.1$age_cat_id_2 > 2,]
adults_motnp1.1 <- adults_motnp1.1[!is.na(adults_motnp1.1$dyad),]
colours <- c('magenta','purple','magenta','purple','grey','blue','purple','blue','purple','magenta','purple','grey','blue','purple')
(edge_vs_demtype <- 
    ggplot(data = adults_motnp1.1, aes(y = mean, x = join))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 fill = colours)+
    theme_classic()+
    labs(x = 'dyad demography (A/P = adult/pubescent, M/F/U = male/female/unknown',
         y = 'mean edge weight')+
    coord_flip())
ggsave(path = 'outputs/motnp_22.03.03/', filename = 'motnp_adults_edge_vs_demtype.pdf',
       plot = edge_vs_demtype,
       width = 20, height = 24, units = 'cm')  

# scatter plot of all adult edges by age difference
adults_motnp1.1$sex_type_names <- ifelse(adults_motnp1.1$sex_type == 'F_F', 'female-female',
                                         ifelse(adults_motnp1.1$sex_type == 'F_M', 'male-female', 'male-male'))
(edge_vs_agediff <- 
    ggplot(adults_motnp1.1[adults_motnp1.1$sex_1 != 'U' & adults_motnp1.1$sex_2 != 'U',],
           aes(x = age_diff, y = mean, fill = sex_type))+
    geom_jitter(alpha = 0.2)+
    theme_classic()+
    theme(legend.position = 'none')+
    facet_grid(~sex_type_names)+
    scale_x_continuous('age category difference between dyad members',
                       expand = c(0.02,0),
                       breaks = c(0:5))+
    scale_y_continuous('mean edge weight',
                       expand = c(0.02,0),
                       limits = c(0,1)))
ggsave(path = 'outputs/motnp_22.03.03/', filename = 'motnp_adults_edge_vs_agediff.pdf',
       plot = edge_vs_agediff,
       width = 30, height = 24, units = 'cm')

rm(draws1_motnp1.1, draws2_motnp1.1, draws3_motnp1.1, draws4_motnp1.1)
rm(edge_vs_agediff, edge_vs_demtype, edge_vs_demtype_all, plot_data_motnp1.1, types, colours, adults_motnp1.1)

################ 8) Create network plots ################
head(summaries)
length(unique(plot_data_motnp1.1$id_1))+1 # number of individuals = 463

### Plot for Sierra and co. ####
counts_df_test <- counts_df[counts_df$id_1 == 'F52' | counts_df$id_1 == 'F60' | counts_df$id_1 == 'F98' | 
                              counts_df$id_1 == 'U17' | counts_df$id_1 == 'U21' | 
                              counts_df$id_1 == 'M40' | counts_df$id_1 == 'M26' | 
                              counts_df$id_1 == 'M15' | counts_df$id_1 == 'F8',]
counts_df_test <- counts_df_test[counts_df_test$id_2 == 'F52' | counts_df_test$id_2 == 'F60' | counts_df_test$id_2 == 'F98' | 
                                   counts_df_test$id_2 == 'U17' | counts_df_test$id_2 == 'U21' | 
                                   counts_df_test$id_2 == 'M40' | counts_df_test$id_2 == 'M26' | 
                                   counts_df_test$id_2 == 'M15' | counts_df_test$id_2 == 'F8',]
sierra_dyads <- counts_df_test$dyad
all_dyads <- colnames(draws_motnp1.1)
draws_test <- draws_motnp1.1[,which(all_dyads == "F52_F60" | all_dyads == "F52_F8" | all_dyads == "F52_F98" | all_dyads == "F52_M15" |
                                      all_dyads == "F52_M26" | all_dyads == "F52_M40" | all_dyads == "F52_U17" | all_dyads == "F52_U21" |
                                      all_dyads == "F60_F8"  | all_dyads == "F60_F98" | all_dyads == "F60_M15" | all_dyads == "F60_M26" |
                                      all_dyads == "F60_M40" | all_dyads == "F60_U17" | all_dyads == "F60_U21" | all_dyads == "F8_F98"  |
                                      all_dyads == "F8_M15"  | all_dyads == "F8_M26"  | all_dyads == "F8_M40"  | all_dyads == "F8_U17"  |
                                      all_dyads == "F8_U21"  | all_dyads == "F98_M15" | all_dyads == "F98_M26" | all_dyads == "F98_M40" |
                                      all_dyads == "F98_U17" | all_dyads == "F98_U21" | all_dyads == "M15_M26" | all_dyads == "M15_M40" |
                                      all_dyads == "M15_U17" | all_dyads == "M15_U21" | all_dyads == "M26_M40" | all_dyads == "M26_U17" |
                                      all_dyads == "M26_U21" | all_dyads == "M40_U17" | all_dyads == "M40_U21" | all_dyads == "U17_U21")]
head(counts_df_test)

# create matrix of edge weights
adj_tensor_test <- array(0, c(9, 9, 4000))
counts_df_test$dyad_id2 <- as.integer(as.factor(counts_df_test$dyad_id))    # create new variable so index variable has no breaks
counts_df_test$node_1_2 <- as.integer(as.factor(counts_df_test$node_1))     # create new variable so node1 variable has no breaks
counts_df_test$node_2_2 <- as.integer(as.factor(counts_df_test$node_2))+1   # create new variable so node2 variable has no breaks
for(dyad_id in 1:nrow(counts_df_test)) {                                    # assign matrix data as draws from chains (all 4 chains included)
  dyad_row <- counts_df_test[counts_df_test$dyad_id2 == dyad_id, ]
  adj_tensor_test[dyad_row$node_1_2, dyad_row$node_2_2, ] <- draws_test[, dyad_id]
}
adj_tensor_test[, , 1] # Print the first sample of the posterior distribution over adjacency matrices

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles_test <- apply(adj_tensor_test, c(1, 2), function(x) quantile(x, probs=c(0.025, 0.5, 0.975))) # identify median and 95% CI
(adj_lower_test <- adj_quantiles_test[1, , ])
(adj_mid_test   <- adj_quantiles_test[2, , ])
(adj_upper_test <- adj_quantiles_test[3, , ])

# Calculate standardised width/range of credible intervals.
adj_range_test <- ((adj_upper_test - adj_lower_test)/adj_mid_test)
adj_range_test[is.nan(adj_range_test)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid_test <- graph_from_adjacency_matrix(adj_mid_test,   mode="undirected", weighted=TRUE)
g_rng_test <- graph_from_adjacency_matrix(adj_range_test, mode="undirected", weighted=TRUE)

# Plot the median graph first and then the standardised width graph to show uncertainty over edges.
coords_test <- igraph::layout_nicely(g_mid_test)
plot(g_mid_test, edge.width = 3*E(g_mid_test)$weight, edge.color = "black",  layout = coords_test)
plot(g_mid_test, edge.width = 3*E(g_rng_test)$weight, edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = c(sort(unique(counts_df_test$id_1)),'U21'),
     vertex.label.dist = 0, vertex.label.color = "black", layout = coords_test, add = TRUE)

# clear environment
rm(adj_lower_test, adj_mid_test, adj_range_test, adj_upper_test, coords_test, counts_df_test, draws_test, dyad_row, g_mid_test, g_rng_test, adj_quantiles_test, adj_tensor_test, dyad_id, i, plot_cols, all_dyads, sierra_dyads)

### Create igraph object for all elephants and then just plot certain individuals ####
# read in draws data -- very slow!
#draws_motnp1.1 <- read_csv('data_processed/motnp_bayesian_edgedistributions_a1.b1_22.03.03.csv')
#draws <- data.matrix(draws_motnp1.1)

# reassign dyad numbers to remove gaps
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1

# create array of draws per dyad (distributions)
adj_tensor <- array(0, c(NROW(unique(counts_df$id_1))+1,
                         NROW(unique(counts_df$id_2))+1,
                         NROW(draws)),
                    dimnames = list(c(unique(counts_df$id_1),'U9'),
                                    c('F1',unique(counts_df$id_2)),
                                    NULL))
N <- nrow(counts_df)

for (i in 1:54000) {            # can cope with jumps of 50000 dyads at a time
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws[, i+1]
}
for (i in 54001:N) {
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws[, i+1]
}
adj_tensor[,,1]

# convert to array of medians and 95% credible intervals
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- ((adj_upper - adj_lower)/adj_mid) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Generate nodes data for plotting characteristics
ele_nodes <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
nodes <- data.frame(id = sort(unique(ele_nodes$id)))
nodes <- left_join(nodes, ele_nodes, by = 'id')
nodes$sex       <- as.factor(nodes$sex)
nodes$age_class <- as.factor(nodes$age_class)
nodes$dem_class <- as.factor(nodes$dem_class)
str(nodes)

# correct age and demographic classes for nodes
unique(nodes$age_category) # "50+"   "10-15" "35-50" "20-35" "15-19" "8-9"   "9-10"  "4-5"   "5-6"   "6-7"   "0-3"   "7-8"   "20-25" "25-40" "40+"   "3-4" "1-2"   NA
nodes$age_category <- ifelse(nodes$age_category == '1-2','0-3',nodes$age_category)
nodes$age_cat_id <- ifelse(nodes$age_category == '0-3', 1,
                           ifelse(nodes$age_category == '3-4', 1,
                                  ifelse(nodes$age_category == '4-5', 1,
                                         ifelse(nodes$age_category == '5-6', 2,
                                                ifelse(nodes$age_category == '6-7', 2,
                                                       ifelse(nodes$age_category == '7-8', 2,
                                                              ifelse(nodes$age_category == '8-9', 2,
                                                                     ifelse(nodes$age_category == '9-10', 2,
                                                                            ifelse(nodes$age_category == '10-15', 3,
                                                                                   ifelse(nodes$age_category == '15-19', 4,
                                                                                          ifelse(nodes$age_category == '20-25', 5,
                                                                                                 ifelse(nodes$age_category == '20-35', 5,
                                                                                                        ifelse(nodes$age_category == '25-40', 6,
                                                                                                               ifelse(nodes$age_category == '35-50', 6,
                                                                                                                      ifelse(nodes$age_category == '40+', 7,
                                                                                                                             ifelse(nodes$age_category == '50+', 7, nodes$age_category))))))))))))))))
nodes[is.na(nodes$age_category),]   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1
nodes$age_cat_id[which(is.na(nodes$age_cat_id))] <- 1

unique(nodes$age_category[nodes$age_class == 'Calf'])      # shouldn't include any ages over 4-5
unique(nodes$age_category[nodes$age_class == 'Juvenile'])  # shouldn't include any ages under 5-6
unique(nodes$age_category[nodes$age_class == 'Pubescent']) # shouldn't include any ages under 9-10 or over 15-19
unique(nodes$age_category[nodes$age_class == 'Adult'])     # shouldn't include any ages under 20-25

nodes$age_class <- ifelse(nodes$age_cat_id == 1, 'Calf',
                          ifelse(nodes$age_cat_id == 2, 'Juvenile',
                                 ifelse(nodes$age_cat_id > 4, 'Adult','Pubescent')))

nodes$dem_class <- ifelse(nodes$age_class == 'Adult', paste('A',nodes$sex, sep = ''),
                          ifelse(nodes$age_class == 'Pubescent', paste('P',nodes$sex, sep = ''),
                                 ifelse(nodes$age_class == 'Juvenile', paste('J',nodes$sex, sep = ''),
                                        paste('C',nodes$sex, sep = ''))))

rm(ele_nodes)

# Plot whole network
coords <- igraph::layout_nicely(g_mid)
plot(g_mid,
     edge.width = E(g_mid)$weight,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = 'black',
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight,
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.size = 8,
     vertex.label = c(sort(unique(counts_df$id_1)),'U9'),
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(nodes$age_class == 'Adult','seagreen1',
                          ifelse(nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)


plot(g_mid,
     edge.width = E(g_mid)$weight*3,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(adj_mid < 0.5,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight*3,
     edge.color = ifelse(adj_mid < 0.5,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 5,
     vertex.label = c(sort(unique(counts_df$id_1)),'U9'),
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(nodes$age_class == 'Adult','seagreen1',
                          ifelse(nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

# create variables for different degrees of node connectedness
nodes$degree_0.1 <- NA
nodes$degree_0.2 <- NA
nodes$degree_0.3 <- NA
nodes$degree_0.4 <- NA
nodes$degree_0.5 <- NA

summaries <- separate(summaries, dyad, c('id_1','id_2'), '_', F)
for(i in 1:NROW(nodes)){
  rows <- summaries[summaries$id_1 == nodes$id[i] | summaries$id_2 == nodes$id[i],]
  nodes$degree_0.1[i] <- length(which(rows$median > 0.1))
  nodes$degree_0.2[i] <- length(which(rows$median > 0.2))
  nodes$degree_0.3[i] <- length(which(rows$median > 0.3))
  nodes$degree_0.4[i] <- length(which(rows$median > 0.4))
  nodes$degree_0.5[i] <- length(which(rows$median > 0.5))
}

which(nodes$degree_0.1 < nodes$degree_0.2)
which(nodes$degree_0.2 < nodes$degree_0.3)
which(nodes$degree_0.3 < nodes$degree_0.4)
which(nodes$degree_0.4 < nodes$degree_0.5)

# plot network with reduced nodes -- only those with degree values > 0.5
g_mid_0.5 <- delete.vertices(graph = g_mid, v = nodes$id[which(nodes$degree_0.5 == 0)])
g_rng_0.5 <- delete.vertices(graph = g_rng, v = nodes$id[which(nodes$degree_0.5 == 0)])

coords_0.5 <- layout_nicely(g_mid_0.5)
plot(g_mid_0.5,
     edge.width = ifelse(E(g_mid_0.5)$weight > 0.5, E(g_rng_0.5)$weight, 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 6,
     vertex.label = NA,
     layout = coords_0.5)
plot(g_mid_0.5,
     edge.width = ifelse(E(g_mid_0.5)$weight > 0.5, E(g_mid_0.5)$weight, 0),
     edge.color = 'black',
     vertex.size = 6,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(nodes[which(nodes$degree_0.5 != 0),]$age_class == 'Adult', 'seagreen1',
                           ifelse(nodes[which(nodes$degree_0.5 != 0),]$age_class == 'Pubescent','skyblue',
                                  ifelse(nodes[which(nodes$degree_0.5 != 0),]$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords_0.5,
     add = TRUE)

#  plot network with reduced nodes -- only those with degree values > 0.3
g_mid_0.3 <- delete.vertices(graph = g_mid, v = nodes$id[which(nodes$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng, v = nodes$id[which(nodes$degree_0.3 == 0)])

coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight*3, 0),
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(nodes[which(nodes$degree_0.3 != 0),]$age_class == 'Adult',
                           'seagreen1',
                           ifelse(nodes[which(nodes$degree_0.3 != 0),]$age_class == 'Pubescent',
                                  'skyblue',
                                  ifelse(nodes[which(nodes$degree_0.3 != 0),]$age_class == 'Juvenile',
                                         'yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

#plot males
g_mid_m <- delete.vertices(graph = g_mid,
                           v = nodes$id[which(nodes$dem_class != 'AM' &
                                                    nodes$dem_class != 'PM')])
g_rng_m <- delete.vertices(graph = g_rng,
                           v = nodes$id[which(nodes$dem_class != 'AM' &
                                                    nodes$dem_class != 'PM')])
males <- nodes[nodes$dem_class == 'AM' | nodes$dem_class == 'PM',]
coords_m <- layout_nicely(g_mid_m)
plot(g_mid_m,
     edge.width = E(g_mid_m)$weight,
     edge.color = 'black',
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m)
plot(g_mid_m,
     edge.color = rgb(0,0,0,0.25),
     edge.width = E(g_rng_m)$weight,
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males$age_class == 'Adult','seagreen1', 'skyblue'),
     layout = coords_m, add = TRUE)

# plot males >0.3
g_mid_m0.3 <- delete.vertices(graph = g_mid,
                           v = nodes$id[which(nodes$dem_class != 'AM' &
                                                    nodes$dem_class != 'PM' |
                                                    nodes$degree_0.3 == 0)])
g_rng_m0.3 <- delete.vertices(graph = g_rng,
                           v = nodes$id[which(nodes$dem_class != 'AM' &
                                                    nodes$dem_class != 'PM' |
                                                    nodes$degree_0.3 == 0)])
males0.3 <- males[males$degree_0.3 > 0,]
coords_m0.3 <- layout_nicely(g_mid_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_rng_m0.3)$weight*3, 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_mid_m0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males0.3$age_class == 'Adult', 'seagreen1', 'skyblue'),
     layout = coords_m0.3, add = TRUE)

m38 <- summaries[summaries$id_1 == 'M38' | summaries$id_2 == 'M38',]
m38 <- m38[m38$median > 0.3,] # leaves on elephants which are >0.3 connected to females/calves
m87 <- summaries[summaries$id_1 == 'M87' | summaries$id_2 == 'M87',]
m87 <- m87[m87$median > 0.3,] # leaves on elephants which are >0.3 connected to females/calves

g_mid_m0.3 <- delete_vertices(graph = g_mid_m0.3, v = c('M240','M71','M99','M233',
                                                        'M158','M148','M87','M38',
                                                        'M128','M10','M111','M95',
                                                        'M7','M70','M17','M117',
                                                        'M221','M156','M196','M80',
                                                        'M114','M177'))
g_rng_m0.3 <- delete_vertices(graph = g_rng_m0.3, v = c('M240','M71','M99','M233',
                                                        'M158','M148','M87','M38',
                                                        'M128','M10','M111','M95',
                                                        'M7','M70','M17','M117',
                                                        'M221','M156','M196','M80',
                                                        'M114','M177'))
males0.3 <- males0.3[males0.3$id != 'M240' & males0.3$id != 'M71' & males0.3$id != 'M99' & males0.3$id != 'M233' & males0.3$id != 'M158' & males0.3$id != 'M148' & males0.3$id != 'M87' & males0.3$id != 'M38' & males0.3$id != 'M128' & males0.3$id != 'M10' & males0.3$id != 'M111' & males0.3$id != 'M95' & males0.3$id != 'M7' & males0.3$id != 'M70' & males0.3$id != 'M17' & males0.3$id != 'M117' & males0.3$id != 'M221' & males0.3$id != 'M156' & males0.3$id != 'M196' & males0.3$id != 'M80' & males0.3$id != 'M114' & males0.3$id != 'M177',]

coords_m0.3 <- layout_nicely(g_mid_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_rng_m0.3)$weight*3, 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_mid_m0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males0.3$age_class == 'Adult', 'seagreen1', 'skyblue'),
     layout = coords_m0.3, add = TRUE)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_rng_m0.3)$weight*3, 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_mid_m0.3)$weight*3, 0),
     edge.color = 'black',
     vertex.size = males0.3$count,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males0.3$age_class == 'Adult', 'seagreen1', 'skyblue'),
     layout = coords_m0.3, add = TRUE)

summary(males0.3$count)
table(males0.3$count)
males0.3$count[which(males0.3$id == 'M157')]
males0.3$count[which(males0.3$id == 'M162')]
summaries$sd[which(summaries$dyad == 'M157_M162')] # 0.2355704
summaries$sd[which(summaries$dyad == 'M101_M16')]  # 0.2371941
summaries$sd[which(summaries$dyad == 'M200_M205')] # 0.1335143
max(summaries$sd)
plot(draws[,which(colnames(draws) == 'M157_M162')], type = 'l')
plot(draws[,which(colnames(draws) == 'M101_M16')],  type = 'l')
plot(draws[,which(colnames(draws) == 'M200_M205')], type = 'l')
g_rng_m0.3[1] # M157 = 20, M162 = 24
g_rng_m0.3[20] # M162 = 1.147962
g_rng_m0.3[24] # M157 = 1.147962
g_rng_m0.3[1]  ; g_rng_m0.3[22] # M101_M16  = 1.166944
g_rng_m0.3[47] ; g_rng_m0.3[48] # M200_M205 = 1.579656
g_mid_m0.3[20] ; g_mid_m0.3[24] # M157_M162 = 0.71228400
g_mid_m0.3[1]  ; g_mid_m0.3[22] # M101_M16  = 0.71306150
g_mid_m0.3[47] ; g_mid_m0.3[48] # M200_M205 = 0.32397950
# is plotting igraph object correctly, igraph is underestimating the uncertainty

adj_range['M157','M162']  # 1.147962
adj_range['M101','M16']   # 1.166944
adj_range['M200','M205']  # 1.579656
# igraph is not making the error, adj_range is
adj_quantiles[1,'M157','M162']  # 0.1706227
adj_quantiles[2,'M157','M162']  # 0.712284
adj_quantiles[3,'M157','M162']  # 0.9882974
adj_quantiles[1,'M101','M16']   # 0.1569084
adj_quantiles[2,'M101','M16']   # 0.7130615
adj_quantiles[3,'M101','M16']   # 0.989011
adj_quantiles[1,'M200','M205']  # 0.1040737
adj_quantiles[2,'M200','M205']  # 0.3239795
adj_quantiles[3,'M200','M205']  # 0.6158498


adj_range2 <- adj_upper - adj_lower ; adj_range2[is.nan(adj_range2)] <- 0
g_rng2 <- graph_from_adjacency_matrix(adj_range2, mode="undirected", weighted=TRUE)
g_rng2_m0.3 <- delete.vertices(graph = g_rng2,
                               v = nodes$id[which(nodes$dem_class != 'AM' &
                                                   nodes$dem_class != 'PM' |
                                                   nodes$degree_0.3 == 0)])

g_rng2_m0.3 <- delete_vertices(graph = g_rng2_m0.3, v = c('M240','M71','M99','M233',
                                                        'M158','M148','M87','M38',
                                                        'M128','M10','M111','M95',
                                                        'M7','M70','M17','M117',
                                                        'M221','M156','M196','M80',
                                                        'M114','M177'))
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, (E(g_rng2_m0.3)$weight*5), 0),
     edge.color = rgb(0,0,0,0.25),
     vertex.size = 1,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_mid_m0.3)$weight*5, 0),
     edge.color = 'black',
     vertex.size = males0.3$count,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males0.3$age_class == 'Adult', 'seagreen1', 'skyblue'),
     layout = coords_m0.3, add = TRUE)
g_mid_m0.3[1]
g_rng_m0.3[1]

plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, (E(g_rng2_m0.3)$weight*10), 0),
     edge.color = 'blue',
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_m0.3)
plot(g_mid_m0.3,
     edge.width = ifelse(E(g_mid_m0.3)$weight > 0.3, E(g_mid_m0.3)$weight*10, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(males0.3$age_class == 'Adult', 'seagreen1', 'skyblue'),
     layout = coords_m0.3, add = TRUE)


### clear environment
rm(coords, coords_0.3, coords_0.5, coords_m, coords_m0.3, dyad_row, edges, ele, ele_nodes, elephants, g_mid, g_mid_0.3, g_mid_0.5, g_mid_m, g_mid_m0.3, g_rng, g_rng_0.3, g_rng_0.5, g_rng_m, g_rng_m0.3, m38, m87, males, males0.3, rows, traceplot_sierra, traceplot_sierra2, i)

################ 9) Just playing with gregariousness ideas a little ################
head(summaries)
length(summaries$dyad) # this is all of the dyads
edges <- left_join(x = summaries, y = counts_df, by = 'dyad')
edges <- edges[,c(1:9,12:41)]
colnames(edges)[2:3] <- c('id_1','id_2')
elephants <- data.frame(id = c(edges$id_1, 'U9'),
                        sex = c(edges$sex_1, 'U'),
                        age_class = c(edges$age_class_1, 'Calf'),
                        age_cat_id = c(edges$age_cat_id_1, 1),
                        mean_edge = rep(NA, nrow(edges)+1),
                        med_edge  = rep(NA, nrow(edges)+1),
                        mean_edge_apm = rep(NA, nrow(edges)+1),
                        mean_edge_apf = rep(NA, nrow(edges)+1),
                        max_edge = rep(NA, nrow(edges)+1),
                        max_edge_dem = rep(NA, nrow(edges)+1))
elephants <- distinct(elephants)
for(i in 1:nrow(elephants)){
  ele <- edges[edges$id_1 == elephants$id[i] | edges$id_2 == elephants$id[i],]
  elephants$mean_edge[i] <- mean(ele$mean)
  elephants$med_edge[i] <- median(ele$median)
  elephants$mean_edge_apm[i] <- mean(ele$mean[which(ele$dem_class_1 == 'AM' | ele$dem_class_1 == 'PM' | 
                                                      ele$dem_class_2 == 'AM' | ele$dem_class_2 == 'PM')])
  elephants$mean_edge_apf[i] <- mean(ele$mean[which(ele$dem_class_1 == 'AF' | ele$dem_class_1 == 'PF' |
                                                      ele$dem_class_2 == 'AF' | ele$dem_class_2 == 'PF')])
  elephants$max_edge[i] <- max(ele$mean)
  elephants$max_edge_dem[i] <- ifelse(ele$id_1[which(ele$mean == elephants$max_edge[i])] == elephants$id[i],
                                      ele$dem_class_2[which(ele$mean == elephants$max_edge[i])],
                                      ele$dem_class_1[which(ele$mean == elephants$max_edge[i])])
}

summary(elephants$max_edge)

ggplot(elephants, aes(x = mean_edge, y = med_edge, colour = sex))+
  geom_point()+
  xlab('mean edge weight')+
  ylab('median edge weight')+
  theme_light()+
  theme(legend.position = c(0.8,0.3))

ggplot(elephants, aes(x = sex, y = mean_edge, fill = sex))+
  geom_boxplot(notch = T)+
  theme_light()+
  theme(legend.position = 'none')+
  ylab('mean edge weight')

ggplot(elephants, aes(x = age_cat_id, y = mean_edge, colour = sex))+
  #geom_point(alpha = 0.7)+
  geom_jitter(width = 0.2)+
  xlab('age category')+
  ylab('median edge weight')+
  theme_light()

ggplot(elephants, aes(x = age_cat_id, y = mean_edge_apm, colour = sex))+
  #geom_point(alpha = 0.7)+
  geom_jitter(width = 0.2)+
  xlab('age category')+
  ylab('median edge weight')+
  theme_light()

ggplot(elephants, aes(x = age_cat_id, y = mean_edge_apm, colour = sex))+
  #geom_point(alpha = 0.7)+
  geom_jitter(width = 0.2)+
  xlab('age category')+
  ylab('median edge weight')+
  theme_light()

ggplot(elephants, aes(x = age_cat_id, y = max_edge, pch = sex, colour = max_edge_dem))+
  #geom_point(alpha = 0.7)+
  geom_jitter(width = 0.2)+
  labs(colour = 'demography of\nclosest associate',
       x = 'age category',
       y = 'maximum edge weight')+
  theme_light()


################ 10) Nodal regression ################
# The final common type of network analysis we will cover here is nodal regression, where a regression is performed to analyse the relationship between a nodal network metric(such as centrality) and nodal traits (such as age and sex).  These analyses are usually used to assess how network position depends on various biological factors, but can also be used where network position is a predictor.  Since node metrics are derivative measures of the network, uncertainty in edge weights should ideally propagate through to node metrics, and on to coefficients in regression analyses, giving an accurate estimate of the total uncertainty in inferred parameters.  The core of the nodal regression is similar to dyadic regression, taking the form:
#X(n)i j∼Poisson(λ(n)i jD(n)i j)
#logÄλ(n)i jä=ωi j+L(n)i j
#Mi(ω)∼Normal(β0+β1xi,σ2)
# where β0 is the intercept parameter, β1 is the slope parameter, and σ is the standard deviation of the residuals. Mi(ω) denotes the node metric estimates when applied to the edge weights estimated in the top part of the model. Calculating node metrics within the model may present a practical challenge when using standard model fitting software, as many node metrics can not be described purely in terms of closed-form equations (Freeman, 1978).  In this case, splitting up the model along the dashed line becomes an important step, as it allows the edge weights to be estimated using a standard piece of software, and leaves the regression part of the model to be fit using a sampler that supports numeric estimates of network centralities. As part of the code we have made available, we have written a custom Metropolis-Hastings sampler that allows numeric estimates of node metrics to be used when fitting the model. Importantly, our sampler samples from the joint distribution of edge weights, rather than sampling edge weights independently. This ensures that the effect of any structure in the observation data (such as location effects) is maintained and propagated through to the node metric estimates, and subsequently to regression coefficient estimates.

### create array of draws per dyad (distributions) -- same as above
adj_tensor <- array(0, c(NROW(unique(counts_df$id_1))+1,
                         NROW(unique(counts_df$id_2))+1,
                         NROW(draws_motnp1.1)),
                    dimnames = list(c(unique(counts_df$id_1),'U9'),
                                    c('F1',unique(counts_df$id_2)),
                                    NULL))
N <- nrow(counts_df)

for (i in 1:54000) {            # can cope with jumps of 50000 dyads at a time
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_motnp1.1[, i+1]
}
for (i in 54001:N) {
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_motnp1.1[, i+1]
}
adj_tensor[,,1]

### create centrality matrix
centrality_matrix <- matrix(0, nrow = NROW(draws_motnp1.1), ncol = length(unique(counts_df$id_1))+1)
for (i in 1:NROW(draws_motnp1.1)) {
  g <- graph_from_adjacency_matrix(adj_tensor[, , i], mode="undirected", weighted=TRUE)
  centrality_matrix[i, ] <- strength(g)
}

colnames(centrality_matrix) <- c(sort(unique(counts_df$id_1)),unique(counts_df$id_2)[462])
head(centrality_matrix)
write_csv(as.data.frame(centrality_matrix), 'data_processed/motnp_bayesian_nodalregression_centralitymatrix_22.03.04.csv')

centrality_quantiles <- t(apply(centrality_matrix, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
centrality_quantiles

cq <- data.frame(id = rownames(centrality_quantiles),
                 strength_lwr = centrality_quantiles[,1],
                 strength_mid = centrality_quantiles[,2],
                 strength_upr = centrality_quantiles[,3])
cq$id[which(cq$strength_mid == max(cq$strength_mid))]
nodes <- left_join(nodes, cq, by = 'id')
head(nodes)

boxplot(nodes$strength_mid ~ nodes$sex, notch = T)
boxplot(nodes$strength_mid ~ nodes$age_class, notch = T)
boxplot(nodes$strength_mid ~ nodes$age_category, notch = T, horizontal = T, las = 1, ylab = '')

nodes$age <- ifelse(nodes$age_category == '0-3', 1,
             ifelse(nodes$age_category == '1-2', 1,
             ifelse(nodes$age_category == '3-4', 1,
             ifelse(nodes$age_category == '4-5', 1,
             ifelse(nodes$age_category == '5-6', 2,
             ifelse(nodes$age_category == '6-7', 2,
             ifelse(nodes$age_category == '7-8', 2,
             ifelse(nodes$age_category == '8-9', 2,
             ifelse(nodes$age_category == '9-10', 3,
             ifelse(nodes$age_category == '10-15', 3,
             ifelse(nodes$age_category == '15-19', 4,
             ifelse(nodes$age_category == '20-25', 5,
             ifelse(nodes$age_category == '20-35', 5,
             ifelse(nodes$age_category == '25-40', 6,
             ifelse(nodes$age_category == '35-50', 6,
             ifelse(nodes$age_category == '40+', 7,
             ifelse(nodes$age_category == '50+', 7, 'error')))))))))))))))))
nodes$age <- ifelse(nodes$id == 'U8', 1, nodes$age)
boxplot(nodes$strength_mid ~ nodes$age, notch = T, horizontal = T, las = 1, ylab = '')

nodes_tidy <- pivot_longer(data = nodes, cols = 14:21,
                           names_to = 'centrality_type', values_to = 'centrality')
nodes_tidy$centrality_type <- factor(nodes_tidy$centrality_type,
                                     levels = c('strength_lwr', 'strength_mid', 'strength_upr',
                                                'degree_0.1','degree_0.2','degree_0.3',
                                                'degree_0.4','degree_0.5'))
ggplot(nodes_tidy, aes(y = centrality, x = age, fill = centrality_type))+
  geom_boxplot()+
  facet_wrap(.~centrality_type, scales = 'free')+ # lwr, mid and upr = very similar patterns. degree = very variable depending on level
  theme_light()+
  theme(legend.position = 'none',
        strip.text = element_text(colour = 'black'))

ggplot(nodes_tidy[nodes_tidy$centrality_type == 'degree_0.2' |
                    nodes_tidy$centrality_type == 'strength_mid',],
       aes(y = centrality, x = age))+
  geom_boxplot()+
  facet_wrap( . ~ centrality_type )

### write nodes data frame to file
write_csv(nodes, 'data_processed/motnp_elenodes_centrality_22.03.03.csv')

# The key challenge to quantifying uncertainty in network analysis is to incorporate uncertainty due to sampling into downstream analyses, commonly regression. This can be achieved by modifying the likelihood function of a regression model to treat the network centralities with uncertainty. We have written a custom MCMC sampler function that samples from the joint distribution of network centralities calculated earlier and treats those samples as the data in the likelihood function. Likelihood functions for the sampler use the index variable to keep track of which data points are being compared internally in the sampler, to ensure that candidate steps in the MCMC are not accepted or rejected because they are being compared to different data points, rather than because of the parameter space. Custom likelihood functions take a similar form to the target += syntax in Stan, but for more specific resources the following document is a good start: https://www.ime.unicamp.br/~cnaber/optim_1.pdf.

loglik <- function(params, Y, X, index) {
  # Define parameters
  intercept <- params[1]
  beta_age <- params[2]
  beta_females <- params[3]
  beta_males <- params[4]
  beta_unks <- params[5]
  sigma <- exp(params[6]) # Exponential keeps underlying value unconstrained, which is much easier for the sampler.
  
  # Sample data according to index
  y <- Y[index %% dim(Y)[1] + 1, ]
  
  # Define model
  target <- 0
  target <- target + sum(dnorm(y, mean = intercept + beta_age*X[,1] + beta_females*X[,2] + beta_males*X[,3] + beta_unks*X[,4],
                               sd = sigma, log = TRUE))                   # Main model
  target <- target + dnorm(intercept,    mean = 0, sd = 2.5, log = TRUE)  # Prior on intercept
  target <- target + dnorm(beta_age,     mean = 0, sd = 2.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_females, mean = 0, sd = 2.5, log = TRUE)  # Prior on female coefficient
  target <- target + dnorm(beta_males,   mean = 0, sd = 2.5, log = TRUE)  # Prior on male coefficient
  target <- target + dnorm(beta_unks,    mean = 0, sd = 2.5, log = TRUE)  # Prior on unknown coefficient
  target <- target + dexp(sigma, 1, log = TRUE)                           # Prior on sigma
  
  return(target)
}

# Now we will prepare data for fitting the model. The predictor matrix is simply a matrix with 3 columns and 463 rows, corresponding to whether each of the 8 nodes is a female (column 1), a male (column 2), or an unknown (column 3).
predictor_matrix <- matrix(0, nrow = length(nodes$id), ncol = 4)
colnames(predictor_matrix) <- c('age','female', 'male','unknown')
predictor_matrix[, 1] <- as.numeric(nodes$age)
predictor_matrix[which(nodes$sex == 'F'), 2] <- 1
predictor_matrix[which(nodes$sex == 'M'), 3] <- 1
predictor_matrix[which(nodes$sex == 'U'), 4] <- 1
predictor_matrix

write_csv(as.data.frame(predictor_matrix), 'data_processed/motnp_bayesian_nodalregression_predictormatrix_22.03.04.csv')

# Since network strength is strictly positive, a Gaussian error is not a reasonable model for the data. The Gaussian family model is much easier to implement as well as interpret than many other models, so we will standardise the centralities by taking z-scores.
centrality_matrix_std <- (centrality_matrix - apply(centrality_matrix, 1, mean))/apply(centrality_matrix, 1, sd)

# Now we’re in a position to fit the model. To do this, we define the target function, which is simply a function that maps candidate parameters and a network centrality index to the log-likelihood of that function for the given sample of the centrality posterior. This means the target function can be written as a function of the data centrality_matrix_std and predictor_matrix.
target <- function(params, index) loglik(params, centrality_matrix_std, predictor_matrix, index)

# define metropolis sampler function
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

# The function metropolis from sampler.R can now be used to fit the model using the provided target function, an initial set of parameters, and some additional MCMC options.
chain <- metropolis(target, c(0, 0, 0, 0), iterations=100000, thin=100, refresh=10000)

colnames(chain) <- c("intercept", "beta_lifeform", "beta_droid", "sigma")
head(chain)

par(mfrow=c(2, 2))
for (i in 1:4) {
  plot(chain[, i], type="l")
}

plot(density(centrality_matrix_std[1, ]), ylim=c(0, 0.7), main="", xlab="Standardised network strength")
sample_ids <- sample(1:1000, size=200)
for (i in sample_ids) {
  pred <- rnorm(8, mean=chain[i, "intercept"] + chain[i, "beta_lifeform"] * predictor_matrix[, 1] + chain[i, "beta_droid"] * predictor_matrix[, 2], sd=exp(chain[i, "sigma"]))
  lines(density(centrality_matrix_std[i, ]), col=rgb(0, 0, 0, 0.25))
  lines(density(pred), col=rgb(0, 0, 1, 0.25))
}

coefficient_quantiles <- t(apply(chain, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
coefficient_quantiles

beta_difference <- chain[, "beta_lifeform"] - chain[, "beta_droid"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))






#
################ 11) Dyadic regression ################
# A common type of analysis in studies of social networks is to consider factors that affectedge weights.  This type of analysis is often called dyadic regression, where explanatoryfactors such as sex or age difference are regressed against edge weight. The edge weightsin dyadic regression are non-independent, as variables such as age difference or sex dif-ference are inherently linked to individual-level attributes.  This means that effects dueto age or sex in individualiaffect all dyads that connecti.  This non-independence canbe controlled for by including node-level effects in the regression (Tranmer et al., 2014).Using the count data example discussed earlier, we propose to conduct dyadic regres-sion using a standard regression model of the form:
#         X_ij ∼ Poisson(lambda_ij)
#         log(lambda_ij) = weight_ij + Location_ij
#         weight_ij ∼ Normal(beta_0 + beta_1*x_ij + u_i + u_j, sigma_sq)
# where beta_0 is the intercept parameter, beta_1 is the slope parameter, ui and uj are parameters accounting for the effect of node membership on edge weight, and sigma is the standard deviation of the residuals.  The dashed line indicates that the model can be fit in two separate parts:the first part being the core edge weight model, and the second part being the dyadic regression model.  It is not strictly necessary to separate out the model like this, but do-ing so makes it possible to fit the edge weight model once and conduct multiple types of analysis on it afterwards without fitting the entire model again.Theβ0andβ1parameters can be interpreted in the same way as usual regression coefficients, and will be accompanied with credible intervals describing the probability distribution of coefficient values.  There are several ways to gain inference from these values.  The overlap of the credible interval with zero is often used to indicate the presence of a ‘significant effect’, though this approach has been criticised for several reasons(Gelman et al., 2012; Ogle et al., 2019).  An alternative option is to use Bayes factorsto quantify support for or against the hypothesis H1:β16=0 (Kass and Raftery, 1995;Jeffreys, 1998). In regression settings, the Bayes factor for a point null hypothesis can be calculated quickly using the Savage-Dickey method (Wagenmakers et al., 2010).The model shown here is intended as a minimal example of a simple linear regression with node-level effects to account for the non-independence of edges. The regression can be extended to support different family distributions, hierarchical effects, non-linearities,and more. We have included examples of diagnostics and some extensions in the example code.

#### simulated data ####
### create data list
simdat_ls <- list(
  n_dyads  = nrow(dyads),       # total number of times one or other of the dyad was observed
  together = dyads$event_count, # count number of sightings seen together
  apart    = dyads$apart,       # count number of sightings seen apart
  age_diff = dyads$age_diff)    # age difference between individuals
  #dem_dyad = dyads$dem_dyad_f)  # demographic pair of dyad

simdat_ls$together <- simdat_ls$together+1
simdat_ls$apart <- simdat_ls$apart+1

#test <- ulam(alist(events ~ dbeta(together, apart)),  data = simdat_ls, chains = 1, cores = 4)
test <- ulam(alist(
  events ~ dbeta(together, apart),
  weight ~ dnorm(mu, sigma),
  mu <- events + b1 * age_diff, # + u_i[id_1] + u_j[id_2], -- should include individual effects but no idea how to do this so that unaffected by which is id_1 and which is id_2
  b1 ~ dnorm(-0.5,2),
  sigma ~ dexp(1)),
  data = simdat_ls, chains = 1, cores = 4)
stancode(test)
precis(test)
traceplot(test)

d_reg_test <- cmdstan_model("models/dyadic_regression_hkm_22.03.03.stan")
d_reg_test

### Fit model
test2 <- d_reg_test$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
test2$summary()

#### real data ###
nodes$age_sex <- paste(nodes$age_cat_id ,nodes$sex, sep = '')
nodes$family[which(nodes$family == "M26 M87 Dead 1 Jul 17 Poached")] <- 'M26 M87'
nodes$family[which(nodes$family == 'n/a')] <- NA
nodes <- separate(nodes, family, into = c('family.1','family.2','family.3'), sep = ' ')

age_sex_family <- nodes[,c('id','age_sex','family.1','family.2','family.3')]
colnames(age_sex_family) <- c('id_1','age_sex_1','family_1.1','family_1.2','family_1.3')
counts_df <- left_join(x = counts_df, y = age_sex_family, by = 'id_1')
colnames(age_sex_family) <- c('id_2','age_sex_2','family_2.1','family_2.2','family_2.3')
counts_df <- left_join(x = counts_df, y = age_sex_family, by = 'id_2')
rm(age_sex_family)
counts_df$dem_dyad <- paste(counts_df$age_sex_1, counts_df$age_sex_2, sep = '_')
barplot(table(counts_df$dem_dyad), las = 1, horiz = T) # some have only a very few pairs

for(i in 1:nrow(counts_df)){
  counts_df$family_1.1[i] <- ifelse(is.na(counts_df$family_1.1[i]) == TRUE, 'NA', counts_df$family_1.1[i])
  counts_df$family_1.2[i] <- ifelse(is.na(counts_df$family_1.2[i]) == TRUE, 'NA', counts_df$family_1.2[i])
  counts_df$family_1.3[i] <- ifelse(is.na(counts_df$family_1.3[i]) == TRUE, 'NA', counts_df$family_1.3[i])
  counts_df$family_2.1[i] <- ifelse(is.na(counts_df$family_2.1[i]) == TRUE, 'NA', counts_df$family_2.1[i])
  counts_df$family_2.2[i] <- ifelse(is.na(counts_df$family_2.2[i]) == TRUE, 'NA', counts_df$family_2.2[i])
  counts_df$family_2.3[i] <- ifelse(is.na(counts_df$family_2.3[i]) == TRUE, 'NA', counts_df$family_2.2[i])
}

counts_df$family1 <- ifelse(counts_df$family_1.1 == counts_df$id_2, 'family', 
                            ifelse(counts_df$family_2.1 == counts_df$id_1, 'family', 'not_family'))
counts_df$family2 <- ifelse(counts_df$family_1.2 == counts_df$id_2, 'family',  
                            ifelse(counts_df$family_2.2 == counts_df$id_1, 'family', 'not_family'))
counts_df$family3 <- ifelse(counts_df$family_1.3 == counts_df$id_2, 'family', 
                            ifelse(counts_df$family_2.3 == counts_df$id_1, 'family', 'not_family'))

counts_df$family <- ifelse(counts_df$family1 == 'family' | counts_df$family2 == 'family' | counts_df$family3 == 'family',
                           'family','not_family')
table(counts_df$family)

family <- counts_df[counts_df$family == 'family',]
unique(family$dem_type)

counts_df$family <- ifelse(counts_df$family == 'family','mother_calf',counts_df$family)
counts_df$family <- ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'AF_PF',
                           'family',)




