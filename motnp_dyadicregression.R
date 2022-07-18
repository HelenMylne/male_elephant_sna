### Bayesian analysis of ALERT data
### Set up ####
# load packages
library(tidyverse)
library(dplyr)
library(rstan)
library(rethinking)
library(igraph)
library(dagitty)
library(cmdstanr)

# information
R.Version()
rstan::stan_version()

# set stan path
set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')

# set seed
set.seed(12345)

################ Create simulated data set ################
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

# clear environment
rm(population, i, N)

################ Dyadic regression -- ulam() ################
### create data list -- can contain no NA values in any column, even if column is not specified in model
simdat_ls <- list(
  n_dyads  = nrow(dyads),          # total number of times one or other of the dyad was observed
  together = dyads$event_count +1, # count number of sightings seen together -- can't work out how in ulam() to add the 1 because it's got a probelm with adding int+int[] again, so for the sake of testing have just adding it first
  apart    = dyads$apart +1,       # count number of sightings seen apart -- can't work out how in ulam() to add the 1 because it's got a probelm with adding int+int[] again, so for the sake of testing have just adding it first
  age_diff = dyads$age_diff)       # age difference between individuals

test <- ulam(alist(
  events ~ dbeta(together, apart),
  weight ~ dnorm(mu, sigma),
  mu <- events + b1 * age_diff, # + u_i[id_1] + u_j[id_2], -- should include individual effects but no idea how to do this so that unaffected by which is id_1 and which is id_2. Also I strongly suspect this is the wrong formula, but I couldn't see how weight and events linked together any other way than if the overall "SRI" became the intercept the age effect?
  b1 ~ dnorm(-0.5,2),
  sigma ~ dexp(1)),
  data = simdat_ls, chains = 1, cores = 4)
stancode(test)
precis(test)
traceplot(test)

d_reg_test <- cmdstan_model("models/dyadic_regression_hkm_22.03.03.stan")
d_reg_test

### create data list -- can contain no NA values in any column, even if column is not specified in model
simdat_ls <- list(
  n_dyads  = nrow(dyads),       # total number of times one or other of the dyad was observed
  together = dyads$event_count, # count number of sightings seen together
  apart    = dyads$apart,       # count number of sightings seen apart
  age_diff = dyads$age_diff)    # age difference between individuals

### Fit model
test2 <- d_reg_test$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4) # entirely fails to run -- generates NaN values in sigma, but not sure why


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
################ 2) Create simulated data set ################
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
dyads$node_1 <- as.integer(as.factor(dyads$id_1)) ; dyads$node_2 <- as.integer(as.factor(dyads$id_2))+1
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

### Assign edge weight to dyad pairs -- no effect of age
for(i in 1:nrow(dyads)){
  if(dyads$pair_type[i] == 'family')  {dyads$edge[i] <- rethinking::rbeta2(1,0.95,50)   # mother-calf: edge weight around 0.95
  } else 
    if(dyads$pair_type[i] == 'group') {dyads$edge[i] <- rethinking::rbeta2(1,0.8,50)    # same family, not mother-calf: edge ~ 0.8
    } else 
      dyads$edge[i] <- rethinking::rbeta2(1,0.20,50)                                    # different family: edge weight around 0.2
}

summary(dyads$edge)
boxplot(dyads$edge ~ dyads$pair_type)                                                   # check simulations are around what I expected
points(y = dyads$edge[which(dyads$pair_type == 'group')],
       x = runif(min = 1.6, max = 2.4, n = length(which(dyads$pair_type == 'group'))),  # had a few problems with this one, but looks good now
       col = col.alpha('black',0.3), pch = 19)

### Assign edge weight to dyad pairs -- dyads containing older individuals are stronger
dyads$edge_age <- ifelse(dyads$pair_type == 'family', dyads$edge*7*7, # retain high strength of mother-calf pairs, despite young age of calves
                         ifelse(dyads$pair_type == 'group',           # retain within-group > between-group, but make older individuals stronger
                                ifelse(dyads$age_1 > dyads$age_2,
                                       dyads$edge*7*dyads$age_1,
                                       dyads$edge*7*dyads$age_2),
                                dyads$edge*dyads$age_1*dyads$age_2))  # outside groups, pairs of older individuals strongest
dyads$edge_age <- dyads$edge_age/max(dyads$edge_age)

summary(dyads$edge_age)
boxplot(dyads$edge_age ~ dyads$pair_type)                                               # check simulations are around what I expected
points(y = dyads$edge_age[which(dyads$pair_type == 'group')],
       x = runif(min = 1.6, max = 2.4, n = length(which(dyads$pair_type == 'group'))),  # had a few problems with this one, but looks good now
       col = col.alpha('black',0.3), pch = 19)

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

################ 3) Import edge weight and nodal data ################
### model draws ####
#draws_motnp2.2 <- read_csv('data_processed/motnp_bayesian_edgedistributions_a2.b2_22.02.07.csv') %>%
#  data.matrix()

# summarise -- look for any anomalies in draw values or chain variation
summaries <- data.frame(dyad = colnames(draws_motnp2.2[,2:106954]),
                        min = rep(NA, ncol(draws_motnp2.2)-1),
                        max = rep(NA, ncol(draws_motnp2.2)-1),
                        mean = rep(NA, ncol(draws_motnp2.2)-1),
                        median = rep(NA, ncol(draws_motnp2.2)-1),
                        sd = rep(NA, ncol(draws_motnp2.2)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_motnp2.2[,i+1])
  summaries$max[i]    <- max(draws_motnp2.2[,i+1])
  summaries$mean[i]   <- mean(draws_motnp2.2[,i+1])
  summaries$median[i] <- median(draws_motnp2.2[,i+1])
  summaries$sd[i]     <- sd(draws_motnp2.2[,i+1])
}

### nodes data ####
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

# create age_sex variable
nodes$age_sex <- paste(nodes$age_cat_id ,nodes$sex, sep = '')
nodes$family[which(nodes$family == "M26 M87 Dead 1 Jul 17 Poached")] <- 'M26 M87'
nodes$family[which(nodes$family == 'n/a')] <- NA
nodes <- separate(nodes, family, into = c('family.1','family.2','family.3'), sep = ' ')

### dyads ####
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

# reassign dyad numbers to remove gaps
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1

# add in age_sex variable
age_sex_family <- nodes[,c('id','age_sex','family.1','family.2','family.3')]
colnames(age_sex_family) <- c('id_1','age_sex_1','family_1.1','family_1.2','family_1.3')
counts_df <- left_join(x = counts_df, y = age_sex_family, by = 'id_1')
colnames(age_sex_family) <- c('id_2','age_sex_2','family_2.1','family_2.2','family_2.3')
counts_df <- left_join(x = counts_df, y = age_sex_family, by = 'id_2')
rm(age_sex_family)
counts_df$dem_dyad <- paste(counts_df$age_sex_1, counts_df$age_sex_2, sep = '_')
barplot(table(counts_df$dem_dyad), las = 1, horiz = T) # some have only a very few pairs

# add in variable for family relationships
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
                           'mother_calf','not_family')
table(counts_df$family)

family <- counts_df[counts_df$family == 'mother_calf',]
unique(family$dem_type)

counts_df$family <- ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'AF_AF','family',
                           ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'AF_PF','family',
                                  ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'AF_PM','family',
                                         ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'JM_CM','family',
                                                ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'CM_CU','family',
                                                       ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'JM_CU','family',
                                                              ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'PM_CM','family',
                                                                     ifelse(counts_df$family == 'mother_calf' & counts_df$dem_type == 'JU_CU','family',counts_df$family))))))))

table(counts_df$family)

# add some summary data
counts_df <- left_join(counts_df, summaries, by = 'dyad')
colnames(counts_df)[2:3] <- c('id_1','id_2')

################ 5) Dyadic regression ################
# A common type of analysis in studies of social networks is to consider factors that affectedge weights.  This type of analysis is often called dyadic regression, where explanatoryfactors such as sex or age difference are regressed against edge weight. The edge weightsin dyadic regression are non-independent, as variables such as age difference or sex dif-ference are inherently linked to individual-level attributes.  This means that effects dueto age or sex in individualiaffect all dyads that connecti.  This non-independence canbe controlled for by including node-level effects in the regression (Tranmer et al., 2014).Using the count data example discussed earlier, we propose to conduct dyadic regres-sion using a standard regression model of the form:
#         X_ij ∼ Poisson(lambda_ij)
#         log(lambda_ij) = weight_ij + Location_ij
#         weight_ij ∼ Normal(beta_0 + beta_1*x_ij + u_i + u_j, sigma_sq)
# where beta_0 is the intercept parameter, beta_1 is the slope parameter, ui and uj are parameters accounting for the effect of node membership on edge weight, and sigma is the standard deviation of the residuals.  The dashed line indicates that the model can be fit in two separate parts:the first part being the core edge weight model, and the second part being the dyadic regression model.  It is not strictly necessary to separate out the model like this, but do-ing so makes it possible to fit the edge weight model once and conduct multiple types of analysis on it afterwards without fitting the entire model again.Theβ0andβ1parameters can be interpreted in the same way as usual regression coefficients, and will be accompanied with credible intervals describing the probability distribution of coefficient values.  There are several ways to gain inference from these values.  The overlap of the credible interval with zero is often used to indicate the presence of a ‘significant effect’, though this approach has been criticised for several reasons(Gelman et al., 2012; Ogle et al., 2019).  An alternative option is to use Bayes factorsto quantify support for or against the hypothesis H1:β16=0 (Kass and Raftery, 1995;Jeffreys, 1998). In regression settings, the Bayes factor for a point null hypothesis can be calculated quickly using the Savage-Dickey method (Wagenmakers et al., 2010).The model shown here is intended as a minimal example of a simple linear regression with node-level effects to account for the non-independence of edges. The regression can be extended to support different family distributions, hierarchical effects, non-linearities,and more. We have included examples of diagnostics and some extensions in the example code.
adults <- counts_df[counts_df$age_cat_id_1 > 3 & counts_df$age_cat_id_2 > 3,]
ggplot(adults, aes(x = all_events, y = median, col = dem_type))+
  geom_point()+
  scale_x_continuous(limits = c(0,20))+
  scale_y_continuous(limits = c(0,0.8))

ggplot(adults, aes(x = all_events, y = median, col = dem_type))+
  geom_point()+
  #scale_x_continuous(limits = c(0,20))+
  #scale_y_continuous(limits = c(0,0.8))+
  facet_wrap(. ~ dem_type, scales = 'free')


#### simulated data ####
### create data list
simdat_ls <- list(
  n_dyads  = nrow(dyads),       # total number of times one or other of the dyad was observed
  together = dyads$event_count, # count number of sightings seen together
  apart    = dyads$apart,       # count number of sightings seen apart
  age_diff = dyads$age_diff)    # age difference between individuals
19:10
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


