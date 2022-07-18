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
draws_motnp1.1 <- read_csv('data_processed/motnp_bayesian_edgedistributions_a1.b1_22.03.03.csv') %>%
  data.matrix()

# summarise -- look for any anomalies in draw values or chain variation
summaries <- data.frame(dyad = colnames(draws_motnp1.1[,2:106954]),
                        min = rep(NA, ncol(draws_motnp1.1)-1),
                        max = rep(NA, ncol(draws_motnp1.1)-1),
                        mean = rep(NA, ncol(draws_motnp1.1)-1),
                        median = rep(NA, ncol(draws_motnp1.1)-1),
                        sd = rep(NA, ncol(draws_motnp1.1)-1))
for(i in 1:nrow(summaries)){
  summaries$min[i]    <- min(draws_motnp1.1[,i+1])
  summaries$max[i]    <- max(draws_motnp1.1[,i+1])
  summaries$mean[i]   <- mean(draws_motnp1.1[,i+1])
  summaries$median[i] <- median(draws_motnp1.1[,i+1])
  summaries$sd[i]     <- sd(draws_motnp1.1[,i+1])
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
colnames(counts_df)

################ 4) Nodal regression ################
# The final common type of network analysis we will cover here is nodal regression, where a regression is performed to analyse the relationship between a nodal network metric(such as centrality) and nodal traits (such as age and sex).  These analyses are usually used to assess how network position depends on various biological factors, but can also be used where network position is a predictor.  Since node metrics are derivative measures of the network, uncertainty in edge weights should ideally propagate through to node metrics, and on to coefficients in regression analyses, giving an accurate estimate of the total uncertainty in inferred parameters.  The core of the nodal regression is similar to dyadic regression, taking the form:
#X(n)i j∼Poisson(λ(n)i jD(n)i j)
#logÄλ(n)i jä=ωi j+L(n)i j
#Mi(ω)∼Normal(β0+β1xi,σ2)
# where β0 is the intercept parameter, β1 is the slope parameter, and σ is the standard deviation of the residuals. Mi(ω) denotes the node metric estimates when applied to the edge weights estimated in the top part of the model. Calculating node metrics within the model may present a practical challenge when using standard model fitting software, as many node metrics can not be described purely in terms of closed-form equations (Freeman, 1978).  In this case, splitting up the model along the dashed line becomes an important step, as it allows the edge weights to be estimated using a standard piece of software, and leaves the regression part of the model to be fit using a sampler that supports numeric estimates of network centralities. As part of the code we have made available, we have written a custom Metropolis-Hastings sampler that allows numeric estimates of node metrics to be used when fitting the model. Importantly, our sampler samples from the joint distribution of edge weights, rather than sampling edge weights independently. This ensures that the effect of any structure in the observation data (such as location effects) is maintained and propagated through to the node metric estimates, and subsequently to regression coefficient estimates.
### simulated data  -- obtain model draws####
# create data list
simdat_ls <- list(
  n_dyads  = nrow(dyads),        # Number of dyads
  together = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart         # Number of sightings of each dyad apart
)

### Fit model
mod_1.1 <- cmdstan_model("models/simpleBetaNet_DWF_22.01.23.stan")
edge_weight_1.1 <- mod_1.1$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
output1 <- as.data.frame(read_cmdstan_csv(edge_weight_1.1$output_files()[1])$post_warmup_draws)
output2 <- as.data.frame(read_cmdstan_csv(edge_weight_1.1$output_files()[2])$post_warmup_draws)
output3 <- as.data.frame(read_cmdstan_csv(edge_weight_1.1$output_files()[3])$post_warmup_draws)
output4 <- as.data.frame(read_cmdstan_csv(edge_weight_1.1$output_files()[4])$post_warmup_draws)

draws_simdat_1.1 <- rbind(output1, output2, output3, output4) %>% 
  data.matrix()
draws_simdat_1.1 <- draws_simdat_1.1[,2:7141]
rm(output1, output2, output3, output4)

### simulated data  -- create centrality matrix -- degree centrality for each individual at every draw from posterior distribution ####
sim_tensor <- array(0, c(NROW(unique(dyads$id_1))+1,
                         NROW(unique(dyads$id_2))+1,
                         NROW(draws_simdat_1.1)),
                    dimnames = list(c(sort(unique(dyads$id_1)),'U9'),
                                    c(sort(unique(dyads$id_1)),'U9'),
                                    NULL))
N <- nrow(dyads)

for (i in 1:N) {
  dyad_row <- dyads[i, ]
  sim_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_simdat_1.1[, i]
}
sim_tensor[,,1]

centrality_matrix <- array(0, c(nrow = NROW(draws_simdat_1.1), ncol = length(unique(dyads$id_1))+1, 3))
for (i in 1:NROW(draws_simdat_1.1)) {
  g <- graph_from_adjacency_matrix(sim_tensor[, , i], mode="undirected", weighted=TRUE)
  centrality_matrix[i, , 1] <- strength(g)
  centrality_matrix[i, , 2] <- betweenness(g)
  centrality_matrix[i, , 3] <- eigen_centrality(g)$vector
}

colnames(centrality_matrix) <- c(sort(unique(dyads$id_1)),unique(dyads$id_2)[119])
head(centrality_matrix)

strength_quantiles <- t(apply(centrality_matrix[,,1], 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
strength_quantiles
btwn_quantiles <- t(apply(centrality_matrix[,,2], 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
btwn_quantiles
eigen_quantiles <- t(apply(centrality_matrix[,,3], 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
eigen_quantiles

cq <- data.frame(id = c(sort(unique(dyads$id_1)),'U9'),
                 strength_lwr = strength_quantiles[,1],
                 strength_mid = strength_quantiles[,2],
                 strength_upr = strength_quantiles[,3],
                 btwn_lwr = btwn_quantiles[,1],
                 btwn_mid = btwn_quantiles[,2],
                 btwn_upr = btwn_quantiles[,3],
                 eigen_lwr = eigen_quantiles[,1],
                 eigen_mid = eigen_quantiles[,2],
                 eigen_upr = eigen_quantiles[,3])

population <- left_join(population, cq, by = 'id')
head(population)

boxplot(population$strength_mid ~ population$sex, notch = T)
boxplot(population$strength_mid ~ population$age, notch = T)

boxplot(population$btwn_mid ~ population$sex, notch = T)
boxplot(population$btwn_mid ~ population$age, notch = T)

boxplot(population$eigen_mid ~ population$sex, notch = T)
boxplot(population$eigen_mid ~ population$age, notch = T)

population$age_jit <- jitter(as.numeric(population$age),1.5)
plot(population$strength_mid ~ population$age_jit, las =  1, pch = 21, col = 'black', ylim = c(20, 40),
     bg = ifelse(population$sex == 'F', 'magenta',
                 ifelse(population$sex == 'M', 'blue', 'yellow')),
     xlab = 'age category (jittered)', ylab = 'network centrality')
for(i in 1:nrow(population)){
  lines(c(population$age_jit[i], population$age_jit[i]),
        c(population$strength_lwr[i], population$strength_upr[i]),
        col = ifelse(population$sex[i] == 'F', 'magenta',
                     ifelse(population$sex[i] == 'M', 'blue', 'yellow')),
        lwd = 1, lty = 1)
}
legend(x = 1, y = 25, legend = c('male','female','unknown sex'),
       pt.bg = c('blue','magenta','yellow'), pch = 21)

plot(population$btwn_mid ~ population$age_jit, las =  1, pch = 21, col = 'black', ylim = c(0, 1400),
     bg = ifelse(population$sex == 'F', 'magenta',
                 ifelse(population$sex == 'M', 'blue', 'yellow')),
     xlab = 'age category (jittered)', ylab = 'network centrality')
for(i in 1:nrow(population)){
  lines(c(population$age_jit[i], population$age_jit[i]),
        c(population$btwn_lwr[i], population$btwn_upr[i]),
        col = ifelse(population$sex[i] == 'F', 'magenta',
                     ifelse(population$sex[i] == 'M', 'blue', 'yellow')),
        lwd = 1, lty = 1)
}
legend(x = 0.8, y = 1350, legend = c('male','female','unknown sex'),
       pt.bg = c('blue','magenta','yellow'), pch = 21)

plot(population$eigen_mid ~ population$age_jit, las =  1, pch = 21, col = 'black', ylim = c(0.5, 1),
     bg = ifelse(population$sex == 'F', 'magenta',
                 ifelse(population$sex == 'M', 'blue', 'yellow')),
     xlab = 'age category (jittered)', ylab = 'network centrality')
for(i in 1:nrow(population)){
  lines(c(population$age_jit[i], population$age_jit[i]),
        c(population$eigen_lwr[i], population$eigen_upr[i]),
        col = ifelse(population$sex[i] == 'F', 'magenta',
                     ifelse(population$sex[i] == 'M', 'blue', 'yellow')),
        lwd = 1, lty = 1)
}
legend(x = 1, y = 0.65, legend = c('male','female','unknown sex'),
       pt.bg = c('blue','magenta','yellow'), pch = 21)

population_tidy <- pivot_longer(data = population, cols = 8:16,
                           names_to = 'centrality_type', values_to = 'centrality')
population_tidy$centrality_type <- factor(population_tidy$centrality_type,
                                     levels = c('strength_lwr', 'strength_mid', 'strength_upr',
                                                'btwn_lwr', 'btwn_mid', 'btwn_upr',
                                                'eigen_lwr', 'eigen_mid', 'eigen_upr'))
ggplot(population_tidy, aes(y = centrality, x = as.factor(age), fill = centrality_type))+
  geom_boxplot(notch = T)+
  facet_wrap(.~centrality_type, scales = 'free')+ # lwr, mid and upr = very similar patterns. degree = very variable depending on level
  theme_light()+
  theme(legend.position = 'none',
        strip.text = element_text(colour = 'black'))

# Since network strength is strictly positive, a Gaussian error is not a reasonable model for the data. The Gaussian family model is much easier to implement as well as interpret than many other models, so we will standardise the centralities by taking z-scores.
centrality_matrix_std <- (centrality_matrix - apply(centrality_matrix, 1, mean))/apply(centrality_matrix, 1, sd)

### define functions and predictor matrix for running model ####
# The key challenge to quantifying uncertainty in network analysis is to incorporate uncertainty due to sampling into downstream analyses, commonly regression. This can be achieved by modifying the likelihood function of a regression model to treat the network centralities with uncertainty. We have written a custom MCMC sampler function that samples from the joint distribution of network centralities calculated earlier and treats those samples as the data in the likelihood function. Likelihood functions for the sampler use the index variable to keep track of which data points are being compared internally in the sampler, to ensure that candidate steps in the MCMC are not accepted or rejected because they are being compared to different data points, rather than because of the parameter space. Custom likelihood functions take a similar form to the target += syntax in Stan, but for more specific resources the following document is a good start: https://www.ime.unicamp.br/~cnaber/optim_1.pdf.

# define metropolis sampler function
metropolis <- function(target, initial, iterations=10000, warmup=iterations/2, thin=4, chain_id=1, refresh=2000) {
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

## loglik function for alternative node age categories
loglik <- function(params, Y, X, index) {
  # Define parameters
  intercept <- params[1]
  beta_age1 <- params[2]
  beta_age2 <- params[3]
  beta_age3 <- params[4]
  beta_age4 <- params[5]
  beta_age5 <- params[6]
  beta_age6 <- params[7]
  beta_age7 <- params[8]
  beta_females <- params[9]
  beta_males <- params[10]
  beta_unks <- params[11]
  sigma <- exp(params[12]) # Exponential keeps underlying value unconstrained, which is much easier for the sampler.
  
  # Sample data according to index
  y <- Y[index %% dim(Y)[1] + 1, ]
  
  # Define model
  target <- 0
  target <- target + sum(dnorm(y, mean = intercept + beta_age1*X[,1] + beta_age2*X[,2] + beta_age3*X[,3] + beta_age4*X[,4] + beta_age5*X[,5] + beta_age1*X[,6] + beta_age7*X[,7] + beta_females*X[,8] + beta_males*X[,9] + beta_unks*X[,10],
                               sd = sigma, log = TRUE))                   # Main model
  target <- target + dnorm(intercept,    mean = 0, sd = 0.5, log = TRUE)  # Prior on intercept
  target <- target + dnorm(beta_age1,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age2,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age3,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age4,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age5,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age6,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age7,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_females, mean = 0, sd = 0.5, log = TRUE)  # Prior on female coefficient
  target <- target + dnorm(beta_males,   mean = 0, sd = 0.5, log = TRUE)  # Prior on male coefficient
  target <- target + dnorm(beta_unks,    mean = 0, sd = 0.5, log = TRUE)  # Prior on unknown coefficient
  target <- target + dexp(sigma, 1, log = TRUE)                           # Prior on sigma
  
  return(target)
}

# predictor matrix
predictor_matrix <- matrix(0, nrow = length(population$id), ncol = 10)
colnames(predictor_matrix) <- c('age1','age2','age3','age4','age5','age6','age7',
                                'female', 'male','unknown')
predictor_matrix[which(population$age == 1), 1] <- 1
predictor_matrix[which(population$age == 2), 2] <- 1
predictor_matrix[which(population$age == 3), 3] <- 1
predictor_matrix[which(population$age == 4), 4] <- 1
predictor_matrix[which(population$age == 5), 5] <- 1
predictor_matrix[which(population$age == 6), 6] <- 1
predictor_matrix[which(population$age == 7), 7] <- 1
predictor_matrix[which(population$sex == 'F'), 8] <- 1
predictor_matrix[which(population$sex == 'M'), 9] <- 1
predictor_matrix[which(population$sex == 'U'), 10] <- 1
predictor_matrix

target <- function(params, index) loglik(params, centrality_matrix_std[,,1], predictor_matrix, index)

## run model ####
chain <- metropolis(target, initial = rep(0,12), iterations=100000, thin=100, refresh=10000)

colnames(chain) <- c("intercept",
                     "beta_age1", "beta_age2", "beta_age3", "beta_age4", "beta_age5", "beta_age6", "beta_age7", 
                     "beta_females", "beta_males", "beta_unks", "sigma")
head(chain)

par(mfrow=c(4,3))
for (i in 1:12) {
  plot(chain[, i], type="l", 
       ylab = colnames(chain)[i])
}

par(mfrow = c(1,1))
plot(density(centrality_matrix_std[1, ,1]), ylim=c(0,40), xlim = c(-0.5,0.1),
     main="", xlab="Standardised network strength")
sample_ids <- sample(1:1000, size=200)
for (i in sample_ids) {
  pred <- rnorm(8, mean = chain[i, "intercept"] + 
                  chain[i,"beta_age1"]*predictor_matrix[,1] + 
                  chain[i,"beta_age2"]*predictor_matrix[,2] + 
                  chain[i,"beta_age3"]*predictor_matrix[,3] + 
                  chain[i,"beta_age4"]*predictor_matrix[,4] + 
                  chain[i,"beta_age5"]*predictor_matrix[,5] + 
                  chain[i,"beta_age6"]*predictor_matrix[,6] + 
                  chain[i,"beta_age7"]*predictor_matrix[,7] + 
                  chain[i,"beta_females"]*predictor_matrix[,8] + 
                  chain[i,"beta_males"]*predictor_matrix[,9] + 
                  chain[i,'beta_unks']*predictor_matrix[,10],
                sd = exp(chain[i, "sigma"]))
  lines(density(centrality_matrix_std[i, ,1]), col=rgb(0, 0, 0, 0.25))
  lines(density(pred), col=rgb(0, 0, 1, 0.25))
}

coefficient_quantiles <- t(apply(chain, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
coefficient_quantiles

beta_difference <- chain[, "beta_females"] - chain[, "beta_males"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

beta_difference <- chain[, "beta_females"] - chain[, "beta_unks"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

beta_difference <- chain[, "beta_males"] - chain[, "beta_unks"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

## adjust loglik function to remove sex effect ####
loglik <- function(params, Y, X, index) {
  # Define parameters
  intercept <- params[1]
  beta_age1 <- params[2]
  beta_age2 <- params[3]
  beta_age3 <- params[4]
  beta_age4 <- params[5]
  beta_age5 <- params[6]
  beta_age6 <- params[7]
  beta_age7 <- params[8]
  sigma <- exp(params[9]) # Exponential keeps underlying value unconstrained, which is much easier for the sampler.
  
  # Sample data according to index
  y <- Y[index %% dim(Y)[1] + 1, ]
  
  # Define model
  target <- 0
  target <- target + sum(dnorm(y, mean = intercept + beta_age1*X[,1] + beta_age2*X[,2] + beta_age3*X[,3] + beta_age4*X[,4] + beta_age5*X[,5] + beta_age1*X[,6] + beta_age7*X[,7],
                               sd = sigma, log = TRUE))                   # Main model
  target <- target + dnorm(intercept,    mean = 0, sd = 0.5, log = TRUE)  # Prior on intercept
  target <- target + dnorm(beta_age1,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age2,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age3,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age4,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age5,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age6,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_age7,    mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dexp(sigma, 1, log = TRUE)                           # Prior on sigma
  return(target)
}

# predictor matrix
predictor_matrix <- matrix(0, nrow = length(population$id), ncol = 7)
colnames(predictor_matrix) <- c('age1','age2','age3','age4','age5','age6','age7')
predictor_matrix[which(population$age == 1), 1] <- 1
predictor_matrix[which(population$age == 2), 2] <- 1
predictor_matrix[which(population$age == 3), 3] <- 1
predictor_matrix[which(population$age == 4), 4] <- 1
predictor_matrix[which(population$age == 5), 5] <- 1
predictor_matrix[which(population$age == 6), 6] <- 1
predictor_matrix[which(population$age == 7), 7] <- 1
predictor_matrix

target <- function(params, index) loglik(params, centrality_matrix_std[,,1], predictor_matrix, index)

### run model ####
chain <- metropolis(target, initial = rep(0,9), iterations=1000000, thin=100, refresh=10000)

colnames(chain) <- c("intercept",
                     "beta_age1", "beta_age2", "beta_age3", "beta_age4", "beta_age5", "beta_age6", "beta_age7", 
                     "sigma")
head(chain)

par(mfrow=c(3,3))
for (i in 1:9) {
  plot(chain[, i], type="l", 
       ylab = colnames(chain)[i])
}

par(mfrow = c(1,1))
plot(density(centrality_matrix_std[1, ,1]), ylim=c(0,40), xlim = c(-0.5,0.1),
     main="", xlab="Standardised network strength")
sample_ids <- sample(1:1000, size=200)
for (i in sample_ids) {
  pred <- rnorm(8, mean = chain[i, "intercept"] + 
                  chain[i,"beta_age1"]*predictor_matrix[,1] + 
                  chain[i,"beta_age2"]*predictor_matrix[,2] + 
                  chain[i,"beta_age3"]*predictor_matrix[,3] + 
                  chain[i,"beta_age4"]*predictor_matrix[,4] + 
                  chain[i,"beta_age5"]*predictor_matrix[,5] + 
                  chain[i,"beta_age6"]*predictor_matrix[,6] + 
                  chain[i,"beta_age7"]*predictor_matrix[,7] + 
                  chain[i,"beta_females"]*predictor_matrix[,8] + 
                  chain[i,"beta_males"]*predictor_matrix[,9] + 
                  chain[i,'beta_unks']*predictor_matrix[,10],
                sd = exp(chain[i, "sigma"]))
  lines(density(centrality_matrix_std[i, ,1]), col=rgb(0, 0, 0, 0.25))
  lines(density(pred), col=rgb(0, 0, 1, 0.25))
}

coefficient_quantiles <- t(apply(chain, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
coefficient_quantiles

beta_difference <- chain[, "beta_females"] - chain[, "beta_males"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

beta_difference <- chain[, "beta_females"] - chain[, "beta_unks"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

beta_difference <- chain[, "beta_males"] - chain[, "beta_unks"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

# repeat for betweenness
target <- function(params, index) loglik(params, centrality_matrix_std[,,2], predictor_matrix, index)

# repeat for eigenvector
target <- function(params, index) loglik(params, centrality_matrix_std[,,3], predictor_matrix, index)


### create centrality matrix -- degree centrality for each individual at every draw from posterior distribution ####
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

centrality_matrix <- matrix(0, nrow = NROW(draws_motnp1.1), ncol = length(unique(counts_df$id_1))+1)
for (i in 1:NROW(draws_motnp1.1)) {
  g <- graph_from_adjacency_matrix(adj_tensor[, , i], mode="undirected", weighted=TRUE)
  centrality_matrix[i, ] <- strength(g)
}

colnames(centrality_matrix) <- c(sort(unique(counts_df$id_1)),unique(counts_df$id_2)[462])
head(centrality_matrix)
write_csv(as.data.frame(centrality_matrix), 'data_processed/motnp_bayesian_nodalregression_centralitymatrix_22.03.14.csv')

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
boxplot(nodes$strength_mid ~ nodes$age_cat_id, notch = T, horizontal = T, las = 1, ylab = '')

ggplot(nodes, aes(y = strength_mid, x = age_cat_id))+
  geom_jitter()+
  theme_light()

nodes$age <- jitter(as.numeric(nodes$age_cat_id),1.5)
plot(nodes$strength_mid ~ nodes$age, las =  1, pch = 21, col = 'black', ylim = c(10, 75),
     bg = ifelse(nodes$sex == 'F', 'magenta',
                 ifelse(nodes$sex == 'M', 'blue', 'yellow')),
     xlab = 'age category (jittered)', ylab = 'network centrality')
for(i in 1:nrow(nodes)){
  lines(c(nodes$age[i], nodes$age[i]),
        c(nodes$strength_lwr[i], nodes$strength_upr[i]),
        col = ifelse(nodes$sex[i] == 'F', 'magenta',
                     ifelse(nodes$sex[i] == 'M', 'blue', 'yellow')),
        lwd = 1, lty = 1)
}
legend(x = 1.5, y = 22, legend = c('male','female','unknown sex'),
       pt.bg = c('blue','magenta','yellow'), pch = 21)

nodes_tidy <- pivot_longer(data = nodes, cols = 15:22,
                           names_to = 'centrality_type', values_to = 'centrality')
nodes_tidy$centrality_type <- factor(nodes_tidy$centrality_type,
                                     levels = c('strength_lwr', 'strength_mid', 'strength_upr',
                                                'degree_0.1','degree_0.2','degree_0.3',
                                                'degree_0.4','degree_0.5'))
ggplot(nodes_tidy, aes(y = centrality, x = age_cat_id, fill = centrality_type))+
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
write_csv(nodes, 'data_processed/motnp_elenodes_centrality_22.03.14.csv')

nodes_tidy$centrality_type <- factor(nodes_tidy$centrality_type,
                                     levels = c('strength_lwr', 'strength_mid', 'strength_upr',
                                                'degree_0.1','degree_0.3','degree_0.5',
                                                'degree_0.2','degree_0.4'))
ggplot(nodes_tidy, aes(y = centrality, x = age_cat_id, fill = centrality_type))+
  geom_boxplot()+
  facet_wrap(.~centrality_type, scales = 'free')+ # lwr, mid and upr = very similar patterns. degree = very variable depending on level
  theme_light()+
  theme(legend.position = 'none',
        strip.text = element_text(colour = 'black'))

# all centrality measures ####
centrality_matrix <- array(0, c(nrow = NROW(draws_motnp1.1), ncol = length(unique(counts_df$id_1))+1, 3))
for (i in 1:NROW(draws_motnp1.1)) {
  g <- graph_from_adjacency_matrix(adj_tensor[, , i], mode="undirected", weighted=TRUE)
  centrality_matrix[i, , 1] <- strength(g)
  centrality_matrix[i, , 2] <- betweenness(g)
  centrality_matrix[i, , 3] <- eigen_centrality(g)$vector
}

colnames(centrality_matrix) <- c(sort(unique(counts_df$id_1)),unique(counts_df$id_2)[462])
head(centrality_matrix)

strength_quantiles <- t(apply(centrality_matrix[,,1], 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
strength_quantiles
btwn_quantiles <- t(apply(centrality_matrix[,,2], 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
btwn_quantiles
eigen_quantiles <- t(apply(centrality_matrix[,,3], 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
eigen_quantiles

cq <- data.frame(id = c(sort(unique(counts_df$id_1)),'U9'),
                 strength_lwr = strength_quantiles[,1],
                 strength_mid = strength_quantiles[,2],
                 strength_upr = strength_quantiles[,3],
                 btwn_lwr = btwn_quantiles[,1],
                 btwn_mid = btwn_quantiles[,2],
                 btwn_upr = btwn_quantiles[,3],
                 eigen_lwr = eigen_quantiles[,1],
                 eigen_mid = eigen_quantiles[,2],
                 eigen_upr = eigen_quantiles[,3])

nodes <- left_join(nodes, cq, by = 'id')
head(nodes)

boxplot(nodes$strength_mid ~ nodes$sex, notch = T)
boxplot(nodes$strength_mid ~ nodes$age_cat_id, notch = T)

boxplot(nodes$btwn_mid ~ nodes$sex, notch = T)
boxplot(nodes$btwn_mid ~ nodes$age_cat_id, notch = T)

boxplot(nodes$eigen_mid ~ nodes$sex, notch = T)
boxplot(nodes$eigen_mid ~ nodes$age_cat_id, notch = T)

nodes$age_jit <- jitter(as.numeric(nodes$age_cat_id),1.5)
plot(nodes$strength_mid ~ nodes$age_jit, las =  1, pch = 21, col = 'black', ylim = c(10, 80),
     bg = ifelse(nodes$sex == 'F', 'magenta',
                 ifelse(nodes$sex == 'M', 'blue', 'yellow')),
     xlab = 'age category (jittered)', ylab = 'degree centrality')
for(i in 1:nrow(nodes)){
  lines(c(nodes$age_jit[i], nodes$age_jit[i]),
        c(nodes$strength_lwr[i], nodes$strength_upr[i]),
        col = ifelse(nodes$sex[i] == 'F', 'magenta',
                     ifelse(nodes$sex[i] == 'M', 'blue', 'yellow')),
        lwd = 1, lty = 1)
}
legend(x = 1.4, y = 25, legend = c('male','female','unknown'),
       pt.bg = c('blue','magenta','yellow'), pch = 21, bg = 'transparent')
text('a)', x = 1, y = 79, cex = 3)

plot(nodes$btwn_mid ~ nodes$age_jit, las =  1, pch = 21, col = 'black', ylim = c(0, 15000),
     bg = ifelse(nodes$sex == 'F', 'magenta',
                 ifelse(nodes$sex == 'M', 'blue', 'yellow')),
     xlab = 'age category (jittered)', ylab = 'betweenness centrality')
for(i in 1:nrow(nodes)){
  lines(c(nodes$age_jit[i], nodes$age_jit[i]),
        c(nodes$btwn_lwr[i], nodes$btwn_upr[i]),
        col = ifelse(nodes$sex[i] == 'F', 'magenta',
                     ifelse(nodes$sex[i] == 'M', 'blue', 'yellow')),
        lwd = 1, lty = 1)
}
legend(x = 1.6, y = 14000, legend = c('male','female','unknown'),
       pt.bg = c('blue','magenta','yellow'), pch = 21)
text('b)', x = 1, y = 14700, cex = 3)

plot(nodes$eigen_mid ~ nodes$age_jit, las =  1, pch = 21, col = 'black', ylim = c(0, 1),
     bg = ifelse(nodes$sex == 'F', 'magenta',
                 ifelse(nodes$sex == 'M', 'blue', 'yellow')),
     xlab = 'age category (jittered)', ylab = 'eigenvector centrality')
for(i in 1:nrow(nodes)){
  lines(c(nodes$age_jit[i], nodes$age_jit[i]),
        c(nodes$eigen_lwr[i], nodes$eigen_upr[i]),
        col = ifelse(nodes$sex[i] == 'F', 'magenta',
                     ifelse(nodes$sex[i] == 'M', 'blue', 'yellow')),
        lwd = 1, lty = 1)
}
legend(x = 1.4, y = 0.2, legend = c('male','female','unknown'),
       pt.bg = c('blue','magenta','yellow'), pch = 21)
text('c)', x = 1, y = 0.98, cex = 3)

nodes_tidy <- pivot_longer(data = nodes, cols = 8:16,
                                names_to = 'centrality_type', values_to = 'centrality')
nodes_tidy$centrality_type <- factor(nodes_tidy$centrality_type,
                                          levels = c('strength_lwr', 'strength_mid', 'strength_upr',
                                                     'btwn_lwr', 'btwn_mid', 'btwn_upr',
                                                     'eigen_lwr', 'eigen_mid', 'eigen_upr'))
ggplot(nodes_tidy, aes(y = centrality, x = as.factor(age), fill = centrality_type))+
  geom_boxplot(notch = T)+
  facet_wrap(.~centrality_type, scales = 'free')+ # lwr, mid and upr = very similar patterns. degree = very variable depending on level
  theme_light()+
  theme(legend.position = 'none',
        strip.text = element_text(colour = 'black'))


### define functions and predictor matrix for running model ####
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
  target <- target + dnorm(intercept,    mean = 0, sd = 0.5, log = TRUE)  # Prior on intercept
  target <- target + dnorm(beta_age,     mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_females, mean = 0, sd = 0.5, log = TRUE)  # Prior on female coefficient
  target <- target + dnorm(beta_males,   mean = 0, sd = 0.5, log = TRUE)  # Prior on male coefficient
  target <- target + dnorm(beta_unks,    mean = 0, sd = 0.5, log = TRUE)  # Prior on unknown coefficient
  target <- target + dexp(sigma, 1, log = TRUE)                           # Prior on sigma
  
  return(target)
}

# Now we will prepare data for fitting the model. The predictor matrix is simply a matrix with 3 columns and 463 rows, corresponding to whether each of the 8 nodes is a female (column 1), a male (column 2), or an unknown (column 3).
predictor_matrix <- matrix(0, nrow = length(nodes$id), ncol = 4)
colnames(predictor_matrix) <- c('age','female', 'male','unknown')
predictor_matrix[, 1] <- as.numeric(nodes$age_cat_id)
predictor_matrix[which(nodes$sex == 'F'), 2] <- 1
predictor_matrix[which(nodes$sex == 'M'), 3] <- 1
predictor_matrix[which(nodes$sex == 'U'), 4] <- 1
predictor_matrix

# Since network strength is strictly positive, a Gaussian error is not a reasonable model for the data. The Gaussian family model is much easier to implement as well as interpret than many other models, so we will standardise the centralities by taking z-scores.
centrality_matrix_std <- (centrality_matrix - apply(centrality_matrix, 1, mean))/apply(centrality_matrix, 1, sd)

# Now we’re in a position to fit the model. To do this, we define the target function, which is simply a function that maps candidate parameters and a network centrality index to the log-likelihood of that function for the given sample of the centrality posterior. This means the target function can be written as a function of the data centrality_matrix_std and predictor_matrix.
target <- function(params, index) loglik(params, centrality_matrix_std, predictor_matrix, index)

### run model ####
# The function metropolis from sampler.R can now be used to fit the model using the provided target function, an initial set of parameters, and some additional MCMC options.
chain <- metropolis(target, initial = c(0, 0, 0, 0, 0, 0), iterations=100000, thin=100, refresh=10000)

colnames(chain) <- c("intercept", "beta_age", "beta_females", "beta_males", "beta_unks", "sigma")
head(chain)

par(mfrow=c(3,2))
for (i in 1:6) {
  plot(chain[, i], type="l", las  = 1,
       ylab = colnames(chain)[i])
}

par(mfrow = c(1,1))
plot(density(centrality_matrix_std[1, ]), ylim=c(0, 1.5), main="", xlab="Standardised network strength")
sample_ids <- sample(1:1000, size=200)
for (i in sample_ids) {
  pred <- rnorm(8, mean=chain[i, "intercept"] + chain[i,"beta_age"]*predictor_matrix[,1] + chain[i,"beta_females"]*predictor_matrix[,2] + chain[i,"beta_males"]*predictor_matrix[,3] + chain[i,'beta_unks']*predictor_matrix[,4], sd=exp(chain[i, "sigma"]))
  lines(density(centrality_matrix_std[i, ]), col=rgb(0, 0, 0, 0.25))
  lines(density(pred), col=rgb(0, 0, 1, 0.25))
}

coefficient_quantiles <- t(apply(chain, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
coefficient_quantiles

beta_difference <- chain[, "beta_females"] - chain[, "beta_males"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

beta_difference <- chain[, "beta_females"] - chain[, "beta_unks"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

beta_difference <- chain[, "beta_males"] - chain[, "beta_unks"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

### define new model ####
## adjust loglik function for alternative node age categories
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
  target <- target + dnorm(intercept,    mean = 0, sd = 0.5, log = TRUE)  # Prior on intercept
  target <- target + dnorm(beta_age,     mean = 0, sd = 0.5, log = TRUE)  # Prior on age coefficient
  target <- target + dnorm(beta_females, mean = 0, sd = 0.5, log = TRUE)  # Prior on female coefficient
  target <- target + dnorm(beta_males,   mean = 0, sd = 0.5, log = TRUE)  # Prior on male coefficient
  target <- target + dnorm(beta_unks,    mean = 0, sd = 0.5, log = TRUE)  # Prior on unknown coefficient
  target <- target + dexp(sigma, 1, log = TRUE)                           # Prior on sigma
  
  return(target)
}

# predictor matrix
predictor_matrix <- matrix(0, nrow = length(nodes$id), ncol = 10)
colnames(predictor_matrix) <- c('age1','age2','age3','age4','age5','age6','age7',
                                'female', 'male','unknown')
predictor_matrix[which(nodes$age_cat_id == '1'), 1] <- 1
predictor_matrix[which(nodes$age_cat_id == '2'), 2] <- 1
predictor_matrix[which(nodes$age_cat_id == '3'), 3] <- 1
predictor_matrix[which(nodes$age_cat_id == '4'), 4] <- 1
predictor_matrix[which(nodes$age_cat_id == '5'), 5] <- 1
predictor_matrix[which(nodes$age_cat_id == '6'), 6] <- 1
predictor_matrix[which(nodes$age_cat_id == '7'), 7] <- 1
predictor_matrix[which(nodes$sex == 'F'), 8] <- 1
predictor_matrix[which(nodes$sex == 'M'), 9] <- 1
predictor_matrix[which(nodes$sex == 'U'), 10] <- 1
predictor_matrix

target <- function(params, index) loglik(params, centrality_matrix_std, predictor_matrix, index)

### run new model ####
chain <- metropolis(target, initial = rep(0,12), iterations=100000, thin=100, refresh=10000)

colnames(chain) <- c("intercept",
                     "beta_age1", "beta_age2", "beta_age3", "beta_age4", "beta_age5", "beta_age6", "beta_age7", 
                     "beta_females", "beta_males", "beta_unks", "sigma")
head(chain)

par(mfrow=c(4,3))
for (i in 1:12) {
  plot(chain[, i], type="l", 
       ylab = colnames(chain)[i])
}

par(mfrow = c(1,1))
plot(density(centrality_matrix_std[1, ]), ylim=c(0, 1.5), main="", xlab="Standardised network strength")
sample_ids <- sample(1:1000, size=200)
for (i in sample_ids) {
  pred <- rnorm(8, mean = chain[i, "intercept"] + 
                  chain[i,"beta_age1"]*predictor_matrix[,1] + 
                  chain[i,"beta_age2"]*predictor_matrix[,2] + 
                  chain[i,"beta_age3"]*predictor_matrix[,3] + 
                  chain[i,"beta_age4"]*predictor_matrix[,4] + 
                  chain[i,"beta_age5"]*predictor_matrix[,5] + 
                  chain[i,"beta_age6"]*predictor_matrix[,6] + 
                  chain[i,"beta_age7"]*predictor_matrix[,7] + 
                  chain[i,"beta_females"]*predictor_matrix[,8] + 
                  chain[i,"beta_males"]*predictor_matrix[,9] + 
                  chain[i,'beta_unks']*predictor_matrix[,10],
                sd = exp(chain[i, "sigma"]))
  lines(density(centrality_matrix_std[i, ]), col=rgb(0, 0, 0, 0.25))
  lines(density(pred), col=rgb(0, 0, 1, 0.25))
}

coefficient_quantiles <- t(apply(chain, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975))))
coefficient_quantiles

beta_difference <- chain[, "beta_females"] - chain[, "beta_males"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

beta_difference <- chain[, "beta_females"] - chain[, "beta_unks"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

beta_difference <- chain[, "beta_males"] - chain[, "beta_unks"]
quantile(beta_difference, probs=c(0.025, 0.5, 0.975))

#
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


