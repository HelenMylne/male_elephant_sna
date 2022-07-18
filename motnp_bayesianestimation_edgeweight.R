## Bayesian analysis of ALERT data
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
R.Version()
rstan::stan_version()

# set stan path
set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')

# set seed
set.seed(12345)

################ 1) Draw DAGS ################
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

### add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$all_events  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$all_events

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

################ 4) Define likelihood distributions, model and parameters to estimate ################
# Binomial model using a beta distribution where shape 1 = times together and shape 2 = times apart

### Code from Dan 2nd February 2022
#data{
#int n_dyads;
#vector[n_dyads] apart;
#vector[n_dyads] together;
#}
#parameters {
#  vector<lower=0,upper=1>[n_dyads] weight; 
#  model {
#    weight ~ beta( 4 + together, 10 +  apart );
#}

################ 5) Run model on simulated data ################
# create data list
simdat_ls <- list(
  n_dyads  = nrow(dyads),        # Number of dyads
  together = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart         # Number of sightings of each dyad apart
)

#### a = 4, b = 10  ####
### Dan model stancode provided 2nd February 2022
# "simpleBetaNet_DWF_22.02.02.stan"
#data{
#int n_dyads;
#vector[n_dyads] apart;
#vector[n_dyads] together;
#}
#parameters {
#  vector<lower=0,upper=1>[n_dyads] weight; 
#  model {
#    weight ~ beta( 4 + together, 10 +  apart );
#}

### Compile Stan model
mod_4.10 <- cmdstan_model("models/simpleBetaNet_DWF_22.02.02.stan")

### Fit model
edge_weight_4.10 <- mod_4.10$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
edge_weight_4.10$summary()
#   variable         mean      median      sd     mad       q5           q95     rhat  ess_bulk  ess_tail
#   <chr>           <dbl>       <dbl>   <dbl>   <dbl>       <dbl>        <dbl>   <dbl>    <dbl>     <dbl>
# 1 lp__         -134935.     -134935    60.5    57.8       -135036    -134835   1.00     1213.     1838.
# 2 weight[1]       0.206       0.200  0.0691  0.0716       0.104        0.329   1.00     7862.     2248.
# 3 weight[2]       0.264       0.260  0.0749  0.0755       0.146        0.397   1.00    11274.     2330.
# 4 weight[3]       0.266       0.261  0.0754  0.0766       0.151        0.395   1.00     8957.     2636.
# 5 weight[4]       0.353       0.350  0.0817  0.0830       0.224        0.495   1.00     8623.     2676.
# 6 weight[5]       0.117       0.110  0.0532  0.0518       0.0438       0.216   1.00     8774.     2404.
# 7 weight[6]       0.177       0.171  0.0650  0.0629       0.0815       0.293   1.00     8541.     2359.
# 8 weight[7]       0.206       0.202  0.0681  0.0684       0.104        0.327   1.00     8776.     2684.
# 9 weight[8]       0.267       0.262  0.0762  0.0783       0.151        0.397   1.00     9888.     2852.
#10 weight[9]       0.118       0.110  0.0555  0.0545       0.0433       0.224   1.00     8904.     2946.
# ... with 7131 more rows
output1 <- read_cmdstan_csv(edge_weight_4.10$output_files()[1])
output2 <- read_cmdstan_csv(edge_weight_4.10$output_files()[2])
output3 <- read_cmdstan_csv(edge_weight_4.10$output_files()[3])
output4 <- read_cmdstan_csv(edge_weight_4.10$output_files()[4])

draws1_4.10 <- as.data.frame(output1$post_warmup_draws)
draws2_4.10 <- as.data.frame(output2$post_warmup_draws)
draws3_4.10 <- as.data.frame(output3$post_warmup_draws)
draws4_4.10 <- as.data.frame(output4$post_warmup_draws)
draws_4.10 <- rbind(draws1_4.10, draws2_4.10, draws3_4.10, draws4_4.10)
rm(output1, output2, output3, output4)

# all chains similar, but not the same
mean(draws1_4.10$`1.weight[1]`) ; mean(draws2_4.10$`1.weight[1]`) ; mean(draws3_4.10$`1.weight[1]`) ; mean(draws4_4.10$`1.weight[1]`)
mean(draws1_4.10$`1.weight[2]`) ; mean(draws2_4.10$`1.weight[2]`) ; mean(draws3_4.10$`1.weight[2]`) ; mean(draws4_4.10$`1.weight[2]`)
mean(draws1_4.10$`1.weight[3]`) ; mean(draws2_4.10$`1.weight[3]`) ; mean(draws3_4.10$`1.weight[3]`) ; mean(draws4_4.10$`1.weight[3]`)

# plot each chain individually to check mixing -- all mixed well
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_4.10$`1.weight[1]`, col = 'purple')
lines(draws2_4.10$`1.weight[1]`, col = 'green')
lines(draws3_4.10$`1.weight[1]`, col = 'blue')
lines(draws4_4.10$`1.weight[1]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_4.10$`1.weight[2]`, col = 'purple')
lines(draws2_4.10$`1.weight[2]`, col = 'green')
lines(draws3_4.10$`1.weight[2]`, col = 'blue')
lines(draws4_4.10$`1.weight[2]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_4.10$`1.weight[3]`, col = 'purple')
lines(draws2_4.10$`1.weight[3]`, col = 'green')
lines(draws3_4.10$`1.weight[3]`, col = 'blue')
lines(draws4_4.10$`1.weight[3]`, col = 'magenta')

# build traceplots -- reasonable, but underestimating family and group relationships
plot(draws1_4.10$`1.weight[1]`, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1_4.10$`1.weight[2]`, col = 'tan')
lines(draws1_4.10$`1.weight[3]`, col = 'orange')
lines(draws1_4.10$`1.weight[4]`, col = 'green')
lines(draws1_4.10$`1.weight[5]`, col = 'chocolate')
lines(draws1_4.10$`1.weight[6]`, col = 'blue')
lines(draws1_4.10$`1.weight[7]`, col = 'red')
lines(draws1_4.10$`1.weight[8]`, col = 'seagreen')
lines(draws1_4.10$`1.weight[9]`, col = 'purple')
lines(draws1_4.10$`1.weight[10]`,col = 'magenta')
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[1]+1], col = 'black')       # traceplot of mother-calf relationships
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[2]+1], col = 'tan')         # +1 because first column is not linked to a dyad
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[3]+1], col = 'orange')
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[4]+1], col = 'green')
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[6]+1], col = 'blue')
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[7]+1], col = 'red')
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[9]+1], col = 'purple')
lines(draws1_4.10[which(dyads$family_1 == dyads$id_2)[10]+1], col = 'magenta')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[1]+1], col = 'black')      # traceplots of elephants in the same group but which are not mother-calf pairs
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[2]+1], col = 'tan')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[3]+1], col = 'orange')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[4]+1], col = 'green')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[6]+1], col = 'blue')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[7]+1], col = 'red')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[9]+1], col = 'purple')
lines(draws1_4.10[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[10]+1], col = 'magenta')

# mean and credible interval for each dyad weight
means_4.10 <- data.frame(dyad = dyads$dyad,               # these are matching up to the right dyads because the traceplots above use the row numbers to specify which ones to plot and they come out at the right approximate edge weight
                         mean = apply(draws_4.10[2:7141], 2, mean), 
                         median = apply(draws_4.10[2:7141], 2, median)) 
draws1_4.10_pi <- apply(draws1_4.10[2:7141], 2, rethinking::PI)
means_4.10$pi_lwr <- draws1_4.10_pi[1,]
means_4.10$pi_upr <- draws1_4.10_pi[2,]

# boxplot comparing types of dyads
plot_data_4.10 <- left_join(x = means_4.10, y = dyads, by = 'dyad')
ggplot(data = plot_data_4.10, aes(x = pair_type, y = mean, fill = pair_type))+
  geom_boxplot()+
  geom_jitter(aes(x = pair_type, y = edge, fill = pair_type),
              pch = 4, colour = col.alpha('red', 0.2))+
  scale_fill_viridis_d()+
  #geom_line(aes(y = mean(edge[pair_type == 'family'])), colour = 'purple', lty = 2)+
  theme_classic()

# clear large objects from environment for this attempt
rm(draws_4.10, draws1_4.10, draws1_4.10_pi, draws2_4.10, draws3_4.10, draws4_4.10, means_4.10, plot_data_4.10)

#### a = b = 5  ####
# Stancode provided 2nd February 2022, adapted to reduce prior strength as currently prevents even the mother-calf pairings exceeding mean of 0.7

### Compile Stan model
mod_5.5 <- cmdstan_model("models/simpleBetaNet_HKM_5.5_22.02.03.stan")

### Fit model
edge_weight_5.5 <- mod_5.5$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
edge_weight_5.5$summary()
#   variable         mean      median      sd     mad          q5          q95    rhat ess_bulk ess_tail
#   <chr>           <dbl>       <dbl>   <dbl>   <dbl>       <dbl>        <dbl>   <dbl>    <dbl>    <dbl>
# 1 lp__         -131826.     -131825    61.0    60.8     -131928      -131728    1.01    1164.    1939.
# 2 weight[1]       0.266       0.261  0.0806  0.0786       0.144        0.405    1.00   10055.    2474.
# 3 weight[2]       0.332       0.328  0.0847  0.0867       0.198        0.480    1.00   10337.    2763.
# 4 weight[3]       0.334       0.330  0.0835  0.0862       0.202        0.477    1.00    9554.    2314.
# 5 weight[4]       0.435       0.434  0.0889  0.0908       0.292        0.583    1.00    9951.    2456.
# 6 weight[5]       0.166       0.160  0.0655  0.0628      0.0706        0.285    1.00    9182.    2303.
# 7 weight[6]       0.232       0.227  0.0745  0.0725       0.118        0.368    1.00   10186.    2420.
# 8 weight[7]       0.267       0.261  0.0791  0.0813       0.145        0.405    1.00    8535.    2602.
# 9 weight[8]       0.333       0.330  0.0834  0.0880       0.204        0.474    1.00   11038.    2936.
#10 weight[9]       0.165       0.159  0.0663  0.0668      0.0717        0.285    1.00    9208.    2743.
# … with 7,131 more rows
output1 <- read_cmdstan_csv(edge_weight_5.5$output_files()[1])
output2 <- read_cmdstan_csv(edge_weight_5.5$output_files()[2])
output3 <- read_cmdstan_csv(edge_weight_5.5$output_files()[3])
output4 <- read_cmdstan_csv(edge_weight_5.5$output_files()[4])

draws1_5.5 <- as.data.frame(output1$post_warmup_draws)
draws2_5.5 <- as.data.frame(output2$post_warmup_draws)
draws3_5.5 <- as.data.frame(output3$post_warmup_draws)
draws4_5.5 <- as.data.frame(output4$post_warmup_draws)
draws_5.5 <- rbind(draws1_5.5, draws2_5.5, draws3_5.5, draws4_5.5)
rm(output1, output2, output3, output4)

# check chain mixing -- well mixed
mean(draws1_5.5$`1.weight[1]`) ; mean(draws2_5.5$`1.weight[1]`) ; mean(draws3_5.5$`1.weight[1]`) ; mean(draws4_5.5$`1.weight[1]`)
mean(draws1_5.5$`1.weight[2]`) ; mean(draws2_5.5$`1.weight[2]`) ; mean(draws3_5.5$`1.weight[2]`) ; mean(draws4_5.5$`1.weight[2]`)
mean(draws1_5.5$`1.weight[3]`) ; mean(draws2_5.5$`1.weight[3]`) ; mean(draws3_5.5$`1.weight[3]`) ; mean(draws4_5.5$`1.weight[3]`)

# plot each chain individually to check mixing -- all mixed well
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_5.5$`1.weight[1]`, col = 'purple')
lines(draws2_5.5$`1.weight[1]`, col = 'green')
lines(draws3_5.5$`1.weight[1]`, col = 'blue')
lines(draws4_5.5$`1.weight[1]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_5.5$`1.weight[2]`, col = 'purple')
lines(draws2_5.5$`1.weight[2]`, col = 'green')
lines(draws3_5.5$`1.weight[2]`, col = 'blue')
lines(draws4_5.5$`1.weight[2]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_5.5$`1.weight[3]`, col = 'purple')
lines(draws2_5.5$`1.weight[3]`, col = 'green')
lines(draws3_5.5$`1.weight[3]`, col = 'blue')
lines(draws4_5.5$`1.weight[3]`, col = 'magenta')

# build traceplots -- reasonable, but underestimating family and group relationships. almost now too strong in both directions, but real data has more points so might be sufficient
plot(draws1_5.5$`1.weight[1]`, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1_5.5$`1.weight[2]`, col = 'tan')
lines(draws1_5.5$`1.weight[3]`, col = 'orange')
lines(draws1_5.5$`1.weight[4]`, col = 'green')
lines(draws1_5.5$`1.weight[5]`, col = 'chocolate')
lines(draws1_5.5$`1.weight[6]`, col = 'blue')
lines(draws1_5.5$`1.weight[7]`, col = 'red')
lines(draws1_5.5$`1.weight[8]`, col = 'seagreen')
lines(draws1_5.5$`1.weight[9]`, col = 'purple')
lines(draws1_5.5$`1.weight[10]`,col = 'magenta')
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[1]+1], col = 'black')       # traceplot of mother-calf relationships
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[2]+1], col = 'tan')         # +1 because first column is not linked to a dyad
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[3]+1], col = 'orange')
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[4]+1], col = 'green')
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[6]+1], col = 'blue')
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[7]+1], col = 'red')
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[9]+1], col = 'purple')
lines(draws1_5.5[which(dyads$family_1 == dyads$id_2)[10]+1], col = 'magenta')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[1]+1], col = 'black')      # traceplots of elephants in the same group but which are not mother-calf pairs
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[2]+1], col = 'tan')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[3]+1], col = 'orange')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[4]+1], col = 'green')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[6]+1], col = 'blue')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[7]+1], col = 'red')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[9]+1], col = 'purple')
lines(draws1_5.5[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[10]+1], col = 'magenta')

# mean and credible interval for each dyad weight
means_5.5 <- data.frame(dyad = dyads$dyad,
                        mean = apply(draws_5.5[2:7141], 2, mean),
                        median = apply(draws_5.5[2:7141], 2, median)) 
draws1_5.5_pi <- apply(draws1_5.5[2:7141], 2, rethinking::PI)
means_5.5$pi_lwr <- draws1_5.5_pi[1,]
means_5.5$pi_upr <- draws1_5.5_pi[2,]

# boxplot comparing types of dyads
plot_data_5.5 <- left_join(x = means_5.5, y = dyads, by = 'dyad')
ggplot(data = plot_data_5.5, aes(x = pair_type, y = mean, fill = pair_type))+
  geom_boxplot()+
  geom_jitter(aes(x = pair_type, y = edge, fill = pair_type),
              pch = 4, colour = col.alpha('red', 0.2))+
  scale_fill_viridis_d()+
  #geom_line(aes(y = mean(edge[pair_type == 'family'])), colour = 'purple', lty = 2)+
  theme_classic()

# clear large objects from environment for this attempt
rm(draws_5.5, draws1_5.5, draws1_5.5_pi, draws2_5.5, draws3_5.5, draws4_5.5, means_5.5, plot_data_5.5)

#### a = b = 2  ####
# As above, but now priors added are just values 2 and 2 -- pretty flat prior, but still more regularising than 1,1. Just enough to help push it away from extremes, but not enough that masses of data are required to drag the distributions toward the true values.
# Having a weak prior might be fine here, given that it's not a hypothesis-connected parameter. However, what's important to look at here is what happens for well-sampled individuals, and what happens for poorly sampled ones. We don't want it to be too enthusiastic about the poorly-sampled individuals.

### Compile Stan model
mod_2.2 <- cmdstan_model("models/simpleBetaNet_HKM_2.2_22.02.03.stan")

### Fit model
edge_weight_2.2 <- mod_2.2$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
edge_weight_2.2$summary()
#   variable        mean      median       sd     mad          q5        q95    rhat ess_bulk ess_tail
#   <chr>          <dbl>       <dbl>    <dbl>   <dbl>       <dbl>      <dbl>   <dbl>    <dbl>    <dbl>
# 1 lp__         -96302.     -96301.     60.4    60.8     -96401.     -96201.   1.00    1097.    1671.
# 2 weight[1]      0.209       0.202   0.0811  0.0814      0.0886      0.358    1.00    9184.    2203.
# 3 weight[2]      0.292       0.287   0.0917  0.0952       0.150      0.450    1.00   11174.    2581.
# 4 weight[3]      0.292       0.287   0.0917  0.0941       0.155      0.452    1.00   10926.    2534.
# 5 weight[4]      0.418       0.416   0.0984   0.101       0.258      0.580    1.00   10820.    2504.
# 6 weight[5]     0.0835      0.0715   0.0565  0.0511      0.0152      0.194    1.00    8490.    1832.
# 7 weight[6]      0.166       0.155   0.0780  0.0770      0.0585      0.311    1.00    9021.    1997.
# 8 weight[7]      0.207       0.199   0.0794  0.0797      0.0915      0.348    1.00    9618.    2459.
# 9 weight[8]      0.291       0.284   0.0929  0.0951       0.150      0.453    1.00    9875.    2330.
#10 weight[9]     0.0832      0.0725   0.0548  0.0484      0.0159      0.190    1.00    7881.    2618.
# … with 7,131 more rows
output1 <- read_cmdstan_csv(edge_weight_2.2$output_files()[1])
output2 <- read_cmdstan_csv(edge_weight_2.2$output_files()[2])
output3 <- read_cmdstan_csv(edge_weight_2.2$output_files()[3])
output4 <- read_cmdstan_csv(edge_weight_2.2$output_files()[4])

draws1_2.2 <- as.data.frame(output1$post_warmup_draws)
draws2_2.2 <- as.data.frame(output2$post_warmup_draws)
draws3_2.2 <- as.data.frame(output3$post_warmup_draws)
draws4_2.2 <- as.data.frame(output4$post_warmup_draws)
draws_2.2 <- rbind(draws1_2.2, draws2_2.2, draws3_2.2, draws4_2.2)
rm(output1, output2, output3, output4)

# check chain mixing -- well mixed
mean(draws1_2.2$`1.weight[1]`) ; mean(draws2_2.2$`1.weight[1]`) ; mean(draws3_2.2$`1.weight[1]`) ; mean(draws4_2.2$`1.weight[1]`)
mean(draws1_2.2$`1.weight[2]`) ; mean(draws2_2.2$`1.weight[2]`) ; mean(draws3_2.2$`1.weight[2]`) ; mean(draws4_2.2$`1.weight[2]`)
mean(draws1_2.2$`1.weight[3]`) ; mean(draws2_2.2$`1.weight[3]`) ; mean(draws3_2.2$`1.weight[3]`) ; mean(draws4_2.2$`1.weight[3]`)

# plot each chain individually to check mixing -- all mixed well
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_2.2$`1.weight[1]`, col = 'purple')
lines(draws2_2.2$`1.weight[1]`, col = 'green')
lines(draws3_2.2$`1.weight[1]`, col = 'blue')
lines(draws4_2.2$`1.weight[1]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_2.2$`1.weight[2]`, col = 'purple')
lines(draws2_2.2$`1.weight[2]`, col = 'green')
lines(draws3_2.2$`1.weight[2]`, col = 'blue')
lines(draws4_2.2$`1.weight[2]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_2.2$`1.weight[3]`, col = 'purple')
lines(draws2_2.2$`1.weight[3]`, col = 'green')
lines(draws3_2.2$`1.weight[3]`, col = 'blue')
lines(draws4_2.2$`1.weight[3]`, col = 'magenta')

# build traceplots -- looks good, but some clique relationships substantially higher than family
plot(draws1_2.2$`1.weight[1]`, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1_2.2$`1.weight[2]`, col = 'tan')
lines(draws1_2.2$`1.weight[3]`, col = 'orange')
lines(draws1_2.2$`1.weight[4]`, col = 'green')
lines(draws1_2.2$`1.weight[5]`, col = 'chocolate')
lines(draws1_2.2$`1.weight[6]`, col = 'blue')
lines(draws1_2.2$`1.weight[7]`, col = 'red')
lines(draws1_2.2$`1.weight[8]`, col = 'seagreen')
lines(draws1_2.2$`1.weight[9]`, col = 'purple')
lines(draws1_2.2$`1.weight[10]`,col = 'magenta')
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[1]+1], col = 'black')       # traceplot of mother-calf relationships
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[2]+1], col = 'tan')         # +1 because first column is not linked to a dyad
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[3]+1], col = 'orange')
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[4]+1], col = 'green')
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[6]+1], col = 'blue')
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[7]+1], col = 'red')
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[9]+1], col = 'purple')
lines(draws1_2.2[which(dyads$family_1 == dyads$id_2)[10]+1], col = 'magenta')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[1]+1], col = 'black')      # traceplots of elephants in the same group but which are not mother-calf pairs
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[2]+1], col = 'tan')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[3]+1], col = 'orange')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[4]+1], col = 'green')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[6]+1], col = 'blue')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[7]+1], col = 'red')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[9]+1], col = 'purple')
lines(draws1_2.2[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[10]+1], col = 'magenta')

# mean and credible interval for each dyad weight
means_2.2 <- data.frame(dyad = dyads$dyad,
                        mean = apply(draws_2.2[2:7141], 2, mean), 
                        median = apply(draws_2.2[2:7141], 2, median))
draws1_2.2_pi <- apply(draws1_2.2[2:7141], 2, rethinking::PI)
means_2.2$pi_lwr <- draws1_2.2_pi[1,]
means_2.2$pi_upr <- draws1_2.2_pi[2,]

# boxplot comparing types of dyads
plot_data_2.2 <- left_join(x = means_2.2, y = dyads, by = 'dyad')
ggplot(data = plot_data_2.2, aes(x = pair_type, y = mean, fill = pair_type))+
  geom_boxplot()+
  geom_jitter(aes(x = pair_type, y = edge, fill = pair_type),
              pch = 4, colour = col.alpha('red', 0.2))+
  scale_fill_viridis_d()+
  #geom_line(aes(y = mean(edge[pair_type == 'family'])), colour = 'purple', lty = 2)+
  theme_classic()

# clear large objects from environment for this attempt
rm(draws_2.2, draws1_2.2, draws1_2.2_pi, draws2_2.2, draws3_2.2, draws4_2.2, means_2.2, plot_data_2.2)

#### pooling by relationship type -- WRITE MODEL ####
# The problem is that the distribution of the network weights has two peaks. One way around this is to model them based on the type of link. That way you could have pooling prior for each type of dyad. 

### Compile Stan model
# mod_pool <- cmdstan_model("models/simpleBetaNet_HKM_pool_22.02.07.stan") -- need to look in Rethinking to see how this will work

### Fit model
edge_weight_pool <- mod_pool$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
edge_weight_pool$summary()
#   variable        mean      median       sd     mad          q5        q95    rhat ess_bulk ess_tail
#   <chr>          <dbl>       <dbl>    <dbl>   <dbl>       <dbl>      <dbl>   <dbl>    <dbl>    <dbl>
# 1 lp__         -96302.     -96301.     60.4    60.8     -96401.     -96201.   1.00    1097.    1671.
# 2 weight[1]      0.209       0.202   0.0811  0.0814      0.0886      0.358    1.00    9184.    2203.
# 3 weight[2]      0.292       0.287   0.0917  0.0952       0.150      0.450    1.00   11174.    2581.
# 4 weight[3]      0.292       0.287   0.0917  0.0941       0.155      0.452    1.00   10926.    2534.
# 5 weight[4]      0.418       0.416   0.0984   0.101       0.258      0.580    1.00   10820.    2504.
# 6 weight[5]     0.0835      0.0715   0.0565  0.0511      0.0152      0.194    1.00    8490.    1832.
# 7 weight[6]      0.166       0.155   0.0780  0.0770      0.0585      0.311    1.00    9021.    1997.
# 8 weight[7]      0.207       0.199   0.0794  0.0797      0.0915      0.348    1.00    9618.    2459.
# 9 weight[8]      0.291       0.284   0.0929  0.0951       0.150      0.453    1.00    9875.    2330.
#10 weight[9]     0.0832      0.0725   0.0548  0.0484      0.0159      0.190    1.00    7881.    2618.
# … with 7,131 more rows
output1 <- read_cmdstan_csv(edge_weight_pool$output_files()[1])
output2 <- read_cmdstan_csv(edge_weight_pool$output_files()[2])
output3 <- read_cmdstan_csv(edge_weight_pool$output_files()[3])
output4 <- read_cmdstan_csv(edge_weight_pool$output_files()[4])

draws1_pool <- as.data.frame(output1$post_warmup_draws)
draws2_pool <- as.data.frame(output2$post_warmup_draws)
draws3_pool <- as.data.frame(output3$post_warmup_draws)
draws4_pool <- as.data.frame(output4$post_warmup_draws)
draws_pool <- rbind(draws1_pool, draws2_pool, draws3_pool, draws4_pool)
rm(output1, output2, output3, output4)

# check chain mixing -- well mixed
mean(draws1_pool$`1.weight[1]`) ; mean(draws2_pool$`1.weight[1]`) ; mean(draws3_pool$`1.weight[1]`) ; mean(draws4_pool$`1.weight[1]`)
mean(draws1_pool$`1.weight[2]`) ; mean(draws2_pool$`1.weight[2]`) ; mean(draws3_pool$`1.weight[2]`) ; mean(draws4_pool$`1.weight[2]`)
mean(draws1_pool$`1.weight[3]`) ; mean(draws2_pool$`1.weight[3]`) ; mean(draws3_pool$`1.weight[3]`) ; mean(draws4_pool$`1.weight[3]`)

# plot each chain individually to check mixing -- all mixed well
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_pool$`1.weight[1]`, col = 'purple')
lines(draws2_pool$`1.weight[1]`, col = 'green')
lines(draws3_pool$`1.weight[1]`, col = 'blue')
lines(draws4_pool$`1.weight[1]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_pool$`1.weight[2]`, col = 'purple')
lines(draws2_pool$`1.weight[2]`, col = 'green')
lines(draws3_pool$`1.weight[2]`, col = 'blue')
lines(draws4_pool$`1.weight[2]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_pool$`1.weight[3]`, col = 'purple')
lines(draws2_pool$`1.weight[3]`, col = 'green')
lines(draws3_pool$`1.weight[3]`, col = 'blue')
lines(draws4_pool$`1.weight[3]`, col = 'magenta')

# build traceplots -- looks good, but some clique relationships substantially higher than family
plot(draws1_pool$`1.weight[1]`, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1_pool$`1.weight[2]`, col = 'tan')
lines(draws1_pool$`1.weight[3]`, col = 'orange')
lines(draws1_pool$`1.weight[4]`, col = 'green')
lines(draws1_pool$`1.weight[5]`, col = 'chocolate')
lines(draws1_pool$`1.weight[6]`, col = 'blue')
lines(draws1_pool$`1.weight[7]`, col = 'red')
lines(draws1_pool$`1.weight[8]`, col = 'seagreen')
lines(draws1_pool$`1.weight[9]`, col = 'purple')
lines(draws1_pool$`1.weight[10]`,col = 'magenta')
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[1]+1], col = 'black')       # traceplot of mother-calf relationships
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[2]+1], col = 'tan')         # +1 because first column is not linked to a dyad
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[3]+1], col = 'orange')
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[4]+1], col = 'green')
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[6]+1], col = 'blue')
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[7]+1], col = 'red')
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[9]+1], col = 'purple')
lines(draws1_pool[which(dyads$family_1 == dyads$id_2)[10]+1], col = 'magenta')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[1]+1], col = 'black')      # traceplots of elephants in the same group but which are not mother-calf pairs
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[2]+1], col = 'tan')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[3]+1], col = 'orange')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[4]+1], col = 'green')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[6]+1], col = 'blue')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[7]+1], col = 'red')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[9]+1], col = 'purple')
lines(draws1_pool[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[10]+1], col = 'magenta')

# mean and credible interval for each dyad weight
means_pool <- data.frame(dyad = dyads$dyad,
                        mean = apply(draws_pool[2:7141], 2, mean), 
                        median = apply(draws_pool[2:7141], 2, median))
draws1_pool_pi <- apply(draws1_pool[2:7141], 2, rethinking::PI)
means_pool$pi_lwr <- draws1_pool_pi[1,]
means_pool$pi_upr <- draws1_pool_pi[2,]

# boxplot comparing types of dyads
plot_data_pool <- left_join(x = means_pool, y = dyads, by = 'dyad')
ggplot(data = plot_data_pool, aes(x = pair_type, y = mean, fill = pair_type))+
  geom_boxplot()+
  geom_jitter(aes(x = pair_type, y = edge, fill = pair_type),
              pch = 4, colour = col.alpha('red', 0.2))+
  scale_fill_viridis_d()+
  #geom_line(aes(y = mean(edge[pair_type == 'family'])), colour = 'purple', lty = 2)+
  theme_classic()




################ 6) Run model on real standardised data -- Binomial model to calculate SRI with uncertainty ################
#### a = b = 2 ####
### create data list -- can contain no NA values in any column, even if column is not specified in model
counts_ls <- list(
  n_dyads  = nrow(counts_df),          # total number of times one or other of the dyad was observed
  together = counts_df$all_events,     # count number of sightings seen together
  apart    = counts_df$apart)          # count number of sightings seen apart

### run model
mod_2.2

### Fit model (slow)
weight_motnp_2.2 <- mod_2.2$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_motnp_2.2
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
weight_motnp_2.2$summary()  # exceeds vector memory on my laptop
output1 <- read_cmdstan_csv(weight_motnp_2.2$output_files()[1])
output2 <- read_cmdstan_csv(weight_motnp_2.2$output_files()[2])
output3 <- read_cmdstan_csv(weight_motnp_2.2$output_files()[3])
output4 <- read_cmdstan_csv(weight_motnp_2.2$output_files()[4])

draws1_motnp2.2 <- as.data.frame(output1$post_warmup_draws)
draws2_motnp2.2 <- as.data.frame(output2$post_warmup_draws)
draws3_motnp2.2 <- as.data.frame(output3$post_warmup_draws)
draws4_motnp2.2 <- as.data.frame(output4$post_warmup_draws)
draws_motnp2.2 <- rbind(draws1_motnp2.2, draws2_motnp2.2, draws3_motnp2.2, draws4_motnp2.2)
rm(output1, output2, output3, output4)

colnames(draws_motnp2.2)[2:106954] <- counts_df$dyad

# Assign random set of columns to check
plot_cols <- sample(x = 2:106954, size = 30, replace = F)

# tidy data -- vector memory exhausted
tidy_draws_2.2 <- pivot_longer(draws_motnp2.2[,2:106954], cols = everything(), names_to = 'dyad', values_to = 'draw')
tidy_draws_2.2$chain <- rep(1:4, each = 106953000)
tidy_draws_2.2$index <- rep(rep(1:1000, each = 106953),4)
head(tidy_draws_2.2, 10)
tail(tidy_draws_2.2, 10)

# subset for Sierra (F52) -- family = U17, F60, U21 and F98
tidy_sierra <- tidy_draws[tidy_draws$dyad == 'F52_M40' | tidy_draws$dyad == 'F52_M26' | tidy_draws$dyad == 'F52_F8' | 
                            tidy_draws$dyad == 'F52_M15' | tidy_draws$dyad == 'F52_F98' | tidy_draws$dyad == 'F52_U17' | 
                            tidy_draws$dyad == 'F52_U21' | tidy_draws$dyad == 'F52_F60',]
tidy_sierra <- separate(tidy_sierra, dyad, sep = '_', remove = F, into = c('F52','id'))
tidy_sierra$ID <- ifelse(tidy_sierra$id == 'U17', 'U17 (calf)',
                         ifelse(tidy_sierra$id == 'U21', 'U21 (relative)',
                                ifelse(tidy_sierra$id == 'F60', 'F60 (relative)',
                                       ifelse(tidy_sierra$id == 'F98', 'F98 (relative)', 
                                              ifelse(tidy_sierra$id == 'F8', 'F8 (unrelated female)',
                                                     ifelse(tidy_sierra$id == 'M26', 'M26 (adult male)',
                                                            ifelse(tidy_sierra$id == 'M15', 'M15 (unrelated pubescent male)', 'M40 (adult male)')))))))
tidy_sierra$ID <- factor(tidy_sierra$ID, levels = c('U17 (calf)', 'F60 (relative)', 'U21 (relative)', 'F98 (relative)', 'M40 (adult male)','M26 (adult male)', 'F8 (unrelated female)','M15 (unrelated pubescent male)'))

#### pooling by relationship type -- WRITE MODEL ####
################ 7) Summarise and plot edge weights --  a = b = 2 ################
### build traceplots ####
# random dyads -- many are very wide, high uncertainty
plot(draws_motnp2.2[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_motnp2.2[,plot_cols[2]], col = 'tan')
lines(draws_motnp2.2[,plot_cols[3]], col = 'orange')
lines(draws_motnp2.2[,plot_cols[4]], col = 'green')
lines(draws_motnp2.2[,plot_cols[5]], col = 'chocolate')
lines(draws_motnp2.2[,plot_cols[6]], col = 'blue')
lines(draws_motnp2.2[,plot_cols[7]], col = 'red')
lines(draws_motnp2.2[,plot_cols[8]], col = 'seagreen')
lines(draws_motnp2.2[,plot_cols[9]], col = 'purple')
lines(draws_motnp2.2[,plot_cols[10]],col = 'magenta')
lines(draws_motnp2.2[,plot_cols[11]],col = 'black')
lines(draws_motnp2.2[,plot_cols[12]], col = 'tan')
lines(draws_motnp2.2[,plot_cols[13]], col = 'orange')
lines(draws_motnp2.2[,plot_cols[14]], col = 'green')
lines(draws_motnp2.2[,plot_cols[15]], col = 'chocolate')
lines(draws_motnp2.2[,plot_cols[16]], col = 'blue')
lines(draws_motnp2.2[,plot_cols[17]], col = 'red')
lines(draws_motnp2.2[,plot_cols[18]], col = 'seagreen')
lines(draws_motnp2.2[,plot_cols[19]], col = 'purple')
lines(draws_motnp2.2[,plot_cols[20]],col = 'magenta')
lines(draws_motnp2.2[,plot_cols[21]],col = 'black')
lines(draws_motnp2.2[,plot_cols[22]], col = 'tan')
lines(draws_motnp2.2[,plot_cols[23]], col = 'orange')
lines(draws_motnp2.2[,plot_cols[24]], col = 'green')
lines(draws_motnp2.2[,plot_cols[25]], col = 'chocolate')
lines(draws_motnp2.2[,plot_cols[26]], col = 'blue')
lines(draws_motnp2.2[,plot_cols[27]], col = 'red')
lines(draws_motnp2.2[,plot_cols[28]], col = 'seagreen')
lines(draws_motnp2.2[,plot_cols[29]], col = 'purple')
lines(draws_motnp2.2[,plot_cols[30]],col = 'magenta')

# Sierra = F52, herd members = F60+U21+F98, calf = U17. This looks really good, and not too wide (all seen a good number of times).
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,4000))
lines(draws_motnp2.2$F52_M40, col = 'black')      # non-herd member, adult male
lines(draws_motnp2.2$F52_M15, col = 'tan')        # non-herd member, pubescent male
lines(draws_motnp2.2$F52_M203,col = 'orange')     # non-herd member, adult male
lines(draws_motnp2.2$F52_M26, col = 'green')      # non-herd member, calf
lines(draws_motnp2.2$F52_F8,  col = 'chocolate')  # non-herd member, adult female
lines(draws_motnp2.2$F52_U9,  col = 'blue')       # non-herd member, calf
lines(draws_motnp2.2$F52_F98, col = 'red')        # herd member most frequently absent from sightings
lines(draws_motnp2.2$F52_U17, col = 'purple')     # calf
lines(draws_motnp2.2$F52_U21, col = 'seagreen')   # sister
lines(draws_motnp2.2$F52_F60, col = 'magenta')    # sister's calf

# plot for Sierra (F52)
(traceplot_sierra <- ggplot(tidy_sierra[tidy_sierra$chain == 1,], aes(x = index, y = draw, colour = ID))+
    geom_line()+
    scale_color_viridis_d()+
    theme_classic()+
    scale_x_continuous('MCMC chain position',  expand = c(0,0))+
    scale_y_continuous('edge weight estimate', expand = c(0.02,0)))
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_sierra_traceplot.pdf',
       plot = traceplot_sierra,
       width = 30, height = 24, units = 'cm')

(traceplot_sierra2 <- ggplot(tidy_sierra, aes(x = index, y = draw, colour = as.factor(chain)))+
    geom_line()+
    scale_color_viridis_d()+
    theme_light()+
    facet_wrap(.~ID)+
    scale_x_continuous('MCMC chain position',  expand = c(0,0))+
    scale_y_continuous('edge weight estimate', expand = c(0.02,0))+
    theme(legend.position = c(0.85, 0.15), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))+
    labs(colour = 'Chain ID'))
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_sierra_traceplot.pdf',
       plot = traceplot_sierra,
       width = 30, height = 24, units = 'cm')

# check chain mixing
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_motnp2.2$`1.weight[1]`, col = 'purple')
lines(draws2_motnp2.2$`1.weight[1]`, col = 'green')
lines(draws3_motnp2.2$`1.weight[1]`, col = 'blue')
lines(draws4_motnp2.2$`1.weight[1]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_motnp2.2$`1.weight[2]`, col = 'purple')
lines(draws2_motnp2.2$`1.weight[2]`, col = 'green')
lines(draws3_motnp2.2$`1.weight[2]`, col = 'blue')
lines(draws4_motnp2.2$`1.weight[2]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_motnp2.2$`1.weight[3]`, col = 'purple')
lines(draws2_motnp2.2$`1.weight[3]`, col = 'green')
lines(draws3_motnp2.2$`1.weight[3]`, col = 'blue')
lines(draws4_motnp2.2$`1.weight[3]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
which(counts_df$dyad == 'F52_U17') # 41192
lines(draws1_motnp2.2$`1.weight[41192]`, col = 'purple')
lines(draws2_motnp2.2$`1.weight[41192]`, col = 'green')
lines(draws3_motnp2.2$`1.weight[41192]`, col = 'blue')
lines(draws4_motnp2.2$`1.weight[41192]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
which(counts_df$dyad == 'F52_F60') # 40896
lines(draws1_motnp2.2$`1.weight[40896]`, col = 'purple')
lines(draws2_motnp2.2$`1.weight[40896]`, col = 'green')
lines(draws3_motnp2.2$`1.weight[40896]`, col = 'blue')
lines(draws4_motnp2.2$`1.weight[40896]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
which(counts_df$dyad == 'F52_U21') # 41197
lines(draws1_motnp2.2$`1.weight[41197]`, col = 'purple')
lines(draws2_motnp2.2$`1.weight[41197]`, col = 'green')
lines(draws3_motnp2.2$`1.weight[41197]`, col = 'blue')
lines(draws4_motnp2.2$`1.weight[41197]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
which(counts_df$dyad == 'F52_F98') # 40937
lines(draws1_motnp2.2$`1.weight[40937]`, col = 'purple')
lines(draws2_motnp2.2$`1.weight[40937]`, col = 'green')
lines(draws3_motnp2.2$`1.weight[40937]`, col = 'blue')
lines(draws4_motnp2.2$`1.weight[40937]`, col = 'magenta')

### density plots ####
# random columns
dens(draws_motnp2.2[,plot_cols[1]], ylim = c(0,50), xlim = c(0,1), las = 1)
dens(add = T, draws_motnp2.2[,plot_cols[2]], col = 'tan')
dens(add = T, draws_motnp2.2[,plot_cols[3]], col = 'orange')
dens(add = T, draws_motnp2.2[,plot_cols[4]], col = 'green')
dens(add = T, draws_motnp2.2[,plot_cols[5]], col = 'chocolate')
dens(add = T, draws_motnp2.2[,plot_cols[6]], col = 'blue')
dens(add = T, draws_motnp2.2[,plot_cols[7]], col = 'red')
dens(add = T, draws_motnp2.2[,plot_cols[8]], col = 'seagreen')
dens(add = T, draws_motnp2.2[,plot_cols[9]], col = 'purple')
dens(add = T, draws_motnp2.2[,plot_cols[10]],col = 'magenta')
dens(add = T, draws_motnp2.2[,plot_cols[11]],col = 'black')
dens(add = T, draws_motnp2.2[,plot_cols[12]],col = 'tan')
dens(add = T, draws_motnp2.2[,plot_cols[13]],col = 'orange')
dens(add = T, draws_motnp2.2[,plot_cols[14]],col = 'green')
dens(add = T, draws_motnp2.2[,plot_cols[15]],col = 'chocolate')
dens(add = T, draws_motnp2.2[,plot_cols[16]],col = 'blue')
dens(add = T, draws_motnp2.2[,plot_cols[17]],col = 'red')
dens(add = T, draws_motnp2.2[,plot_cols[18]],col = 'seagreen')
dens(add = T, draws_motnp2.2[,plot_cols[19]],col = 'purple')
dens(add = T, draws_motnp2.2[,plot_cols[20]],col = 'magenta')
dens(add = T, draws_motnp2.2[,plot_cols[21]],col = 'black')
dens(add = T, draws_motnp2.2[,plot_cols[22]],col = 'tan')
dens(add = T, draws_motnp2.2[,plot_cols[23]],col = 'orange')
dens(add = T, draws_motnp2.2[,plot_cols[24]],col = 'green')
dens(add = T, draws_motnp2.2[,plot_cols[25]],col = 'chocolate')
dens(add = T, draws_motnp2.2[,plot_cols[26]],col = 'blue')
dens(add = T, draws_motnp2.2[,plot_cols[27]],col = 'red')
dens(add = T, draws_motnp2.2[,plot_cols[28]],col = 'seagreen')
dens(add = T, draws_motnp2.2[,plot_cols[29]],col = 'purple')
dens(add = T, draws_motnp2.2[,plot_cols[30]],col = 'magenta')
dens(add = T, draws_motnp2.2[,which(counts_df$dyad == 'F52_U17')+1],col = 'magenta')
dens(add = T, draws_motnp2.2[,which(counts_df$dyad == 'F52_F60')+1],col = 'purple')
dens(add = T, draws_motnp2.2[,which(counts_df$dyad == 'F52_U21')+1],col = 'blue')
dens(add = T, draws_motnp2.2[,which(counts_df$dyad == 'F52_F98')+1],col = 'red')

# density plots
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
means_motnp2.2 <- data.frame(dyad = counts_df$dyad,
                             mean = apply(draws_motnp2.2[,2:106954], 2, mean),
                             median = apply(draws_motnp2.2[,2:106954], 2, median))
draws_pi_motnp2.2 <- apply(draws_motnp2.2[2:106954], 2, rethinking::PI)
means_motnp2.2$pi_lwr <- draws_pi_motnp2.2[1,] ; means_motnp2.2$pi_upr <- draws_pi_motnp2.2[2,]

# organise dem_class
plot_data_motnp2.2 <- left_join(x = means_motnp2.2, y = counts_df, by = 'dyad')
head(plot_data_motnp2.2)
plot_data_motnp2.2$dem_class_1_cat <- as.integer(as.factor(plot_data_motnp2.2$dem_class_1))
plot_data_motnp2.2$dem_class_2_cat <- as.integer(as.factor(plot_data_motnp2.2$dem_class_2))
plot_data_motnp2.2$dem_type_short <- ifelse(plot_data_motnp2.2$dem_class_1_cat <= plot_data_motnp2.2$dem_class_2_cat,
                                   paste(plot_data_motnp2.2$dem_class_1, plot_data_motnp2.2$dem_class_2, sep = '_'),
                                   paste(plot_data_motnp2.2$dem_class_2, plot_data_motnp2.2$dem_class_1, sep = '_'))
sort(unique(plot_data_motnp2.2$dem_type_short))
types <- data.frame(all = sort(unique(plot_data_motnp2.2$dem_type_short)), type1 = NA, type2 = NA)
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
plot_data_motnp2.2 <- left_join(plot_data_motnp2.2, types, by = 'dem_type_short')
plot_data_motnp2.2$join <- paste(plot_data_motnp2.2$join, ' ', sep = ' ')
which(is.na(plot_data_motnp2.2$join))

# boxplot edge weights by demographic type of dyad -- all types
colours <- c('magenta','purple','grey','grey','magenta','purple','grey','blue','purple','purple','purple','blue','purple','grey','grey','grey','grey','grey','magenta','purple','grey','purple','purple','blue','purple','grey','grey')
(edge_vs_demtype_all <- 
    ggplot(data = plot_data_motnp2.2, aes(y = mean, x = join))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 fill = colours)+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown',
         y = 'mean edge weight')+
    coord_flip())
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_all_edge_vs_demtype.pdf',
       plot = edge_vs_demtype_all,
       width = 20, height = 24, units = 'cm')  
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_all_edge_vs_demtype.pdf',
       plot = edge_vs_demtype_all,
       width = 20, height = 24, units = 'cm')  

# boxplot edge weights by demographic type of dyad -- only adults/pubescents
adults_motnp2.2 <- plot_data_motnp2.2[plot_data_motnp2.2$age_cat_id_1 > 2 & plot_data_motnp2.2$age_cat_id_2 > 2,]
adults_motnp2.2 <- adults_motnp2.2[!is.na(adults_motnp2.2$dyad),]
colours <- c('magenta','purple','magenta','purple','grey','blue','purple','blue','purple','magenta','purple','grey','blue','purple')
(edge_vs_demtype <- 
    ggplot(data = adults_motnp2.2, aes(y = mean, x = join))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 fill = colours)+
    theme_classic()+
    labs(x = 'dyad demography (A/P = adult/pubescent, M/F/U = male/female/unknown',
         y = 'mean edge weight')+
    coord_flip())
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_adults_edge_vs_demtype.pdf',
       plot = edge_vs_demtype,
       width = 20, height = 24, units = 'cm')  
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_adults_edge_vs_demtype.pdf',
       plot = edge_vs_demtype,
       width = 20, height = 24, units = 'cm')

# scatter plot of all adult edges by age difference
adults_motnp2.2$sex_type_names <- ifelse(adults_motnp2.2$sex_type == 'F_F', 'female-female',
                                         ifelse(adults_motnp2.2$sex_type == 'F_M', 'male-female', 'male-male'))
(edge_vs_agediff <- 
    ggplot(adults_motnp2.2[adults_motnp2.2$sex_1 != 'U' & adults_motnp2.2$sex_2 != 'U',],
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
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_adults_edge_vs_agediff.pdf',
       plot = edge_vs_agediff,
       width = 30, height = 24, units = 'cm')
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_adults_edge_vs_agediff.pdf',
       plot = edge_vs_agediff,
       width = 30, height = 24, units = 'cm')

#### save data ####
write_csv(draws_motnp2.2, 'data_processed/motnp_bayesian_edgedistributions_a2.b2_22.02.07.csv')
write_csv(tidy_draws, 'data_processed/motnp_bayesian_edgedistributions_tidy_22.02.03.csv')

################ 7) Summarise and plot edge weights --  pooling by relationship type -- WRITE MODEL ################
### build traceplots ####
# random dyads -- many are very wide, high uncertainty
plot(draws_motnp_pool[,plot_cols[1]], type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws_motnp_pool[,plot_cols[2]], col = 'tan')
lines(draws_motnp_pool[,plot_cols[3]], col = 'orange')
lines(draws_motnp_pool[,plot_cols[4]], col = 'green')
lines(draws_motnp_pool[,plot_cols[5]], col = 'chocolate')
lines(draws_motnp_pool[,plot_cols[6]], col = 'blue')
lines(draws_motnp_pool[,plot_cols[7]], col = 'red')
lines(draws_motnp_pool[,plot_cols[8]], col = 'seagreen')
lines(draws_motnp_pool[,plot_cols[9]], col = 'purple')
lines(draws_motnp_pool[,plot_cols[10]],col = 'magenta')
lines(draws_motnp_pool[,plot_cols[11]],col = 'black')
lines(draws_motnp_pool[,plot_cols[12]], col = 'tan')
lines(draws_motnp_pool[,plot_cols[13]], col = 'orange')
lines(draws_motnp_pool[,plot_cols[14]], col = 'green')
lines(draws_motnp_pool[,plot_cols[15]], col = 'chocolate')
lines(draws_motnp_pool[,plot_cols[16]], col = 'blue')
lines(draws_motnp_pool[,plot_cols[17]], col = 'red')
lines(draws_motnp_pool[,plot_cols[18]], col = 'seagreen')
lines(draws_motnp_pool[,plot_cols[19]], col = 'purple')
lines(draws_motnp_pool[,plot_cols[20]],col = 'magenta')
lines(draws_motnp_pool[,plot_cols[21]],col = 'black')
lines(draws_motnp_pool[,plot_cols[22]], col = 'tan')
lines(draws_motnp_pool[,plot_cols[23]], col = 'orange')
lines(draws_motnp_pool[,plot_cols[24]], col = 'green')
lines(draws_motnp_pool[,plot_cols[25]], col = 'chocolate')
lines(draws_motnp_pool[,plot_cols[26]], col = 'blue')
lines(draws_motnp_pool[,plot_cols[27]], col = 'red')
lines(draws_motnp_pool[,plot_cols[28]], col = 'seagreen')
lines(draws_motnp_pool[,plot_cols[29]], col = 'purple')
lines(draws_motnp_pool[,plot_cols[30]],col = 'magenta')

# Sierra = F52, herd members = F60+U21+F98, calf = U17. This looks really good, and not too wide (all seen a good number of times).
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,4000))
lines(draws_motnp_pool$F52_M40, col = 'black')      # non-herd member, adult male
lines(draws_motnp_pool$F52_M15, col = 'tan')        # non-herd member, pubescent male
lines(draws_motnp_pool$F52_M203,col = 'orange')     # non-herd member, adult male
lines(draws_motnp_pool$F52_M26, col = 'green')      # non-herd member, calf
lines(draws_motnp_pool$F52_F8,  col = 'chocolate')  # non-herd member, adult female
lines(draws_motnp_pool$F52_U9,  col = 'blue')       # non-herd member, calf
lines(draws_motnp_pool$F52_F98, col = 'red')        # herd member most frequently absent from sightings
lines(draws_motnp_pool$F52_U17, col = 'purple')     # calf
lines(draws_motnp_pool$F52_U21, col = 'seagreen')   # sister
lines(draws_motnp_pool$F52_F60, col = 'magenta')    # sister's calf

# plot for Sierra (F52)
(traceplot_sierra <- ggplot(tidy_sierra[tidy_sierra$chain == 1,], aes(x = index, y = draw, colour = ID))+
    geom_line()+
    scale_color_viridis_d()+
    theme_classic()+
    scale_x_continuous('MCMC chain position',  expand = c(0,0))+
    scale_y_continuous('edge weight estimate', expand = c(0.02,0)))
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_sierra_traceplot.pdf',
       plot = traceplot_sierra,
       width = 30, height = 24, units = 'cm')

(traceplot_sierra2 <- ggplot(tidy_sierra, aes(x = index, y = draw, colour = as.factor(chain)))+
    geom_line()+
    scale_color_viridis_d()+
    theme_light()+
    facet_wrap(.~ID)+
    scale_x_continuous('MCMC chain position',  expand = c(0,0))+
    scale_y_continuous('edge weight estimate', expand = c(0.02,0))+
    theme(legend.position = c(0.85, 0.15), legend.text = element_text(size = 12),
          legend.title = element_text(size = 14))+
    labs(colour = 'Chain ID'))
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_sierra_traceplot.pdf',
       plot = traceplot_sierra,
       width = 30, height = 24, units = 'cm')

# check chain mixing
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_motnp_pool$`1.weight[1]`, col = 'purple')
lines(draws2_motnp_pool$`1.weight[1]`, col = 'green')
lines(draws3_motnp_pool$`1.weight[1]`, col = 'blue')
lines(draws4_motnp_pool$`1.weight[1]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_motnp_pool$`1.weight[2]`, col = 'purple')
lines(draws2_motnp_pool$`1.weight[2]`, col = 'green')
lines(draws3_motnp_pool$`1.weight[2]`, col = 'blue')
lines(draws4_motnp_pool$`1.weight[2]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
lines(draws1_motnp_pool$`1.weight[3]`, col = 'purple')
lines(draws2_motnp_pool$`1.weight[3]`, col = 'green')
lines(draws3_motnp_pool$`1.weight[3]`, col = 'blue')
lines(draws4_motnp_pool$`1.weight[3]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
which(counts_df$dyad == 'F52_U17') # 41192
lines(draws1_motnp_pool$`1.weight[41192]`, col = 'purple')
lines(draws2_motnp_pool$`1.weight[41192]`, col = 'green')
lines(draws3_motnp_pool$`1.weight[41192]`, col = 'blue')
lines(draws4_motnp_pool$`1.weight[41192]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
which(counts_df$dyad == 'F52_F60') # 40896
lines(draws1_motnp_pool$`1.weight[40896]`, col = 'purple')
lines(draws2_motnp_pool$`1.weight[40896]`, col = 'green')
lines(draws3_motnp_pool$`1.weight[40896]`, col = 'blue')
lines(draws4_motnp_pool$`1.weight[40896]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
which(counts_df$dyad == 'F52_U21') # 41197
lines(draws1_motnp_pool$`1.weight[41197]`, col = 'purple')
lines(draws2_motnp_pool$`1.weight[41197]`, col = 'green')
lines(draws3_motnp_pool$`1.weight[41197]`, col = 'blue')
lines(draws4_motnp_pool$`1.weight[41197]`, col = 'magenta')
plot(NULL, ylim = c(0,1), las = 1, ylab = 'edge weight', xlim = c(0,1000))
which(counts_df$dyad == 'F52_F98') # 40937
lines(draws1_motnp_pool$`1.weight[40937]`, col = 'purple')
lines(draws2_motnp_pool$`1.weight[40937]`, col = 'green')
lines(draws3_motnp_pool$`1.weight[40937]`, col = 'blue')
lines(draws4_motnp_pool$`1.weight[40937]`, col = 'magenta')

### density plots ####
# random columns
dens(draws_motnp_pool[,plot_cols[1]], ylim = c(0,50), xlim = c(0,1), las = 1)
dens(add = T, draws_motnp_pool[,plot_cols[2]], col = 'tan')
dens(add = T, draws_motnp_pool[,plot_cols[3]], col = 'orange')
dens(add = T, draws_motnp_pool[,plot_cols[4]], col = 'green')
dens(add = T, draws_motnp_pool[,plot_cols[5]], col = 'chocolate')
dens(add = T, draws_motnp_pool[,plot_cols[6]], col = 'blue')
dens(add = T, draws_motnp_pool[,plot_cols[7]], col = 'red')
dens(add = T, draws_motnp_pool[,plot_cols[8]], col = 'seagreen')
dens(add = T, draws_motnp_pool[,plot_cols[9]], col = 'purple')
dens(add = T, draws_motnp_pool[,plot_cols[10]],col = 'magenta')
dens(add = T, draws_motnp_pool[,plot_cols[11]],col = 'black')
dens(add = T, draws_motnp_pool[,plot_cols[12]],col = 'tan')
dens(add = T, draws_motnp_pool[,plot_cols[13]],col = 'orange')
dens(add = T, draws_motnp_pool[,plot_cols[14]],col = 'green')
dens(add = T, draws_motnp_pool[,plot_cols[15]],col = 'chocolate')
dens(add = T, draws_motnp_pool[,plot_cols[16]],col = 'blue')
dens(add = T, draws_motnp_pool[,plot_cols[17]],col = 'red')
dens(add = T, draws_motnp_pool[,plot_cols[18]],col = 'seagreen')
dens(add = T, draws_motnp_pool[,plot_cols[19]],col = 'purple')
dens(add = T, draws_motnp_pool[,plot_cols[20]],col = 'magenta')
dens(add = T, draws_motnp_pool[,plot_cols[21]],col = 'black')
dens(add = T, draws_motnp_pool[,plot_cols[22]],col = 'tan')
dens(add = T, draws_motnp_pool[,plot_cols[23]],col = 'orange')
dens(add = T, draws_motnp_pool[,plot_cols[24]],col = 'green')
dens(add = T, draws_motnp_pool[,plot_cols[25]],col = 'chocolate')
dens(add = T, draws_motnp_pool[,plot_cols[26]],col = 'blue')
dens(add = T, draws_motnp_pool[,plot_cols[27]],col = 'red')
dens(add = T, draws_motnp_pool[,plot_cols[28]],col = 'seagreen')
dens(add = T, draws_motnp_pool[,plot_cols[29]],col = 'purple')
dens(add = T, draws_motnp_pool[,plot_cols[30]],col = 'magenta')
dens(add = T, draws_motnp_pool[,which(counts_df$dyad == 'F52_U17')+1],col = 'magenta')
dens(add = T, draws_motnp_pool[,which(counts_df$dyad == 'F52_F60')+1],col = 'purple')
dens(add = T, draws_motnp_pool[,which(counts_df$dyad == 'F52_U21')+1],col = 'blue')
dens(add = T, draws_motnp_pool[,which(counts_df$dyad == 'F52_F98')+1],col = 'red')

# density plots
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
means <- data.frame(dyad = counts_df$dyad,
                    mean = apply(draws[,2:106954], 2, mean),
                    median = apply(draws[,2:106954], 2, median))
draws_pi <- apply(draws[2:106954], 2, rethinking::PI)
means$pi_lwr <- draws_pi[1,] ; means$pi_upr <- draws_pi[2,]

# organise dem_class
plot_data <- left_join(x = means, y = counts_df, by = 'dyad')
head(plot_data)
plot_data$dem_class_1_cat <- as.integer(as.factor(plot_data$dem_class_1))
plot_data$dem_class_2_cat <- as.integer(as.factor(plot_data$dem_class_2))
plot_data$dem_type_short <- ifelse(plot_data$dem_class_1_cat <= plot_data$dem_class_2_cat,
                                   paste(plot_data$dem_class_1, plot_data$dem_class_2, sep = '_'),
                                   paste(plot_data$dem_class_2, plot_data$dem_class_1, sep = '_'))
sort(unique(plot_data$dem_type_short))
types <- data.frame(all = sort(unique(plot_data$dem_type_short)), type1 = NA, type2 = NA)
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
plot_data <- left_join(plot_data, types, by = 'dem_type_short')
plot_data$join <- paste(plot_data$join, ' ', sep = ' ')
which(is.na(plot_data$join))

# boxplot edge weights by demographic type of dyad -- all types
colours <- c('magenta','purple','grey','grey','magenta','purple','grey','blue','purple','purple','purple','blue','purple','grey','grey','grey','grey','grey','magenta','purple','grey','purple','purple','blue','purple','grey','grey')
(edge_vs_demtype_all <- 
    ggplot(data = plot_data, aes(y = mean, x = join))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 fill = colours)+
    theme_classic()+
    labs(x = 'dyad demography (A/P/J/C = adult/pubescent/juvenile/calf, M/F/U = male/female/unknown',
         y = 'mean edge weight')+
    coord_flip())
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_all_edge_vs_demtype.pdf',
       plot = edge_vs_demtype_all,
       width = 20, height = 24, units = 'cm')  
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_all_edge_vs_demtype.pdf',
       plot = edge_vs_demtype_all,
       width = 20, height = 24, units = 'cm')  

# boxplot edge weights by demographic type of dyad -- only adults/pubescents
adults <- plot_data[plot_data$age_cat_id_1 > 2 & plot_data$age_cat_id_2 > 2,]
adults <- adults[!is.na(adults$dyad),]
colours <- c('magenta','purple','magenta','purple','grey','blue','purple','blue','purple','magenta','purple','grey','blue','purple')
(edge_vs_demtype <- 
    ggplot(data = adults, aes(y = mean, x = join))+
    geom_boxplot(notch = T, outlier.shape = 4,
                 fill = colours)+
    theme_classic()+
    labs(x = 'dyad demography (A/P = adult/pubescent, M/F/U = male/female/unknown',
         y = 'mean edge weight')+
    coord_flip())
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_adults_edge_vs_demtype.pdf',
       plot = edge_vs_demtype,
       width = 20, height = 24, units = 'cm')  
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_adults_edge_vs_demtype.pdf',
       plot = edge_vs_demtype,
       width = 20, height = 24, units = 'cm')

# scatter plot of all adult edges by age difference
adults$sex_type_names <- ifelse(adults$sex_type == 'F_F', 'female-female',
                                ifelse(adults$sex_type == 'F_M', 'male-female', 'male-male'))
(edge_vs_agediff <- 
    ggplot(adults[adults$sex_1 != 'U' & adults$sex_2 != 'U',],
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
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_adults_edge_vs_agediff.pdf',
       plot = edge_vs_agediff,
       width = 30, height = 24, units = 'cm')
ggsave(path = 'outputs/motnp_22.02.01/', filename = 'motnp_adults_edge_vs_agediff.pdf',
       plot = edge_vs_agediff,
       width = 30, height = 24, units = 'cm')

#### save data ####
write_csv(draws, 'data_processed/motnp_bayesian_edgedistributions_pooled_22.02.03.csv')
write_csv(tidy_draws, 'data_processed/motnp_bayesian_edgedistributions_tidy_22.02.03.csv')

################ 8) Create network plots ################
head(means_motnp2.2)
length(unique(plot_data_motnp2.2$id_1))+1 # number of individuals

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
dyads <- colnames(draws_motnp2.2)
#which(dyads == "F52_F60" | dyads == "F52_F8" | dyads == "F52_F98" | dyads == "F52_M15" | dyads == "F52_M26" | dyads == "F52_M40" |
#        dyads == "F52_U17" | dyads == "F52_U21" | dyads == "F60_F8"  | dyads == "F60_F98" | dyads == "F60_M15" | dyads == "F60_M26" |
#        dyads == "F60_M40" | dyads == "F60_U17" | dyads == "F60_U21" | dyads == "F8_F98"  | dyads == "F8_M15"  | dyads == "F8_M26"  |
#        dyads == "F8_M40"  | dyads == "F8_U17"  | dyads == "F8_U21"  | dyads == "F98_M15" | dyads == "F98_M26" | dyads == "F98_M40" | 
#        dyads == "F98_U17" | dyads == "F98_U21" | dyads == "M15_M26" | dyads == "M15_M40" | dyads == "M15_U17" | dyads == "M15_U21" |
#        dyads == "M26_M40" | dyads == "M26_U17" | dyads == "M26_U21" | dyads == "M40_U17" | dyads == "M40_U21" | dyads == "U17_U21")
#draws_test <- draws[,c(40897,40918,40938,40992,41104,41120,41193,41198,44140,44160,44214,44326,44342,44415,44420,51363,51417,51529,
#                       51545,51618,51623,57867,57979,57995,58068,58073,73396,73412,73485,73490,96092,96165,96170,98381,98386,105248)]

draws_test <- draws_motnp2.2[,which(dyads == "F52_F60" | dyads == "F52_F8" | dyads == "F52_F98" | dyads == "F52_M15" |
                                      dyads == "F52_M26" | dyads == "F52_M40" | dyads == "F52_U17" | dyads == "F52_U21" |
                                      dyads == "F60_F8"  | dyads == "F60_F98" | dyads == "F60_M15" | dyads == "F60_M26" |
                                      dyads == "F60_M40" | dyads == "F60_U17" | dyads == "F60_U21" | dyads == "F8_F98"  |
                                      dyads == "F8_M15"  | dyads == "F8_M26"  | dyads == "F8_M40"  | dyads == "F8_U17"  |
                                      dyads == "F8_U21"  | dyads == "F98_M15" | dyads == "F98_M26" | dyads == "F98_M40" |
                                      dyads == "F98_U17" | dyads == "F98_U21" | dyads == "M15_M26" | dyads == "M15_M40" |
                                      dyads == "M15_U17" | dyads == "M15_U21" | dyads == "M26_M40" | dyads == "M26_U17" |
                                      dyads == "M26_U21" | dyads == "M40_U17" | dyads == "M40_U21" | dyads == "U17_U21")]
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
adj_tensor_test[, 1, ] # Print the first sample of the posterior distribution over adjacency matrices
adj_tensor_test[1, , ] # Print the first sample of the posterior distribution over adjacency matrices

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

# work out how to filter data for plotting
apply(draws_test, 2, median)
# F52_F60 = 0.863223000 ; 52_F8 = 0.033248300 ; F52_F98 = 0.623663500 ; F52_M15 = 0.009397185 ; F52_M26 = 0.034015850 ; F52_M40 = 0.054451600 ; F52_U17 = 0.935521500 ; F52_U21 = 0.837885000 ; F60_F8 = 0.035207550 
adj_mid
#      [,1]     [,2]       [,3]      [,4]        [,5]       [,6]       [,7]        [,8]       [,9]
# [1,]    0 0.863223 0.03324830 0.6236635 0.009397185 0.03401585 0.05445160 0.935521500 0.83788500
# [2,]    0 0.000000 0.03520755 0.5926630 0.009502835 0.03659575 0.05712505 0.820356000 0.95435200
## 1 = F52, 2 = F60, 3 = F8, 4 = F98, 5 = M15, 6 = M26, 7 = M40, 8 = U17, 9 = U21 --> They go in alphanumeric order

### Plot for all... this isn't going to work! There are too many to see anything! ####
unique(counts_df$node_1)         # there are 462 individuals in node_1, but the numbers go up to 471
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))     # new variable reassigning dyad number 1-462
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1   # new variable reassigning dyad number 2-463

adj_tensor <- array(0, c(463, 463, 4000))
for (dyad_id in 1:nrow(counts_df)) {
  dyad_row <- counts_df[counts_df$dyad_id == dyad_id, ]
  adj_tensor[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws[, dyad_id+1]
}
adj_tensor[, , 1]

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])

# Calculate standardised width/range of credible intervals.
adj_range <- ((adj_upper - adj_lower)/adj_mid)
adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Plot the median graph first and then the standardised width graph to show uncertainty over edges.
coords <- igraph::layout_nicely(g_mid)
plot(g_mid, edge.width = E(g_mid)$weight, edge.color = "black",  layout = coords)
plot(g_mid, edge.width = E(g_rng)$weight, edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = c(unique(counts_df$id_1),'U9'),
     vertex.label.dist = 0, vertex.label.color="black", layout=coords, add=TRUE)

### plot males only  -- still too many to see anything! ####
counts_df_male <- counts_df[counts_df$dem_class_1 == 'AM' | counts_df$dem_class_1 == 'PM',]
counts_df_male <- counts_df_male[counts_df_male$dem_class_2 == 'AM' | counts_df_male$dem_class_2 == 'PM',]
male_dyads <- counts_df_male$dyad
dyads <- colnames(draws_motnp2.2)
male_dyads  # 23653 dyads -- can't subset by hand!
draws_male <- draws_motnp2.2[,male_dyads]
colnames(draws_male)

adj_tensor_male <- array(0, c(length(unique(counts_df_male$id_1))+1,
                              length(unique(counts_df_male$id_1))+1,
                              length(draws_male$M1_M10)))
counts_df_male$dyad_id2 <- as.integer(as.factor(counts_df_male$dyad_id))
counts_df_male$node_1_nogaps <- as.integer(as.factor(counts_df_male$node_1))
counts_df_male$node_2_nogaps <- as.integer(as.factor(counts_df_male$node_2))+1
for (dyad_id in 1:nrow(counts_df_male)) {
  dyad_row <- counts_df_male[counts_df_male$dyad_id2 == dyad_id, ]
  adj_tensor_male[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws_male[, dyad_id]
}
adj_tensor_male[, , 1]

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles_male <- apply(adj_tensor_male, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower_male <- adj_quantiles_male[1, , ])
(adj_mid_male   <- adj_quantiles_male[2, , ])
(adj_upper_male <- adj_quantiles_male[3, , ])

# Calculate standardised width/range of credible intervals.
(adj_range_male <- ((adj_upper_male - adj_lower_male)/adj_mid_male))
(adj_range_male[is.nan(adj_range_male)] <- 0)

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid_m <- graph_from_adjacency_matrix(adj_mid_male,   mode="undirected", weighted=TRUE)
g_rng_m <- graph_from_adjacency_matrix(adj_range_male, mode="undirected", weighted=TRUE)

# Plot the median graph first and then the standardised width graph to show uncertainty over edges.
coords <- igraph::layout_nicely(g_mid_m)
ages <- counts_df_male[,c('id_1','age_cat_id_1')] %>% 
  distinct()
M99 <- data.frame(id_1 = 'M99', age_cat_id_1 = 6)
ages <- rbind(ages, M99) ; head(ages)

unique(ages$age_cat_id_1) # 2-7

plot(g_mid_m, edge.width = E(g_mid_m)$weight, edge.color = "black",  layout = coords)
plot(g_mid_m, edge.width = E(g_rng_m)$weight, edge.color = rgb(0, 0, 0, 0.05), 
     vertex.label = ages$id_1,
     vertex.color = ifelse(ages$age_cat_id_1 == 2, 'yellow',
                           ifelse(ages$age_cat_id_1 == 3, 'orange',
                                  ifelse(ages$age_cat_id_1 == 4, 'red',
                                         ifelse(ages$age_cat_id_1 == 5, 'purple',
                                                ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen'))))),
     vertex.label.dist = 0,
     vertex.label.color = ifelse(ages$age_cat_id_1 == 6, 'white','black'),
     layout = coords, add = TRUE)

ages <- separate(ages, id_1, c('m','num'), 1, F)
plot(g_mid_m, edge.width = E(g_mid_m)$weight, layout = coords,
     edge.color = 'black',
     vertex.size = 7)
plot(g_mid_m, edge.width = E(g_rng_m)$weight, edge.color = rgb(0, 0, 0, 0.05), 
     vertex.size = 7,
     vertex.color = ifelse(ages$age_cat_id_1 == 2, 'yellow',
                           ifelse(ages$age_cat_id_1 == 3, 'orange',
                                  ifelse(ages$age_cat_id_1 == 4, 'red',
                                         ifelse(ages$age_cat_id_1 == 5, 'purple',
                                                ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen'))))),
     vertex.label = ages$num,
     vertex.label.dist = 0,
     vertex.label.color = ifelse(ages$age_cat_id_1 == 6, 'white','black'),
     layout = coords, add = TRUE)

plot(g_mid_m, edge.width = ifelse(E(g_mid_m)$weight<0.3, 0,
                                  E(g_mid_m)$weight),
     layout = coords,
     edge.color = 'black',
     vertex.size = 7)
plot(g_mid_m, edge.width = ifelse(E(g_mid_m)$weight<0.3, 0,
                                  E(g_rng_m)$weight),
     edge.color = rgb(0, 0, 0, 0.25), 
     vertex.size = 7,
     vertex.color = ifelse(ages$age_cat_id_1 == 2, 'yellow',
                           ifelse(ages$age_cat_id_1 == 3, 'orange',
                                  ifelse(ages$age_cat_id_1 == 4, 'red',
                                         ifelse(ages$age_cat_id_1 == 5, 'purple',
                                                ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen'))))),
     vertex.label = ages$num,
     vertex.label.dist = 0,
     vertex.label.color = ifelse(ages$age_cat_id_1 == 6, 'white','black'),
     layout = coords, add = TRUE)

### plot ADULT males only  -- STILL too many to see anything! ####
counts_df_am <- counts_df[counts_df$dem_class_1 == 'AM' & counts_df$dem_class_2 == 'AM',]
am_dyads <- counts_df_am$dyad
am_dyads  # 23653 dyads -- can't subset by hand!
draws_motnp2.2_am <- draws_motnp2.2[,am_dyads]
colnames(draws_motnp2.2_am)
View(draws_motnp2.2_am)
meds_am <- as.data.frame(apply(draws_motnp2.2_am, 2, median))
summary(meds_am$`apply(draws_motnp2.2_am, 2, median)`)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.01995 0.06660 0.10206 0.12150 0.15869 0.62023 -- much more as expected

adj_tensor_am <- array(0, c(length(unique(counts_df_am$id_1))+1,
                            length(unique(counts_df_am$id_1))+1,
                            length(draws_motnp2.2_am$M1_M10)))
counts_df_am$dyad_id2 <- as.integer(as.factor(counts_df_am$dyad_id))
counts_df_am$node_1_nogaps <- as.integer(as.factor(counts_df_am$node_1))
counts_df_am$node_2_nogaps <- as.integer(as.factor(counts_df_am$node_2))+1
for (dyad_id in 1:nrow(counts_df_am)) {
  dyad_row <- counts_df_am[counts_df_am$dyad_id2 == dyad_id, ]
  adj_tensor_am[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws_motnp2.2_am[, dyad_id+1]
}
adj_tensor_am[, , 1]

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles_am <- apply(adj_tensor_am, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower_am <- adj_quantiles_am[1, , ])
(adj_mid_am   <- adj_quantiles_am[2, , ])
(adj_upper_am <- adj_quantiles_am[3, , ])

max(adj_mid_am)

# Calculate standardised width/range of credible intervals.
(adj_range_am <- ((adj_upper_am - adj_lower_am)/adj_mid_am))
adj_range_am[is.nan(adj_range_am)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid_am <- graph_from_adjacency_matrix(adj_mid_am,   mode="undirected", weighted=TRUE)
g_rng_am <- graph_from_adjacency_matrix(adj_range_am, mode="undirected", weighted=TRUE)

# Plot the median graph first and then the standardised width graph to show uncertainty over edges.
coords <- igraph::layout_nicely(g_mid_am)
ages <- counts_df_am[,c('id_1','age_cat_id_1')] %>% 
  distinct()
M99 <- data.frame(id_1 = 'M99', age_cat_id_1 = 6)
ages <- rbind(ages, M99) ; head(ages)

unique(ages$age_cat_id_1) # 5-7

plot(g_mid_am, edge.width = E(g_mid_am)$weight, edge.color = "black",  layout = coords)
plot(g_mid_am, edge.width = E(g_rng_am)$weight, edge.color = rgb(0, 0, 0, 0.05), 
     vertex.label = ages$id_1,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     vertex.label.color = 'white',
     layout = coords, add = TRUE)

plot(g_mid_am, edge.width = E(g_mid_am)$weight, edge.color = "black",  layout = coords, vertex.size = 3, vertex.label.color = 'transparent')
plot(g_mid_am, edge.width = E(g_rng_am)$weight, edge.color = rgb(0, 0, 0, 0.05), 
     vertex.label = ages$id_1,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'red',
                           ifelse(ages$age_cat_id_1 == 6, 'magenta', 'green')),
     vertex.label.dist = 0.5,
     vertex.label.color = ifelse(ages$age_cat_id_1 == 5, 'red',
                                 ifelse(ages$age_cat_id_1 == 6, 'magenta', 'green')),
     vertex.size = 3,
     layout = coords, add = TRUE)

ages <- separate(ages, id_1, into = c('M', 'num'), remove = F, sep = 1)
summary(E(g_mid_am)$weight)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01995 0.06661 0.10207 0.12151 0.15871 0.62023 
plot(g_mid_am, edge.width = E(g_mid_am)$weight,
     edge.color = ifelse(E(g_mid_am)$weight < 0.16, 'transparent','black'),
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords)
plot(g_mid_am, edge.width = E(g_rng_am)$weight,
     edge.color = ifelse(E(g_mid_am)$weight < 0.16, 'transparent',
                         rgb(0, 0, 0, 0.05)), 
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords, add = TRUE)

plot(g_mid_am, edge.width = 3*E(g_mid_am)$weight,
     edge.color = ifelse(E(g_mid_am)$weight < 0.3, 'transparent','black'),
     layout = coords, vertex.size = 3, vertex.label.color = 'transparent')
plot(g_mid_am, edge.width = 3*E(g_rng_am)$weight,
     edge.color = ifelse(E(g_mid_am)$weight < 0.3, 'transparent',
                         rgb(0, 0, 0, 0.25)), 
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords, add = TRUE)

#### Create igraph object for all elephants and then just plot certain individuals ####
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))     # new variable reassigning dyad number 1-462
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1   # new variable reassigning dyad number 2-463

adj_tensor <- array(0, c(463, 463, 4000))
for (dyad_id in 1:nrow(counts_df)) {
  dyad_row <- counts_df[counts_df$dyad_id == dyad_id, ]
  adj_tensor[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws_motnp2.2[, dyad_id+1]
}
adj_tensor[, , 1]

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])

# Calculate standardised width/range of credible intervals.
adj_range <- ((adj_upper - adj_lower)/adj_mid)
adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode = "undirected", weighted = TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode = "undirected", weighted = TRUE)

# Plot all
coords <- igraph::layout_nicely(g_mid)
plot(g_mid, edge.width = E(g_mid)$weight, edge.color = "black", layout = coords)
plot(g_mid, edge.width = E(g_rng)$weight, edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = c(unique(counts_df$id_1),'U9'),
     vertex.label.dist = 0, vertex.label.color="black", layout=coords, add=TRUE)

# Plot males
which(counts_df$sex_1 == 'M' & counts_df$sex_2 == 'M')
plot(g_mid[which(counts_df$sex_1 == 'M' & counts_df$sex_2 == 'M')],
     edge.width = E(g_mid)$weight, edge.color = "black", layout = coords)

eles <- rbind(counts_df[,c('id_1','sex_1','age_cat_id_1','node_1','node_1_nogaps')],
              c('U9','U','1','472','463')) %>%
  distinct()
which(eles$sex_1 == 'M')
plot(g_mid[which(eles$sex_1 == 'M')],
     edge.width = E(g_mid)$weight, edge.color = "black", layout = coords)
warnings()
plot(g_mid, edge.width = E(g_rng)$weight, edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = c(unique(counts_df$id_1),'U9'),
     vertex.label.dist = 0, vertex.label.color="black",
     layout=coords, add=TRUE)

delete_test <- g_mid_am %>% delete.edges(1:1000)
plot(delete_test)
delete_test <- g_mid %>% delete.edges(1:5000)
plot(delete_test)
length(counts_df$dyad_id)
delete_test <- g_mid %>% delete.edges(102750)   # this works
delete_test <- g_mid %>% delete.edges(102751)   # this doesn't -- only seems to think there are 102750 edges, not 106953
# length(which(g_mid == 0))
delete_test <- g_mid %>% delete.edges(1:102750) # this removes all edges
plot(delete_test, vertex.size = 2)


m <- g_mid %>% delete.edges(which(counts_df$sex_1 == 'F' | counts_df$sex_2 == 'F' | counts_df$sex_1 == 'U' | counts_df$sex_2 == 'U')) # this doesn't work because includes values > 102750
length(g_mid)

test_delete <- g_mid_am %>% delete.edges(which(E(g_mid_am)$weight < 0.16))
test_delete_rng <- g_rng_am %>% delete.edges(which(E(g_mid_am)$weight < 0.16))

plot(test_delete, edge.width = E(test_delete)$weight,
     edge.color = 'black',
     vertex.color = 'transparent',
     vertex.size = 8,
     layout = coords)
plot(test_delete, edge.width = E(test_delete_rng)$weight,
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     edge.color = rgb(0, 0, 0, 0.05), 
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords, add = TRUE)

test_delete <- g_mid_am %>% delete.edges(which(E(g_mid_am)$weight < 0.4))
test_delete_rng <- g_rng_am %>% delete.edges(which(E(g_mid_am)$weight < 0.4))
test_delete2 <- test_delete %>% igraph::delete.vertices(which(E(test_delete)$weight < 0.4))
test_delete_rng2 <- test_delete %>% igraph::delete.vertices(which(E(test_delete)$weight < 0.16))

plot(test_delete2, edge.width = E(test_delete2)$weight,
     edge.color = 'black',
     vertex.color = 'transparent',
     vertex.size = 8,
     layout = coords)
plot(test_delete2, edge.width = E(test_delete_rng)$weight,
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     edge.color = rgb(0, 0, 0, 0.05), 
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords, add = TRUE)

test_delete <- g_mid_am %>% delete.edges(which(E(g_mid_am)$weight < 0.16))
test_delete_rng <- g_rng_am %>% delete.edges(which(E(g_mid_am)$weight < 0.16))
test_delete2 <- test_delete %>% delete_vertices(104)
test_delete_rng2 <- test_delete %>% delete_vertices(104)
plot(test_delete, edge.width = E(test_delete)$weight,
     edge.color = 'black',
     vertex.color = 'white',
     vertex.size = 8,
     layout = coords)
plot(test_delete, edge.width = E(test_delete_rng)$weight,
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     edge.color = rgb(0, 0, 0, 0.05), 
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords, add = TRUE)
### can get the delete.edges function to work but not delete_vertices -- this seems to add the edges back in and also not delete the ones I tell it to delete

which(counts_df$sex_1 == 'M' & counts_df$sex_2 == 'M')
eles <- data.frame(id = c('F1',counts_df$id_2[1:462]),
                   sex = c('F',counts_df$sex_2[1:462]))
eles
plot(g_mid[which(eles$sex == 'M'),which(eles$sex == 'M')],
     edge.width = E(g_mid)$weight[which(counts_df$sex_1 == 'M' & counts_df$sex_2 == 'M')], edge.color = "black", layout = coords)
warnings()

plot(g_mid_test,
     vertex.color = 'transparent',
     vertex.label.color = 'transparent',
     edge.color = ifelse(E(g_mid_test)$weight < 0.5, 'transparent','black'),
     edge.weight = E(g_mid_test)$weight,
     lty = 2)

Isolated <- which(degree(g_mid_test) < 2)
G2 <- delete.vertices(g_mid_test, Isolated)
plot(G2)

E(g_mid_test)$weight <- ifelse(E(g_mid_test)$weight < 0.5, 0, E(g_mid_test)$weight)
coords <- layout_nicely(g_mid_test)
plot(g_mid_test,
     vertex.color = 'transparent',
     vertex.label.color = 'transparent',
     edge.color = ifelse(E(g_mid_test)$weight < 0.5, 'transparent','black'),
     edge.weight = E(g_mid_test)$weight,
     lty = 2,
     layout = coords)

# Recalculate medians and distributions
adj_quantiles_test <- apply(adj_tensor_test, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower_test <- adj_quantiles_test[1, , ])
(adj_mid_test   <- adj_quantiles_test[2, , ])
(adj_upper_test <- adj_quantiles_test[3, , ])
adj_range_test <- ((adj_upper_test - adj_lower_test)/adj_mid_test)
adj_range_test[is.nan(adj_range_test)] <- 0

# FILTER SMALL WEIGHTS
adj_mid_test_0 <- ifelse(adj_mid_test < 0.4, 0, adj_mid_test)
adj_rng_test_0 <- ifelse(adj_mid_test < 0.4, 0, adj_range_test)

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid_test_0 <- graph_from_adjacency_matrix(adj_mid_test_0,   mode = "undirected", weighted = TRUE)
g_rng_test_0 <- graph_from_adjacency_matrix(adj_rng_test_0, mode = "undirected", weighted = TRUE)

plot(g_mid_test_0)

Isolated <- which(degree(g_mid_test_0)==0)
G2 <- delete.vertices(g_mid_test_0, Isolated)
plot(G2)

Isolated <- which(sum(E(g_mid_test)$weight) < 0.5)
G2 <- delete.vertices(g_mid_test, Isolated)
plot(G2)

### SO...
# it works to go each time right from the initial data
# it works to filter each time at the graph_from_adjacency_matrix() stage
# delete_vertices by connection works after filtering at  graph_from_adjacency_matrix() stage 
# delete_vertices by connection doesn't work without filtering at  graph_from_adjacency_matrix() stage 

### start over...
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))     # new variable reassigning dyad number 1-462
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1   # new variable reassigning dyad number 2-463

adj_tensor <- array(0, c(463, 463, 4000))
for (dyad_id in 1:nrow(counts_df)) {
  dyad_row <- counts_df[counts_df$dyad_id == dyad_id, ]
  adj_tensor[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws_motnp2.2[, dyad_id+1]
}
adj_tensor[, , 1]

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])

# Calculate standardised width/range of credible intervals.
adj_range <- ((adj_upper - adj_lower)/adj_mid)
adj_range[is.nan(adj_range)] <- 0

# FILTER WEAK RELATIONSHIPS
adj_mid_0.3 <- ifelse(adj_mid < 0.3, 0, adj_mid)
adj_rng_0.3 <- ifelse(adj_mid < 0.3, 0, adj_range)

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid_0.3 <- graph_from_adjacency_matrix(adj_mid_0.3,   mode = "undirected", weighted = TRUE)
g_rng_0.3 <- graph_from_adjacency_matrix(adj_rng_0.3, mode = "undirected", weighted = TRUE)

# FILTER UNCONNECTED NODES
Isolated_0.3 <- which(degree(g_mid_0.3)==0)
filt_mid_0.3 <- delete.vertices(g_mid_0.3, Isolated_0.3)
length(Isolated_0.3)

# plot
coords_0.3 <- layout_nicely(filt_mid_0.3)
eles_0.3 <- eles[-Isolated_0.3,]
eles_0.3 <- separate(eles_0.3, id, c('sex','num'), 1, F)
eles_0.3
plot(filt_mid_0.3,
     edge.width = 2*E(g_mid_0.3)$weight,
     edge.color = 'black',
     vertex.size = 8,
     vertex.color = ifelse(eles_0.3$sex == 'M','blue',
                           ifelse(eles_0.3$sex == 'F','magenta','grey')),
     layout = coords_0.3)
plot(filt_mid_0.3,
     edge.width = 2*E(g_rng_0.3)$weight,
     edge.color = rgb(0,0,0,0.2),
     vertex.size = 8,
     vertex.label.color = 'white',
     vertex.label = eles_0.3$num,
     vertex.label.family = 'Helvetica',
     vertex.color = ifelse(eles_0.3$sex == 'M','blue',
                           ifelse(eles_0.3$sex == 'F','magenta','grey')),
     layout = coords_0.3,
     add = T)

coords2 <- layout.kamada.kawai(filt_mid_0.3)
plot(filt_mid_0.3,
     edge.width = 2*E(g_mid_0.3)$weight,
     edge.color = 'black',
     vertex.size = 8,
     vertex.color = ifelse(eles_0.3$sex == 'M','blue',
                           ifelse(eles_0.3$sex == 'F','magenta','grey')),
     layout = coords2)
plot(filt_mid_0.3,
     edge.width = 2*E(g_rng_0.3)$weight,
     edge.color = rgb(0,0,0,0.2),
     vertex.size = 8,
     vertex.label.color = 'white',
     vertex.label = eles_0.3$num,
     vertex.label.family = 'Helvetica',
     vertex.color = ifelse(eles_0.3$sex == 'M','blue',
                           ifelse(eles_0.3$sex == 'F','magenta','grey')),
     layout = coords2,
     add = T)

### Filter higher
# FILTER WEAK RELATIONSHIPS
adj_mid_0.5 <- array(0, c(463, 463))
for (i in 1:nrow(counts_df)) {
  dyad_row <- counts_df[counts_df$dyad_id == i, ]
  adj_mid_0.5[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps] <- ifelse(adj_mid[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps] < 0.5, 0, adj_mid[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps])
}

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid_0.5 <- graph_from_adjacency_matrix(adj_mid_0.5,   mode = "undirected", weighted = TRUE)
g_rng_0.5 <- graph_from_adjacency_matrix(adj_rng_0.5, mode = "undirected", weighted = TRUE)

# FILTER UNCONNECTED NODES
Isolated_0.5 <- which(degree(g_mid_0.5)==0)
filt_mid_0.5 <- delete.vertices(g_mid_0.5, Isolated_0.5)
length(Isolated_0.5)

# plot
eles_0.5 <- eles[-Isolated_0.5,]
eles_0.5 <- separate(eles_0.5, id, c('sex','num'), 1, F)
eles_0.5
coords_0.5 <- layout.kamada.kawai(filt_mid_0.5)
plot(filt_mid_0.5,
     edge.width = 2*E(g_mid_0.5)$weight,
     edge.color = 'black',
     vertex.size = 8,
     vertex.color = ifelse(eles_0.5$sex == 'M','blue',
                           ifelse(eles_0.5$sex == 'F','magenta','grey')),
     layout = coords_0.5)
plot(filt_mid_0.5,
     edge.width = 2*E(g_rng_0.5)$weight,
     edge.color = rgb(0,0,0,0.2),
     vertex.size = 8,
     #vertex.label.color = 'white',
     #vertex.label = eles_0.5$num,
     vertex.label.family = 'Helvetica',
     #vertex.color = ifelse(eles_0.5$sex == 'M','blue',
      #                     ifelse(eles_0.5$sex == 'F','magenta','grey')),
     layout = coords_0.5,
     add = T)

counts_df[counts_df$dyad == 'F1_F112',]

## check that these values actually make sense -- no they don't...
counts_df[counts_df$dyad == 'F105_M100',] # never observed together
counts_df[counts_df$dyad == 'F105_M60',]  # never observed together
counts_df[counts_df$dyad == 'F102_U8',]   # together once
counts_df[counts_df$dyad == 'F102_U63',]  # never observed together
counts_df[counts_df$dyad == 'F102_F21',]  # never observed together
counts_df[counts_df$dyad == 'F1_F57',]    # never observed together
counts_df[counts_df$dyad == 'F1_U41',]    # never observed together

which(E(g_mid)$weight > 0.5) # should be 306 connections still

# start over again... ####
adj_tensor <- array(0, c(463, 463, 4000))
for (dyad_id in 1:nrow(counts_df)) {
  dyad_row <- counts_df[counts_df$dyad_id == dyad_id, ]
  adj_tensor[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws_motnp2.2[, dyad_id+1]
}
adj_tensor[, , 1]

adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- ((adj_upper - adj_lower)/adj_mid)
adj_range[is.nan(adj_range)] <- 0

ele_links_mid <- as.data.frame(adj_mid)
colnames(ele_links_mid) <- c(sort(unique(counts_df$id_1)),'U9')
rownames(ele_links_mid) <- c(sort(unique(counts_df$id_1)),'U9')
ele_links_long_m <- pivot_longer(ele_links_mid, cols = everything(), names_to = 'id_2', values_to = 'mid_weight')
ele_links_long_m$id_1 <- rep(c(sort(unique(counts_df$id_1)),'U9'), each = 463)

ele_links_rng <- as.data.frame(adj_range)
colnames(ele_links_rng) <- c(sort(unique(counts_df$id_1)),'U9')
rownames(ele_links_rng) <- c(sort(unique(counts_df$id_1)),'U9')
ele_links_long_r <- pivot_longer(ele_links_rng, cols = everything(), names_to = 'id_2', values_to = 'rng_weight') # produces something 2 times longer (+463) than it should be (should be 106953)
ele_links_long_r$id_1 <- rep(c(sort(unique(counts_df$id_1)),'U9'), each = 463)

ele_links <- left_join(ele_links_long_m, ele_links_long_r)
ele_links <- ele_links[ele_links$id_2 != ele_links$id_1,]
ele_links$node_1 <- as.integer(as.factor(ele_links$id_1))
ele_links$node_2 <- as.integer(as.factor(ele_links$id_2))
ele_links$dyad_sort <- ifelse(ele_links$node_1 > ele_links$node_2,
                              paste(ele_links$id_2, ele_links$id_1, sep = '_'),
                              paste(ele_links$id_1, ele_links$id_2, sep = '_'))
ele_links$dyad_unsort <- paste(ele_links$id_1, ele_links$id_2, sep = '_')
table(ele_links$dyad)

links <- data.frame(dyad = sort(unique(ele_links$dyad))) %>% 
  separate(col = 'dyad', into = c('id_1','id_2'), sep = '_', remove = F)
links$mid <- NA ; links$rng <- NA
for(i in 1:nrow(links)){
  dyad_row <- ele_links[ele_links$dyad == links$dyad[i],]
  links$mid[i] <- dyad_row$mid_weight[which(dyad_row$dyad == dyad_row$dyad_unsort)]
  links$rng[i] <- dyad_row$rng_weight[which(dyad_row$dyad == dyad_row$dyad_unsort)]
}

weights <- left_join(links, counts_df, by = 'dyad')
min(weights$mid)
length(which(weights$mid == 0)) # 4203 -- everything from M76_U55 onwards

# traceplots -- all show reasonabe mixing and definitely above 0.00!!!
plot(draws_motnp2.2$M76_U55[1:1000], type = 'l', ylim = c(0,1))
lines(draws_motnp2.2$M76_U55[1001:2000], col = 'blue')
lines(draws_motnp2.2$M76_U55[2001:3000], col = 'red')
lines(draws_motnp2.2$M76_U55[3001:4000], col = 'green')
plot(draws_motnp2.2$M99_U21[1:1000], type = 'l', ylim = c(0,1))
lines(draws_motnp2.2$M99_U21[1001:2000], col = 'blue')
lines(draws_motnp2.2$M99_U21[2001:3000], col = 'red')
lines(draws_motnp2.2$M99_U21[3001:4000], col = 'green')
plot(draws_motnp2.2$U8_U9[1:1000], type = 'l', ylim = c(0,1))
lines(draws_motnp2.2$U8_U9[1001:2000], col = 'blue')
lines(draws_motnp2.2$U8_U9[2001:3000], col = 'red')
lines(draws_motnp2.2$U8_U9[3001:4000], col = 'green')

at1 <- as.data.frame(adj_tensor[,,1])
colnames(at1) <- c(sort(unique(counts_df$id_1)),'U9')
rownames(at1) <- c(sort(unique(counts_df$id_1)),'U9')
View(at1) # even right back at this stage, there is nothing but 0s in the last section of the matrix

# start over again again... ####
# RE-WRITTEN
adj_tensor <- array(0, c(NROW(unique(counts_df$id_1)), NROW(unique(counts_df$id_2)), NROW(draws_motnp2.2)),
                              dimnames = list(unique(counts_df$id_1), unique(counts_df$id_2), NULL))
for (i in 1:nrow(counts_df)) {
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_motnp2.2[, i+1]
}
View(adj_tensor[, , 1])

at1 <- as.data.frame(adj_tensor[,,1])
colnames(at1) <- sort(unique(counts_df$id_2))
rownames(at1) <- sort(unique(counts_df$id_1))
View(at1)

adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- ((adj_upper - adj_lower)/adj_mid)
adj_range[is.nan(adj_range)] <- 0

ele_links_mid <- as.data.frame(adj_mid)
colnames(ele_links_mid) <- sort(unique(counts_df$id_2))
rownames(ele_links_mid) <- sort(unique(counts_df$id_1))
ele_links_long_m <- pivot_longer(ele_links_mid, cols = everything(), names_to = 'id_2', values_to = 'mid_weight')
ele_links_long_m$id_1 <- rep(sort(unique(counts_df$id_1)), each = 462)

ele_links_rng <- as.data.frame(adj_range)
colnames(ele_links_rng) <- sort(unique(counts_df$id_2))
rownames(ele_links_rng) <- sort(unique(counts_df$id_1))
ele_links_long_r <- pivot_longer(ele_links_rng, cols = everything(), names_to = 'id_2', values_to = 'rng_weight')
ele_links_long_r$id_1 <- rep(sort(unique(counts_df$id_1)), each = 462)

ele_links <- left_join(ele_links_long_m, ele_links_long_r)
ele_links <- ele_links[ele_links$id_2 != ele_links$id_1,]
ele_links$node_1 <- as.integer(as.factor(ele_links$id_1))
ele_links$node_2 <- as.integer(as.factor(ele_links$id_2))
ele_links$dyad_sort <- ifelse(ele_links$node_1 > ele_links$node_2,
                              paste(ele_links$id_2, ele_links$id_1, sep = '_'),
                              paste(ele_links$id_1, ele_links$id_2, sep = '_'))
ele_links$dyad_unsort <- paste(ele_links$id_1, ele_links$id_2, sep = '_')
table(ele_links$dyad_sort)

links <- data.frame(dyad = sort(unique(counts_df$dyad))) %>% 
  separate(col = 'dyad', into = c('id_1','id_2'), sep = '_', remove = F)
links$mid <- NA ; links$rng <- NA
for(i in 1:nrow(links)){
  dyad_row <- ele_links[ele_links$dyad_sort == links$dyad[i],]
  links$mid[i] <- dyad_row$mid_weight[which(dyad_row$dyad_sort == dyad_row$dyad_unsort)]
  links$rng[i] <- dyad_row$rng_weight[which(dyad_row$dyad_sort == dyad_row$dyad_unsort)]
}

weights <- left_join(links, counts_df, by = 'dyad')
min(weights$mid)
which(weights$id_1.x != weights$id_1.y)
which(weights$id_2.x != weights$id_2.y)
weights <- weights[,c(6,1:5,37,38,11:14,19:36)]
colnames(weights)[3:4] <- c('id_1','id_2')

ele_nodes <- read_csv('data_processed/motnp_elenodes_22.01.13.csv')
ele_nodes$sex       <- as.factor(ele_nodes$sex)
ele_nodes$age_class <- as.factor(ele_nodes$age_class)
ele_nodes$dem_class <- as.factor(ele_nodes$dem_class)
ele_nodes <- ele_nodes[,c(12,2:4,5,6,8,9,10,13)]
str(ele_nodes)

nodes_id <- c(sort(unique(ele_nodes$id)))
links_id <- c(sort(unique(ele_links$id_1)),'U9')
which(nodes_id != links_id) # none

freq_links <- read_csv('data_processed/motnp_6thJan2022/motnp_elelinks_22.01.06.csv')
freq_links$age.dyad   <- as.factor(freq_links$age.dyad)
freq_links$sex1       <- as.factor(freq_links$sex1)
freq_links$sex2       <- as.factor(freq_links$sex2)
freq_links$dem_class1 <- as.factor(freq_links$dem_class1)
freq_links$dem_class2 <- as.factor(freq_links$dem_class2)
freq_links <- freq_links[,1:4]
str(freq_links)

mid_links <- weights[,c('id_1','id_2','dem_type','mid')]
str(mid_links)
rng_links <- weights[,c('id_1','id_2','dem_type','rng')]
str(rng_links)

mid_net <- igraph::graph_from_data_frame(d = mid_links,         # data frame of network edges
                                         vertices = ele_nodes,  # node ID and attributes of nodes
                                         directed = F)          # direction of edge strengths
rng_net <- igraph::graph_from_data_frame(d = rng_links,         # data frame of network edges
                                         vertices = ele_nodes,  # node ID and attributes of nodes
                                         directed = F)          # direction of edge strengths

set.seed(1) ; coords <- layout_with_lgl(mid_net)

# test plot
plot(mid_net,
     edge.width = mid_links$mid,
     edge.color = ifelse(mid_links$mid < 0.5,'transparent','black'),
     vertex.size = 5,
     vertex.label = NA,
     layout = coords)
plot(mid_net,
     edge.width = 5*rng_links$rng,
     edge.color = ifelse(mid_links$mid < 0.5,'transparent',rgb(0,0,0,0.25)),
     vertex.color= ifelse(ele_nodes$age_class == 'Adult','seagreen1',
                          ifelse(ele_nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(ele_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 5,
     vertex.label = ele_nodes$sex,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
     vertex.label.cex = 0.5,
     vertex.label.family = 'Helvetica',
     add = T,
     layout = coords)

# plot males
ele_nodes_m <- ele_nodes[ele_nodes$dem_class == 'AM' |
                           ele_nodes$dem_class == 'PM',]
mid_links_m <- mid_links[mid_links$dem_type == 'AM_AM' | 
                           mid_links$dem_type == 'AM_PM' |
                           mid_links$dem_type == 'PM_PM',]
str(mid_links_m)
rng_links_m <- rng_links[rng_links$dem_type == 'AM_AM' | 
                           rng_links$dem_type == 'AM_PM' |
                           rng_links$dem_type == 'PM_PM',]
str(rng_links)

mid_net_m <- igraph::graph_from_data_frame(d = mid_links,
                                           vertices = ele_nodes,
                                           directed = F)
rng_net_m <- igraph::graph_from_data_frame(d = rng_links,
                                           vertices = ele_nodes,
                                           directed = F)
coords_m <- layout_with_lgl(mid_net_m)
plot(mid_net_m,
     edge.width = mid_links_m$mid,
     edge.color = ifelse(mid_links_m$mid < 0.2,'transparent','black'),
     vertex.size = 5,
     vertex.label = NA,
     layout = coords_m)
plot(mid_net_m,
     edge.width = rng_links_m$rng,
     edge.color = ifelse(mid_links_m$mid < 0.2,'transparent',rgb(0,0,0,0.25)),
     vertex.color= ifelse(ele_nodes_m$age_class == 'Adult','seagreen1',
                   ifelse(ele_nodes_m$age_class == 'Pubescent','skyblue',
                   ifelse(ele_nodes_m$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 5,
     vertex.label = ele_nodes_m$sex,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
     vertex.label.cex = 0.5,
     vertex.label.family = 'Helvetica',
     add = T,
     layout = coords_m)

plot(mid_net_m,
     edge.width = mid_links_m$mid,
     edge.color = ifelse(mid_links_m$mid < 0.3,'transparent','black'),
     vertex.size = 5,
     vertex.label = NA,
     layout = coords_m)
plot(mid_net_m,
     edge.width = rng_links_m$rng,
     edge.color = ifelse(mid_links_m$mid < 0.3,'transparent',rgb(0,0,0,0.25)),
     vertex.color= ifelse(ele_nodes_m$age_class == 'Adult','seagreen1',
                          ifelse(ele_nodes_m$age_class == 'Pubescent','skyblue',
                                 ifelse(ele_nodes_m$age_class == 'Juvenile', 'yellow','magenta'))),
     vertex.size = 5,
     vertex.label = ele_nodes_m$sex,
     edge.color = rgb(0.1,0.1,0.1,alpha=0.5),
     vertex.label.cex = 0.5,
     vertex.label.family = 'Helvetica',
     add = T,
     layout = coords_m)

# work out how to plot all from a single igraph object rather than redefining it each time ####
# DOES need the initial column in after all, else it insists on plotting self-links
adj_tensor <- array(0, c(NROW(unique(counts_df$id_1))+1,
                         NROW(unique(counts_df$id_2))+1,
                         NROW(draws_motnp2.2)),
                    dimnames = list(c(unique(counts_df$id_1),'U9'),
                                    c('F1',unique(counts_df$id_2)),
                                    NULL))
for (i in 1:nrow(counts_df)) {
  dyad_row <- counts_df[i, ]
  adj_tensor[dyad_row$id_1, dyad_row$id_2, ] <- draws_motnp2.2[, i+1]
}
at1 <- as.data.frame(adj_tensor[,,1])
View(at1)
adj_quantiles <- apply(adj_tensor, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower <- adj_quantiles[1, , ])
(adj_mid   <- adj_quantiles[2, , ])
(adj_upper <- adj_quantiles[3, , ])
adj_range <- ((adj_upper - adj_lower)/adj_mid) ; adj_range[is.nan(adj_range)] <- 0

# Generate two igraph objects, one from the median and one from the standardised width.
g_mid <- graph_from_adjacency_matrix(adj_mid,   mode="undirected", weighted=TRUE)
g_rng <- graph_from_adjacency_matrix(adj_range, mode="undirected", weighted=TRUE)

# Plot all
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
     vertex.color= ifelse(ele_nodes$age_class == 'Adult','seagreen1',
                          ifelse(ele_nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(ele_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

plot(g_mid,
     edge.width = E(g_mid)$weight,
     vertex.label = NA,
     vertex.size = 5,
     edge.color = ifelse(mid_links$mid < 0.5,'transparent','black'),
     layout = coords)
plot(g_mid,
     edge.width = E(g_rng)$weight,
     edge.color = ifelse(mid_links$mid < 0.5,'transparent',rgb(0,0,0,0.25)),
     vertex.size = 5,
     vertex.label = c(sort(unique(counts_df$id_1)),'U9'),
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(ele_nodes$age_class == 'Adult','seagreen1',
                          ifelse(ele_nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(ele_nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords, add = TRUE)

#at1$sum <- rowSums(at1)
#at1$degree_0.1 <- NA
#at1$degree_0.2 <- NA
#at1$degree_0.3 <- NA
#at1$degree_0.4 <- NA
#at1$degree_0.5 <- NA
#for(i in 1:nrow(at1)){
#  at1$degree_0.1[i] <- length(which(at1[i,]>0.1))
#  at1$degree_0.2[i] <- length(which(at1[i,]>0.2))
#  at1$degree_0.3[i] <- length(which(at1[i,]>0.3))
#  at1$degree_0.4[i] <- length(which(at1[i,]>0.4))
#  at1$degree_0.5[i] <- length(which(at1[i,]>0.5))
#}
#summary(at1$degree_0.5) # minimum 4 -- doesn't seem right. U9 sh

ele_nodes$degree_0.1 <- NA
ele_nodes$degree_0.2 <- NA
ele_nodes$degree_0.3 <- NA
ele_nodes$degree_0.4 <- NA
ele_nodes$degree_0.5 <- NA
#tapply(X = links$mid, INDEX = links$id_1, FUN = sum)
count <- function(x, a){
  length(which(x > a))
}
count(test$mid, a = 0.1)
for(i in 1:NROW(ele_nodes)){
  ele_rows <- links[links$id_1 == ele_nodes$id[i] | links$id_2 == ele_nodes$id[i],]
  ele_nodes$degree_0.1[i] <- count(ele_rows$mid, a = 0.1)
  ele_nodes$degree_0.2[i] <- count(ele_rows$mid, a = 0.2)
  ele_nodes$degree_0.3[i] <- count(ele_rows$mid, a = 0.3)
  ele_nodes$degree_0.4[i] <- count(ele_rows$mid, a = 0.4)
  ele_nodes$degree_0.5[i] <- count(ele_rows$mid, a = 0.5)
}
test <- links[links$id_1 == 'F1',] ; sum(test$mid)
length(which(test$mid > 0.1)) # YAY!
table(test$mid) # actually counting = also 14

#g_mid_0.5 <- delete.vertices(g_mid, v = which(ele_nodes$degree_0.5 == 0)) -- this doesn't work, because ele_nodes order = F1, F2, F3... and g_mid order = F1, F10, F100...
nodes <- data.frame(id = sort(ele_nodes$id))      # reorder so that in same order
nodes <- left_join(nodes, ele_nodes, by = 'id')
g_mid_0.5 <- delete.vertices(g_mid, v = which(nodes$degree_0.5 == 0))
g_mid_0.5

coords_0.5 <- layout_nicely(g_mid_0.5)
plot(g_mid_0.5,
     edge.width = 3*E(g_mid_0.5)$weight,
     vertex.label = NA,
     vertex.shape = ifelse(nodes$degree_0.5 == 0, 'none', 'circle'),
     vertex.size = 5,
     edge.color = ifelse(mid_links$mid < 0.5,'transparent','black'),
     layout = coords_0.5)
plot(g_mid_0.5,
     edge.width = E(g_rng)$weight,
     edge.color = ifelse(mid_links$mid < 0.5,'transparent',rgb(0,0,0,0.25)),
     vertex.shape = ifelse(nodes$degree_0.5 == 0, 'none', 'circle'),
     vertex.size = 5,
     vertex.label = nodes$id[which(nodes$degree_0.5 > 0)],
     vertex.label.dist = 0,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.color= ifelse(nodes$age_class == 'Adult','seagreen1',
                          ifelse(nodes$age_class == 'Pubescent','skyblue',
                                 ifelse(nodes$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords_0.5, add = TRUE)


test <- delete.vertices(graph = g_mid, v = 'F1')
test
plot(test)

test <- delete.vertices(graph = g_mid, v = 2)
test
plot(test)

which(ele_nodes$degree_0.2 == 0) # 154 155 159 162 165 168 170 180 191 192 193 194 201 206 208 216 234 287
ele_nodes$id[which(ele_nodes$degree_0.2 == 0)] # "M3"   "M4"   "M8"   "M11"  "M15"  "M18"  "M20"  "M31"  "M42"  "M43"  "M44"  "M45"  "M52" "M57"  "M59"  "M67"  "M85"  "M141"

test <- delete.vertices(graph = g_mid, v = c("M3","M4","M8","M11","M15","M18","M20","M31","M42","M43","M44","M45","M52","M57","M59","M67","M85","M141"))
plot(test)

test <- delete.vertices(graph = g_mid, v = ele_nodes$id[which(ele_nodes$degree_0.2 == 0)])
plot(test)

test <- delete.vertices(graph = g_mid, v = ele_nodes$id[which(ele_nodes$degree_0.5 == 0)])
plot(test)
test2 <- delete.vertices(graph = g_rng, v = ele_nodes$id[which(ele_nodes$degree_0.5 == 0)])
plot(test2)

coords_0.5 <- layout_nicely(test)
plot(test,
     edge.width = ifelse(E(test)$weight > 0.5, E(test)$weight, 0),
     edge.color = 'black',
     vertex.size = 5,
     vertex.label = NA,
     layout = coords_0.5)
plot(test,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(test)$weight > 0.5, E(test2)$weight, 0),
     vertex.size = 5,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Adult', 'seagreen1',
                           ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Pubescent','skyblue',
                                  ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords_0.5, add = TRUE)

# plot 0.3
g_mid_0.3 <- delete.vertices(graph = g_mid,
                             v = ele_nodes$id[which(ele_nodes$degree_0.5 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng,
                             v = ele_nodes$id[which(ele_nodes$degree_0.5 == 0)])
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.5, E(g_mid_0.3)$weight, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.5, E(g_rng_0.3)$weight, 0),
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Adult',
                           'seagreen1',
                    ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Pubescent',
                           'skyblue',
                    ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Juvenile',
                           'yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

coords_0.3 <- layout_in_circle(g_mid_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.5, E(g_mid_0.3)$weight, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.5, E(g_rng_0.3)$weight, 0),
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Adult',
                           'seagreen1',
                           ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Pubescent',
                                  'skyblue',
                                  ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Juvenile',
                                         'yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

coords <- layout_in_circle(g_mid)
plot(g_mid,
     edge.width = ifelse(E(g_mid)$weight < 0.05, E(g_mid)$weight, 0),
     edge.color = 'black',
     vertex.size = 3,
     vertex.label = NA,
     layout = coords)
plot(g_mid,
     edge.color = rgb(0,0,0,0.05),
     edge.width = ifelse(E(g_mid)$weight < 0.05, E(g_rng)$weight, 0), # OK this is definitely working -- they looked to not be adding anything, but actually they're just very narrow
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Adult',
                           'seagreen1',
                    ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Pubescent',
                           'skyblue',
                    ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Juvenile',
                           'yellow','magenta'))),
     layout = coords, add = TRUE)

coords <- layout_in_circle(g_mid)
plot(g_mid,
     edge.width = ifelse(E(g_mid)$weight > 0.5, E(g_mid)$weight, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label = NA,
     layout = coords)
plot(g_mid,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid)$weight > 0.5, E(g_rng)$weight, 0),
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Adult',
                           'seagreen1',
                           ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Pubescent',
                                  'skyblue',
                                  ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Juvenile',
                                         'yellow','magenta'))),
     layout = coords, add = TRUE)



g_mid_0.5 <- delete.vertices(graph = g_mid, v = ele_nodes$id[which(ele_nodes$degree_0.5 == 0)])
plot(g_mid_0.5)
g_rng_0.5 <- delete.vertices(graph = g_rng, v = ele_nodes$id[which(ele_nodes$degree_0.5 == 0)])
plot(g_rng_0.5)

coords_0.5 <- layout_nicely(g_mid_0.5)
plot(g_mid_0.5,
     edge.width = ifelse(E(g_mid_0.5)$weight > 0.5, E(g_mid_0.5)$weight, 0),
     edge.color = 'black',
     vertex.size = 5,
     vertex.label = NA,
     layout = coords_0.5)
plot(g_mid_0.5,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.5)$weight > 0.5, E(g_rng_0.5)$weight, 0),
     vertex.size = 5,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Adult', 'seagreen1',
                           ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Pubescent','skyblue',
                                  ifelse(ele_nodes[which(ele_nodes$degree_0.5 != 0),]$age_class == 'Juvenile', 'yellow','magenta'))),
     layout = coords_0.5, add = TRUE)

# plot 0.3
g_mid_0.3 <- delete.vertices(graph = g_mid,
                             v = ele_nodes$id[which(ele_nodes$degree_0.3 == 0)])
g_rng_0.3 <- delete.vertices(graph = g_rng,
                             v = ele_nodes$id[which(ele_nodes$degree_0.3 == 0)])
coords_0.3 <- layout_nicely(g_mid_0.3)
plot(g_mid_0.3,
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_mid_0.3)$weight, 0),
     edge.color = 'black',
     vertex.size = 7,
     vertex.label = NA,
     layout = coords_0.3)
plot(g_mid_0.3,
     edge.color = rgb(0,0,0,0.25),
     edge.width = ifelse(E(g_mid_0.3)$weight > 0.3, E(g_rng_0.3)$weight, 0),
     vertex.size = 7,
     vertex.label.color = 'black',
     vertex.label.family = 'Helvetica',
     vertex.label.cex = 0.5,
     vertex.label.dist = 0,
     vertex.color = ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Adult',
                           'seagreen1',
                    ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Pubescent',
                           'skyblue',
                    ifelse(ele_nodes[which(ele_nodes$degree_0.3 != 0),]$age_class == 'Juvenile',
                           'yellow','magenta'))),
     layout = coords_0.3, add = TRUE)

#plot males
g_mid_m <- delete.vertices(graph = g_mid,
                           v = ele_nodes$id[which(ele_nodes$dem_class != 'AM' &
                                                    ele_nodes$dem_class != 'PM')])
g_rng_m <- delete.vertices(graph = g_rng,
                             v = ele_nodes$id[which(ele_nodes$dem_class != 'AM' &
                                                      ele_nodes$dem_class != 'PM')])
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
     vertex.color = ifelse(ele_nodes[which(ele_nodes$dem_class == 'AM' | 
                                             ele_nodes$dem_class == 'PM'),]$age_class == 'Adult',
                           'seagreen1', 'skyblue'),
     layout = coords_m, add = TRUE)

# plot males >0.3 -- OK I THINK this is working, but it keeps only showing a handful of the labels. If I then change the size of the plot window, they normally appear so I think it's just that my laptop is struggling, but there might be a problem with these plots still.
g_mid_m <- delete.vertices(graph = g_mid,
                           v = ele_nodes$id[which(ele_nodes$dem_class != 'AM' &
                                                    ele_nodes$dem_class != 'PM' |
                                                    ele_nodes$degree_0.3 == 0)])
g_rng_m <- delete.vertices(graph = g_rng,
                           v = ele_nodes$id[which(ele_nodes$dem_class != 'AM' &
                                                    ele_nodes$dem_class != 'PM' |
                                                    ele_nodes$degree_0.3 == 0)])

#anti_m0.3 <- ele_nodes[which(ele_nodes$dem_class != 'AM' & ele_nodes$dem_class != 'PM' | ele_nodes$degree_0.3 == 0),]
#View(anti_m0.3)

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
     vertex.color = ifelse(ele_nodes[which(ele_nodes$dem_class == 'AM' | 
                                             ele_nodes$dem_class == 'PM'),]$age_class == 'Adult',
                           'seagreen1', 'skyblue'),
     layout = coords_m, add = TRUE)



###
################ 9) Just playing with gregariousness ideas a little ################
head(means_motnp2.2)
length(means_motnp2.2$dyad) # this is all of the dyads
edges <- left_join(x = means_motnp2.2, y = counts_df, by = 'dyad')
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

ggplot(elephants, aes(x = age_cat_id, y = mean_edge_apf, colour = sex))+
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
       y = 'median edge weight')+
  theme_light()
