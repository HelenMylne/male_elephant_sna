## Bayesian analysis of ALERT data
#### Set up ####
library(tidyverse)
library(dplyr)
library(rstan)
library(rethinking)
library(igraph)
library(dagitty)
library(cmdstanr)

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

counts_df$age_cat_id_1 <- as.numeric(counts_df$age_cat_id_1)
counts_df$age_cat_id_2 <- as.numeric(counts_df$age_cat_id_2)

counts_df$age_diff <- abs(counts_df$age_cat_id_1 - counts_df$age_cat_id_2)
summary(counts_df$age_diff) # 0-6

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
# X_ij ∼ Binomial(D_ij, p_ij)
# logit(p_ij) = w_ij

# weight ~ Beta(a + T_ij, b + A_ij) --> T = sightings Together, A = sightings Apart

################ 5) Define priors and hyperpriors, run prior predictive checks ################
##### run without a and b priors #####
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

### Compile Stan model
set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')
mod <- cmdstan_model("simpleBetaNet_DWF_22.01.23.stan")

### Fit model
edge_weight <- mod$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4
)

### check model
edge_weight$summary()
edge_weight$output_files()
output1 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpwMJ25r/simpleBetaNet_DWF_22.01.23-202201282116-1-251293.csv")
output2 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpwMJ25r/simpleBetaNet_DWF_22.01.23-202201282116-2-251293.csv")
output3 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpwMJ25r/simpleBetaNet_DWF_22.01.23-202201282116-3-251293.csv")
output4 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpwMJ25r/simpleBetaNet_DWF_22.01.23-202201282116-4-251293.csv")
draws1 <- as.data.frame(output1$post_warmup_draws)
draws2 <- as.data.frame(output1$post_warmup_draws)
draws3 <- as.data.frame(output1$post_warmup_draws)
draws4 <- as.data.frame(output1$post_warmup_draws)

# all chains give perfectly identical outputs -- am I understanding the output wrong and these aren't the draws?
mean(draws1$`1.weight[1]`) ; mean(draws2$`1.weight[1]`) ; mean(draws3$`1.weight[1]`) ; mean(draws4$`1.weight[1]`)
mean(draws1$`1.weight[2]`) ; mean(draws2$`1.weight[2]`) ; mean(draws3$`1.weight[2]`) ; mean(draws4$`1.weight[2]`)
mean(draws1$`1.weight[3]`) ; mean(draws2$`1.weight[3]`) ; mean(draws3$`1.weight[3]`) ; mean(draws4$`1.weight[3]`)

# build traceplots -- all look pretty good, overlap nicely and related/within-group elephants have a higher association index
plot(draws1$`1.weight[1]`, type = 'l',
     ylim = c(0,1), las = 1, ylab = 'edge weight')
lines(draws1$`1.weight[2]`, col = 'tan')
lines(draws1$`1.weight[3]`, col = 'orange')
lines(draws1$`1.weight[4]`, col = 'green')
lines(draws1$`1.weight[5]`, col = 'chocolate')
lines(draws1$`1.weight[6]`, col = 'blue')
lines(draws1$`1.weight[7]`, col = 'red')
lines(draws1$`1.weight[8]`, col = 'seagreen')
lines(draws1$`1.weight[9]`, col = 'purple')
lines(draws1$`1.weight[10]`,col = 'magenta')
lines(draws1[which(dyads$family_1 == dyads$id_2)[1]+1], col = 'black')       # traceplot of mother-calf relationships
lines(draws1[which(dyads$family_1 == dyads$id_2)[2]+1], col = 'tan')         # +1 because first column is not linked to a dyad
lines(draws1[which(dyads$family_1 == dyads$id_2)[3]+1], col = 'orange')
lines(draws1[which(dyads$family_1 == dyads$id_2)[4]+1], col = 'green')
lines(draws1[which(dyads$family_1 == dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1[which(dyads$family_1 == dyads$id_2)[6]+1], col = 'blue')
lines(draws1[which(dyads$family_1 == dyads$id_2)[7]+1], col = 'red')
lines(draws1[which(dyads$family_1 == dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1[which(dyads$family_1 == dyads$id_2)[9]+1], col = 'purple')
lines(draws1[which(dyads$family_1 == dyads$id_2)[10]+1], col = 'magenta')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[1]+1], col = 'black')      # traceplots of elephants in the same group but which are not mother-calf pairs
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[2]+1], col = 'tan')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[3]+1], col = 'orange')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[4]+1], col = 'green')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[5]+1], col = 'chocolate')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[6]+1], col = 'blue')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[7]+1], col = 'red')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[8]+1], col = 'seagreen')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[9]+1], col = 'purple')
lines(draws1[which(dyads$clique_1 == dyads$clique_2 & dyads$family_1 != dyads$id_2)[10]+1], col = 'magenta')

# plot each chain individually to check mixing -- every chain has identical outputs, with a perfect match of the draws made -- this doesn't seem right...
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

# mean and credible interval for each dyad weight
means <- data.frame(dyad = dyads$dyad, mean = colMeans(draws1[2:7141])) # these are matching up to the right dyads because the traceplots above use the row numbers to specify which ones to plot and they come out at the right approximate edge weight
draws1_pi <- apply(draws1[2:7141], 2, rethinking::PI)
means$pi_lwr <- draws1_pi[1,]
means$pi_upr <- draws1_pi[2,]

# boxplot comparing types of dyads
plot_data <- left_join(x = means, y = dyads, by = 'dyad')

boxplot(plot_data$mean ~ plot_data$pair_type,
        las = 1, xlab = 'dyad relationship', ylab = 'mean edge weight', ylim = c(0,1),
        col = c('turquoise1','orchid','firebrick1')) 
abline(h = mean(dyads[dyads$pair_type ==   'family',]$edge), lty = 2, col = 'blue')     # plot true mean for the group
abline(h = mean(dyads[dyads$pair_type ==    'group',]$edge), lty = 2, col = 'purple3')  # plot true mean for the group
abline(h = mean(dyads[dyads$pair_type == 'no_group',]$edge), lty = 2, col = 'darkred')  # plot true mean for the group
points(x = rnorm(nrow(dyads[dyads$pair_type == 'family',]), 1,0.15),y = dyads[dyads$pair_type == 'family',]$edge,
       pch = 4, cex = 0.8, col = 'blue')                      # true values
points(x = rnorm(nrow(dyads[dyads$pair_type == 'group',]),  2,0.15), y = dyads[dyads$pair_type == 'group',]$edge,
       pch = 4, cex = 0.8, col = 'purple3')                   # true values
points(x = rnorm(nrow(dyads[dyads$pair_type == 'no_group',]), 3,0.15), y = dyads[dyads$pair_type == 'no_group',]$edge,
       pch = 4, cex = 0.8, col = col.alpha('darkred', 0.2))   # true values

ggplot(data = plot_data, aes(x = pair_type, y = mean, fill = pair_type))+
  geom_boxplot()+
  geom_point(aes(x = pair_type, y = edge, fill = pair_type),
             pch = 4, colour = col.alpha('red', 0.2))+
  theme_classic()+
  scale_fill_viridis_d()

# plot output: weight ~ age_difference
plot_data$age_diff_jitter <- jitter(plot_data$age_diff,1) # make plot slightly more visible
plot(mean ~ age_diff_jitter, data = plot_data,
     xlab = 'age category difference', ylab = 'edge weight',
     ylim = c(0,1), pch = 19, col = col.alpha(rangi2, 0.2), las = 1)
#for(i in 1:nrow(plot_data)){                             # adding error bars actually just makes it impossible to see anything
#  lines(x = c(plot_data$age_diff_jitter[i], plot_data$age_diff_jitter[i]),
#        y = c(plot_data$pi_lwr[i]  , plot_data$pi_upr[i]  ),
#        lty = ifelse(plot_data$sex_diff[i] == 0, 1, 2),
#        col = ifelse(plot_data$sex_dyad[i] == 'M_M' | plot_data$sex_dyad[i] == 'F_M' | 
#                       plot_data$sex_dyad[i] == 'M_U', 'red', 'green'))
#}

##### add a and b priors #####
# Compile Stan model
mod2 <- cmdstan_model("simpleBetaNet_HKM_22.01.25.stan")

# create data list
simdat_ls <- list(
  n_dyads  = nrow(dyads),        # Number of dyads
  together = dyads$event_count,  # Number of sightings of each dyad together
  apart    = dyads$apart         # Number of sightings of each dyad apart
)
summary(simdat_ls$together)      # definitely nothing in this that could be classed as a non-negative integer

# Fit model
edge_weight_ab <- mod2$sample(
  data = simdat_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

a <- rep(min(rnorm(1000, 10, 0.01)),1000) ; min(a) # smallest this ever reaches is >9.9 so cannot be creating a negative number when added to together, even when together is 0

################ 8) Run model on real standardised data -- Binomial model to calculate SRI with uncertainty ################
### create data list -- can contain no NA values in any column, even if column is not specified in model
counts_ls <- list(
  n_dyads  = nrow(counts_df),          # total number of times one or other of the dyad was observed
  together = counts_df$all_events,     # count number of sightings seen together
  apart    = counts_df$apart)          # count number of sightings seen apart

### run model
mod

### Fit model (slow)
weight_motnp <- mod$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4, 
  parallel_chains = 4)

### check model
weight_motnp
weight_motnp$summary()
weight_motnp$output_files()
output1 <- read_cmdstan_csv("")
output2 <- read_cmdstan_csv("")
output3 <- read_cmdstan_csv("")
output4 <- read_cmdstan_csv("")

draws1 <- as.data.frame(output1$post_warmup_draws)
draws2 <- as.data.frame(output2$post_warmup_draws)
draws3 <- as.data.frame(output3$post_warmup_draws)
draws4 <- as.data.frame(output4$post_warmup_draws)

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
# OK so every chain has identical outputs, with a perfect match of the draws made -- this doesn't seem right... -- thats because I was a plank and assigned everything as output1 rather than outputs1-4

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

################ 9) Plot results so far to send to ALERT ################
### create data list
counts_ls <- list(
  n_dyads  = nrow(counts_df),          # total number of times one or other of the dyad was observed
  together = counts_df$all_events,     # count number of sightings seen together
  apart    = counts_df$apart)          # count number of sightings seen apart

### run model
mod

### Fit model (slow)
weight_motnp <- mod$sample(
  data = counts_ls, 
  seed = 12345, 
  chains = 4,
  parallel_chains = 4)

### check model
weight_motnp$summary()
weight_motnp$output_files()
output1 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpGT4xWN/simpleBetaNet_DWF_22.01.23-202202011550-1-53e232.csv")
output2 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpGT4xWN/simpleBetaNet_DWF_22.01.23-202202011550-2-53e232.csv")
output3 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpGT4xWN/simpleBetaNet_DWF_22.01.23-202202011550-3-53e232.csv")
output4 <- read_cmdstan_csv("/var/folders/yt/j5ms6jn178n94npylplm3qzr0000gn/T/RtmpGT4xWN/simpleBetaNet_DWF_22.01.23-202202011550-4-53e232.csv")

draws1 <- as.data.frame(output1$post_warmup_draws)
draws2 <- as.data.frame(output2$post_warmup_draws)
draws3 <- as.data.frame(output3$post_warmup_draws)
draws4 <- as.data.frame(output4$post_warmup_draws)
plot(draws1$`1.weight[1]`, pch = 19)
points(draws2$`1.weight[1]`, pch = 19)
points(draws3$`1.weight[1]`, pch = 19)
points(draws4$`1.weight[1]`, pch = 19)

plot(draws1$`1.weight[1]`, type = 'l', ylim = c(0,0.25), las = 1)
lines(draws2$`1.weight[1]`, col = 'red')
lines(draws3$`1.weight[1]`, col = 'green')
lines(draws4$`1.weight[1]`, col = 'blue')

plot(draws1$`1.weight[2]`, type = 'l', ylim = c(0,0.25), las = 1)
lines(draws2$`1.weight[2]`, col = 'red')
lines(draws3$`1.weight[2]`, col = 'green')
lines(draws4$`1.weight[2]`, col = 'blue')

plot(draws1$`1.weight[3]`, type = 'l', ylim = c(0,0.25), las = 1)
lines(draws2$`1.weight[3]`, col = 'red')
lines(draws3$`1.weight[3]`, col = 'green')
lines(draws4$`1.weight[3]`, col = 'blue')

weight_motnp$sampler_diagnostics()
weight_motnp$cmdstan_diagnose()
draws <- rbind(draws1, draws2, draws3, draws4)
rm(output1, output2, output3, output4)
colnames(draws)[2:106954] <- counts_df$dyad
tidy_draws <- pivot_longer(draws[,2:106954], cols = everything(), names_to = 'dyad', values_to = 'draw')
tidy_draws$chain <- rep(1:4, each = 106953000)
tidy_draws$index <- rep(rep(1:1000, each = 106953),4)
head(tidy_draws, 10)
tail(tidy_draws, 10)

# plot for Sierra (F52), herd members = F60+U21+F98, calf = U17
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

# density plots
dens(tidy_sierra[tidy_sierra$id == 'M40',]$draw, xlab = 'edge weight estimation',
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
    ggplot(data = tidy_sierra[tidy_sierra$chain == 1], aes(x = draw, fill = ID))+
    geom_density()+  # why different values to the rethinking version?
    scale_fill_viridis_d(alpha = 0.4)+
    scale_x_continuous('edge weight estimation', expand = c(0,0), limits = c(0,1.1))+
    scale_y_continuous(expand = c(0,0), limits = c(0,55))+
    theme_classic()+
    theme(legend.position = c(0.5, 0.6)))
ggsave(path = 'outputs/motnp_22.01.31/', filename = 'motnp_sierra_density.pdf',
       plot = sierra_density,
       width = 20, height = 24, units = 'cm')

(sierra_density <-
    ggplot(data = tidy_sierra, aes(x = draw, col = as.factor(chain)))+
    geom_density()+  # why different values to the rethinking version?
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

# boxplot of all edge types
means <- data.frame(dyad = counts_df$dyad, mean = colMeans(draws[2:106954]),
                    median = apply(draws[2:106954], 2, FUN = median))
draws_pi <- apply(draws[2:106954], 2, rethinking::PI)
means$pi_lwr <- draws_pi[1,]
means$pi_upr <- draws_pi[2,]
plot_data <- left_join(x = means, y = counts_df, by = 'dyad')
head(plot_data)
plot_data$dem_class_1_cat <- as.integer(as.factor(plot_data$dem_class_1))
plot_data$dem_class_2_cat <- as.integer(as.factor(plot_data$dem_class_2))
plot_data$dem_type_short <- ifelse(plot_data$dem_class_1_cat <= plot_data$dem_class_2_cat,
                                   paste(plot_data$dem_class_1, plot_data$dem_class_2, sep = '_'),
                                   paste(plot_data$dem_class_2, plot_data$dem_class_1, sep = '_'))
sort(unique(plot_data$dem_type_short))
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
plot_data <- left_join(plot_data, types, by = 'dem_type_short')
plot_data$join <- paste(plot_data$join, ' ', sep = ' ')
which(is.na(plot_data$join))

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

# only adults/pubescents
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

#### create network plots ####
head(means)
head(plot_data)
plot_data$ci_range <- plot_data$pi_upr - plot_data$pi_lwr

# Hart example: "A more useful format for network data is usually adjacency matrices, rather than edge lists, so instead we’ll convert the distribution of edge lists to a distribution of adjacency matrices, and store the result in an 8 x 8 x 4000 tensor, as there are 8 nodes and 4000 samples from the posterior." -- I need a 463*463*1000 tensor
length(unique(plot_data$id_1))+1 # number of individuals

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
dyads <- colnames(draws)
sierra_dyads
which(dyads == "F52_F60" | dyads == "F52_F8" | dyads == "F52_F98" | dyads == "F52_M15" | dyads == "F52_M26" | dyads == "F52_M40" |
        dyads == "F52_U17" | dyads == "F52_U21" | dyads == "F60_F8"  | dyads == "F60_F98" | dyads == "F60_M15" | dyads == "F60_M26" |
        dyads == "F60_M40" | dyads == "F60_U17" | dyads == "F60_U21" | dyads == "F8_F98"  | dyads == "F8_M15"  | dyads == "F8_M26"  |
        dyads == "F8_M40"  | dyads == "F8_U17"  | dyads == "F8_U21"  | dyads == "F98_M15" | dyads == "F98_M26" | dyads == "F98_M40" | 
        dyads == "F98_U17" | dyads == "F98_U21" | dyads == "M15_M26" | dyads == "M15_M40" | dyads == "M15_U17" | dyads == "M15_U21" |
        dyads == "M26_M40" | dyads == "M26_U17" | dyads == "M26_U21" | dyads == "M40_U17" | dyads == "M40_U21" | dyads == "U17_U21")
draws_test <- draws[,c(40897,40918,40938,40992,41104,41120,41193,41198,44140,44160,44214,44326,44342,44415,44420,51363,51417,51529,
                       51545,51618,51623,57867,57979,57995,58068,58073,73396,73412,73485,73490,96092,96165,96170,98381,98386,105248)]
head(counts_df_test)
# logit_p_samples <- extract(weight_motnp)$logit_p
adj_tensor_test <- array(0, c(9, 9, 4000))
counts_df_test$dyad_id2 <- as.integer(as.factor(counts_df_test$dyad_id))
counts_df_test$node_1_2 <- as.integer(as.factor(counts_df_test$node_1))
counts_df_test$node_2_2 <- as.integer(as.factor(counts_df_test$node_2))+1
for(dyad_id in 1:nrow(counts_df_test)) {
  dyad_row <- counts_df_test[counts_df_test$dyad_id2 == dyad_id, ]
  adj_tensor_test[dyad_row$node_1_2, dyad_row$node_2_2, ] <- draws_test[, dyad_id]
}
adj_tensor_test[, , 1] # Print the first sample of the posterior distribution over adjacency matrices
adj_tensor_test[, , 2]
adj_tensor_test[, , 3]
adj_tensor_test[, , 4]
adj_tensor_test[, , 5]
adj_tensor_test[, 1, ]
adj_tensor_test[, 2, ]
adj_tensor_test[, 3, ]
adj_tensor_test[, 4, ]
adj_tensor_test[, 5, ]
adj_tensor_test[1, , ]
adj_tensor_test[2, , ]
adj_tensor_test[3, , ]
adj_tensor_test[4, , ]
adj_tensor_test[5, , ]

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles <- apply(adj_tensor_test, c(1, 2), function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
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
plot(g_mid, edge.width = 3*E(g_mid)$weight, edge.color = "black",  layout = coords)
plot(g_mid, edge.width = 3*E(g_rng)$weight, edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = c(sort(unique(counts_df_test$id_1)),'U21'),
     vertex.label.dist = 0, vertex.label.color="black", layout=coords, add=TRUE)

check <- read_csv('data_processed/motnp_6thJan2022/motnp_id_22.01.06.csv')
check[check$id_no == 'F0008',]$family # M26
check[check$id_no == 'F0008',]$herd   # B3
check[check$id_no == 'M0015',]$herd   # B3

plot(g_mid, edge.width = E(g_mid)$weight, edge.color = "black",  layout = coords)
plot(g_mid, edge.width = E(g_rng)$weight, edge.color = rgb(0, 0, 0, 0.25), 
     vertex.label = c(sort(unique(counts_df_test$id_1)),'U21'),
     vertex.label.dist = 0, vertex.label.color="black", layout=coords, add=TRUE)
unique(counts_df_test$id_1)

# work out how to filter data for plotting
apply(draws_test, 2, median)
# F52_F60 = 0.863223000      F52_F8 = 0.033248300    F52_F98 = 0.623663500    F52_M15 = 0.009397185    F52_M26 = 0.034015850    F52_M40 = 0.054451600    F52_U17 = 0.935521500     F52_U21 = 0.837885000    F60_F8 = 0.035207550     F60_F98 = 0.592663000    F60_M15 = 0.009502835    F60_M26 = 0.036595750     F60_M40 = 0.057125050    F60_U17 = 0.820356000    F60_U21 = 0.954352000     F8_F98 = 0.015987200     F8_M15 = 0.036715550      F8_M26 = 0.864075500     F8_M40 = 0.016578300     F8_U17 = 0.033878750     F8_U21 = 0.036029300     F98_M15 = 0.025888000     F98_M26 = 0.015917100    F98_M40 = 0.027839400    F98_U17 = 0.589748500    F98_U21 = 0.567639500    M15_M26 = 0.038345600     M15_M40 = 0.077398150    M15_U17 = 0.009718925    M15_U21 = 0.010018430    M26_M40 = 0.016436100    M26_U17 = 0.033737150     M26_U21 = 0.036970250    M40_U17 = 0.055941700    M40_U21 = 0.058043900    U17_U21 = 0.795998000 
adj_mid
#      [,1]     [,2]       [,3]      [,4]        [,5]       [,6]       [,7]        [,8]       [,9]
# [1,]    0 0.863223 0.03324830 0.6236635 0.009397185 0.03401585 0.05445160 0.935521500 0.83788500
# [2,]    0 0.000000 0.03520755 0.5926630 0.009502835 0.03659575 0.05712505 0.820356000 0.95435200
# [3,]    0 0.000000 0.00000000 0.0159872 0.036715550 0.86407550 0.01657830 0.033878750 0.03602930
# [4,]    0 0.000000 0.00000000 0.0000000 0.025888000 0.01591710 0.02783940 0.589748500 0.56763950
# [5,]    0 0.000000 0.00000000 0.0000000 0.000000000 0.03834560 0.07739815 0.009718925 0.01001843
# [6,]    0 0.000000 0.00000000 0.0000000 0.000000000 0.00000000 0.01643610 0.033737150 0.03697025
# [7,]    0 0.000000 0.00000000 0.0000000 0.000000000 0.00000000 0.00000000 0.055941700 0.05804390
# [8,]    0 0.000000 0.00000000 0.0000000 0.000000000 0.00000000 0.00000000 0.000000000 0.79599800
# [9,]    0 0.000000 0.00000000 0.0000000 0.000000000 0.00000000 0.00000000 0.000000000 0.00000000

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

#plot(g_mid, edge.width = E(g_mid)$weight, edge.color = "black",  layout = coords)
#plot(g_mid, edge.width = E(g_rng)$weight, edge.color = rgb(0, 0, 0, 0.25), 
#     vertex.label = c(unique(counts_df$id_1),'U9'),
#     vertex.color= ifelse(counts_df$age_class_1 == 'Adult','seagreen1',
#                          ifelse(counts_df$age_class_1 == 'Pubescent','skyblue',
#                                 ifelse(counts_df$age_class_1 == 'Juvenile', 'yellow',
#                                        'magenta'))),
#     vertex.label.dist = 0, vertex.label.color="black", layout=coords, add=TRUE)

### plot males only  -- still too many to see anything! ####
counts_df_male <- counts_df[counts_df$dem_class_1 == 'AM' | counts_df$dem_class_1 == 'PM',]
counts_df_male <- counts_df_male[counts_df_male$dem_class_2 == 'AM' | counts_df_male$dem_class_2 == 'PM',]
male_dyads <- counts_df_male$dyad
dyads <- colnames(draws)
male_dyads  # 23653 dyads -- can't subset by hand!
draws_male <- draws[,male_dyads]
colnames(draws_male)

adj_tensor_male <- array(0, c(length(unique(counts_df_male$id_1))+1,
                              length(unique(counts_df_male$id_1))+1,
                              length(draws_male$M1_M10)))
counts_df_male$dyad_id2 <- as.integer(as.factor(counts_df_male$dyad_id))
counts_df_male$node_1_nogaps <- as.integer(as.factor(counts_df_male$node_1))
counts_df_male$node_2_nogaps <- as.integer(as.factor(counts_df_male$node_2))+1
for (dyad_id in 1:nrow(counts_df_male)) {
  dyad_row <- counts_df_male[counts_df_male$dyad_id2 == dyad_id, ]
  adj_tensor_male[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws[, dyad_id+1]
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


plot(g_mid_m, edge.width = E(g_mid_m)$weight, layout = coords,
     edge.color = 'black',
     vertex.color = 'transparent')
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


### plot ADULT males only  -- still too many to see anything! ####
counts_df_am <- counts_df[counts_df$dem_class_1 == 'AM' & counts_df$dem_class_2 == 'AM',]
am_dyads <- counts_df_am$dyad
am_dyads  # 23653 dyads -- can't subset by hand!
draws_am <- draws[,am_dyads]
colnames(draws_am)

adj_tensor_am <- array(0, c(length(unique(counts_df_am$id_1))+1,
                            length(unique(counts_df_am$id_1))+1,
                            length(draws_am$M1_M10)))
counts_df_am$dyad_id2 <- as.integer(as.factor(counts_df_am$dyad_id))
counts_df_am$node_1_nogaps <- as.integer(as.factor(counts_df_am$node_1))
counts_df_am$node_2_nogaps <- as.integer(as.factor(counts_df_am$node_2))+1
for (dyad_id in 1:nrow(counts_df_am)) {
  dyad_row <- counts_df_am[counts_df_am$dyad_id2 == dyad_id, ]
  adj_tensor_am[dyad_row$node_1_nogaps, dyad_row$node_2_nogaps, ] <- draws[, dyad_id+1]
}
adj_tensor_am[, , 1]

# Calculate lower, median, and upper quantiles of edge weights. Lower and upper give credible intervals.
adj_quantiles_am <- apply(adj_tensor_am, c(1, 2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
(adj_lower_am <- adj_quantiles_am[1, , ])
(adj_mid_am   <- adj_quantiles_am[2, , ])
(adj_upper_am <- adj_quantiles_am[3, , ])

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
#0.01020 0.03851 0.06417 0.08625 0.11038 0.97140
plot(g_mid_am, edge.width = E(g_mid_am)$weight,
     edge.color = ifelse(E(g_mid_am)$weight < 0.11, 'transparent','black'),
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords)
plot(g_mid_am, edge.width = E(g_rng_am)$weight,
     edge.color = ifelse(E(g_mid_am)$weight < 0.11, 'transparent',
                         rgb(0, 0, 0, 0.05)), 
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords, add = TRUE)

plot(g_mid_am, edge.width = 3*E(g_mid_am)$weight,
     edge.color = ifelse(E(g_mid_am)$weight < 0.4, 'transparent','black'),
     layout = coords, vertex.size = 3, vertex.label.color = 'transparent')
plot(g_mid_am, edge.width = 3*E(g_rng_am)$weight,
     edge.color = ifelse(E(g_mid_am)$weight < 0.4, 'transparent',
                         rgb(0, 0, 0, 0.25)), 
     vertex.label = ages$num,
     vertex.color = ifelse(ages$age_cat_id_1 == 5, 'purple',
                           ifelse(ages$age_cat_id_1 == 6, 'blue', 'seagreen')),
     vertex.label.dist = 0,
     vertex.label.color = 'white',
     vertex.size = 8,
     layout = coords, add = TRUE)

which(E(g_mid_am)$weight == 0.9713955) # 352

which(adj_mid_am == 0.9713955) # 8553
adj_mid_am[151] # 150 rows, 8553 = elephant 3 with elephant 58

sort(unique(counts_df_am$id_1))[c(3,58)]
counts_df[counts_df$dyad == 'M100_M2',] # ONLY SEEN TOGETHER ONCE


raw <- read_csv('data_processed/motnp_eles_long_22.01.06.csv')
length(which(raw$elephant == 'M2'))   # 17 sightings
length(which(raw$elephant == 'M100')) # 12 sightings
raw$encounter[which(raw$elephant == 'M2')] # 289 696 344 345 496  17 382 746 822  90 515 146 738 748 844 186 203
raw$encounter[which(raw$elephant == 'M100')] # 342 705 463 563 792 830 253 793 129 843 829 738
(17+12) - length(unique(c(raw$encounter[which(raw$elephant == 'M2')], raw$encounter[which(raw$elephant == 'M100')]))) # 1 == correct count in the input data for number of sightings together, but must be incorrect calculation for number of times seen apart


counts_ls$together[which(counts_df$dyad == 'M100_M2')]  # 1
counts_ls$apart[which(counts_df$dyad == 'M100_M2')]     # 27 -- should be 28, but that's another issue, if the right data are going in, why is it giving such a high reading here?

colnames(draws)
draws[,'M100_M2']
median(draws[,'M100_M2']) # 0.057 -- this is definitely more like correct, it's the plot which is the issue

means$dyad[which(means$median == 0.9713955)] # F1_M6 NOT M100_M2




