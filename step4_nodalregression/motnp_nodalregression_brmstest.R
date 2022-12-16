#### Set up ####
# load packages
library(tidyverse)
library(dplyr)
library(rstan)
library(rethinking)
library(igraph)
library(dagitty)
library(cmdstanr)
library(bisonR)

# information
sessionInfo()
R.Version()
rstan::stan_version()

# set stan path
#set_cmdstan_path('/Users/helen/.cmdstanr/cmdstan-2.28.2')

# set seed
set.seed(12345)

#### create data lists ####
### import data for aggregated model (binomial)
#counts_df <- read_csv('../../../../Google Drive/Shared drives/Helen PhD/chapter1_age/data_processed/motnp_bayesian_binomialpairwiseevents.csv')
counts_df <- read_csv('../data_processed/motnp_bayesian_binomialpairwiseevents.csv')


# correct sex_1, which has loaded in as a logical vector not a character/factor
unique(counts_df$sex_1) # FALSE or NA
sex_1 <- data.frame(sex_1 = counts_df$id_1)
sex_1 <- sex_1 %>% separate(sex_1, c("sex", "number"), sep = 1, remove = FALSE) ; unique(sex_1$sex) # F, M, U
counts_df$sex_1 <- as.character(sex_1$sex) ; rm(sex_1)
str(counts_df)  # sex_1 still comes up as logical, but now contains the right levels

# create variable for age difference
unique(counts_df$age_category_1) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
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
counts_df$age_cat_id_1[which(is.na(counts_df$age_cat_id_1))] <- 1 # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1

counts_df$age_class_1 <- ifelse(counts_df$age_cat_id_1 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_1 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_1 > 4, 'Adult','Pubescent')))

unique(counts_df$age_category_2) # "0-3","1-2","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-15","15-19","20-25","20-35","25-40","35-50","40+","50+",NA 
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
counts_df$age_cat_id_2[which(is.na(counts_df$age_cat_id_2))] <- 1   # U8 doesn't have an age but not a problem -- won't be part of the main analysis. Is a calf so age_cat_id = 1

counts_df$age_class_2 <- ifelse(counts_df$age_cat_id_2 == 1, 'Calf',
                                ifelse(counts_df$age_cat_id_2 == 2, 'Juvenile',
                                       ifelse(counts_df$age_cat_id_2 > 4, 'Adult','Pubescent')))

### correct dem_class with corrected age classes
counts_df$dem_class_1 <- ifelse(counts_df$age_class_1 == 'Adult', paste('A',counts_df$sex_1, sep = ''),
                                ifelse(counts_df$age_class_1 == 'Pubescent', paste('P',counts_df$sex_1, sep = ''),
                                       ifelse(counts_df$age_class_1 == 'Juvenile', paste('J',counts_df$sex_1, sep = ''),
                                              paste('C',counts_df$sex_1, sep = ''))))
counts_df$dem_class_2 <- ifelse(counts_df$age_class_2 == 'Adult', paste('A',counts_df$sex_2, sep = ''),
                                ifelse(counts_df$age_class_2 == 'Pubescent', paste('P',counts_df$sex_2, sep = ''),
                                       ifelse(counts_df$age_class_2 == 'Juvenile', paste('J',counts_df$sex_2, sep = ''),
                                              paste('C',counts_df$sex_2, sep = ''))))

### combined dem_class of dyad
counts_df$age_class_id_1 <- ifelse(counts_df$age_class_1 == 'Adult',4,
                                   ifelse(counts_df$age_class_1 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_1 == 'Juvenile',2,1)))
counts_df$age_class_id_2 <- ifelse(counts_df$age_class_2 == 'Adult',4,
                                   ifelse(counts_df$age_class_2 == 'Pubescent',3,
                                          ifelse(counts_df$age_class_2 == 'Juvenile',2,1)))
counts_df$dem_type <- ifelse(counts_df$age_class_id_1 > counts_df$age_class_id_2,
                             paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'), # when 1 is older: dc1_dc2
                             ifelse(counts_df$age_class_id_1 < counts_df$age_class_id_2,
                                    paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_'), # when 2 is older: dc2_dc1
                                    ifelse(counts_df$sex_1 == 'F',                                  # when age1 = age2...
                                           paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'), # ...when 1 is F: dc1_dc2
                                           ifelse(counts_df$sex_2 == 'F',
                                                  paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_'), # ...when 2 is F: dc2_dc1
                                                  ifelse(counts_df$sex_1 == 'M', # ...when neither are F...
                                                         paste(counts_df$dem_class_1, counts_df$dem_class_2, sep = '_'), # ...when 1 is M: dc1_dc2
                                                         paste(counts_df$dem_class_2, counts_df$dem_class_1, sep = '_')))))) # ...when 2 is M or both are U: dc2_dc1
sort(unique(counts_df$dem_type))

### add column for age difference between dyad
counts_df$age_diff <- abs(as.numeric(counts_df$age_cat_id_1) - as.numeric(counts_df$age_cat_id_2))

### add column for total number of sightings per pair
counts_df$count_dyad <- (counts_df$count_1 + counts_df$count_2) - counts_df$event_count  # maximum possible sightings of pairing = sum of times see node_1 and times see node_2, but this includes the times they were seen together twice, so then subtract once the count of paired sightings.

### add column for total number of sightings per pair where they were NOT together
counts_df$apart <- counts_df$count_dyad - counts_df$event_count

# reassign dyad numbers to remove gaps
counts_df$node_1_nogaps <- as.integer(as.factor(counts_df$node_1))
counts_df$node_2_nogaps <- as.integer(as.factor(counts_df$node_2))+1

############### EDGE WEIGHT ESTIMATION -- WORKING ######
#### test simulated data ####
sim_data <- simulate_bison_model("binary", aggregated = TRUE)
df <- sim_data$df_sim
priors <- get_default_priors("binary")
prior_check(priors, "binary")
fit_edge <- bison_model(
  (event | duration) ~ dyad(node_1_id, node_2_id), 
  data=df, 
  model_type="binary",
  priors=priors
)
plot_trace(fit_edge, par_ids=2)
plot_predictions(fit_edge, num_draws=20, type="density")
plot_predictions(fit_edge, num_draws=20, type="point")
summary(fit_edge)
plot_network(fit_edge, lwd=5)
fit_null <- bison_model(
  (event | duration) ~ 1, 
  data=df, 
  model_type="binary",
  priors=priors
)
model_comparison(list(non_random_model=fit_edge, random_model=fit_null))

cv_samples <- extract_metric(fit_edge, "global_cv")
head(cv_samples)
plot(density(cv_samples))

#### recreate data structure -- 10 elephants only ####
str(df)
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
#counts_df_males <- counts_df[counts_df$dem_type == 'AM_AM' | counts_df$dem_type == 'AM_PM' | counts_df$dem_type == 'PM_PM',]

counts_df <- counts_df %>% 
  filter(id_1 %in% unique(motnp_ages$id)) %>% 
  filter(id_2 %in% unique(motnp_ages$id))

random_eles <- sample(counts_df$id_1, 10, replace = F)
cdf_test <- counts_df[counts_df$id_1 %in% random_eles,]
cdf_test <- cdf_test[cdf_test$id_2 %in% random_eles,]

if(nrow(cdf_test) == 36) { print('redraw random eles') }

cdf_test$age_1 <- NA
cdf_test$age_2 <- NA
cdf_test$age_difference <- NA
#cdf_test$age_1 <- list(as.numeric(rep(NA, 8000)))
#cdf_test$age_2 <- list(as.numeric(rep(NA, 8000)))
#cdf_test$age_difference <- list(as.numeric(rep(NA, 8000)))
for(i in 1:nrow(cdf_test)){
  individual_age1 <- motnp_ages[cdf_test$id_1[i] == motnp_ages$id,]
  individual_age1 <- individual_age1[!is.na(individual_age1$id),]
  individual_age2 <- motnp_ages[cdf_test$id_2[i] == motnp_ages$id,]
  individual_age2 <- individual_age2[!is.na(individual_age2$id),]
  cdf_test$age_1[i] <- mean(individual_age1$age)
  cdf_test$age_2[i] <- mean(individual_age2$age)
  cdf_test$age_difference[i] <- cdf_test$age_1[i] - cdf_test$age_2[i]
  #cdf_test$age_1[i] <- list(individual_age1$age)
  #cdf_test$age_2[i] <- list(individual_age2$age)
}
#for(i in 1:nrow(cdf_test)){
#  cdf_test$age_difference[i] <- list(abs(unlist(cdf_test$age_1[i]) - unlist(cdf_test$age_2[i])))
#}
rm(individual_age1, individual_age2)

str(cdf_test)
str(df)

cdf_test$node_1_id <- factor(as.integer(as.factor(cdf_test$node_1)), levels = 1:10)
cdf_test$node_2_id <- factor(as.integer(as.factor(cdf_test$node_2))+1, levels = 1:10)
counts_df_model <- cdf_test[, c('node_1_id','node_2_id','event_count','count_dyad','age_difference','age_1','age_2')] %>% distinct() #,'location_id'
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration','age_diff','age_1','age_2') #,'location'

str(counts_df_model)
str(df)

#### run edge weight model ####
priors <- get_default_priors('binary')
#priors$edge <- 'beta(2,2)'
prior_check(priors, 'binary')

fit_edge_test <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)
plot_trace(fit_edge_test, par_ids=2)
plot_predictions(fit_edge_test, num_draws=20, type="density")
plot_predictions(fit_edge_test, num_draws=20, type="point")
summary(fit_edge_test)
plot_network(fit_edge_test, lwd=5)

rm(fit_edge, fit_null, sim_data)

#### run without unnecessary variables ####
counts_df_short <- counts_df_model[,1:4]
fit_edge_short <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_short, 
  model_type = "binary",
  priors = priors
)

############### NODAL REGRESSION -- 10 elephants only = WORKING ######
#### recreate data structure ####
df_nodal <- distinct(counts_df_model[,c('node_1_id','age_1')])
colnames(df_nodal) <- c('node_2_id', 'age_2')
df_nodal <- rbind(df_nodal, counts_df_model[nrow(counts_df_model),c('node_2_id','age_2')])
colnames(df_nodal) <- c('node', 'age')

#### nodal regression with single age value -- WORKING ####
fit_nodal <- bison_brm(
  bison(node_eigen(node)) ~ age,   # eigenvector centrality ~ age
  fit_edge_short,                   # model to obtain edge weights -- do I need to have used this model format to make this work?
  #chains = 4,
  df_nodal                         # data frame containing node information?
)
summary(fit_nodal)

#### redo but with ages as posterior distribution -- WORKING ####
df_nodal$id <- sort(random_eles)
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
motnp_ages <- motnp_ages[motnp_ages$id %in% random_eles,]
motnp_ages$node <- as.factor(as.integer(as.factor(motnp_ages$id)))

fit_nodal <- bison_brm(
  bison(node_eigen(node)) ~ age,   # eigenvector centrality ~ age
  fit_edge_short,                  # model to obtain edge weights
  #chains = 4,
  motnp_ages                       # data frame containing node information
)
summary(fit_nodal)

fit_null <- bison_brm(
  bison(node_eigen(node)) ~ 1,
  fit_edge_test,
  motnp_ages,
  chains = 4
)
model_comparison(list(non_random_model = fit_nodal, random_model = fit_null))

############### DYADIC REGRESSION -- 10 elephants only = WORKING ######
#### test simulated data ####
df_dyadic <- df %>% distinct(node_1_id, node_2_id, age_diff)
df_dyadic

fit_dyadic <- bison_brm (
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff,
  fit_edge_short,
  df_dyadic,
  num_draws=5, # Small sample size for demonstration purposes -- what does this mean?? chains only 5 draws long?
  refresh=0
)
summary(fit_dyadic)
rm(fit_dyadic)

#### recreate data structure ####
str(df_dyadic)
counts_df_dyadic <- counts_df_model[,c('node_1_id','node_2_id', 'age_diff')]
str(counts_df_dyadic)

#### dyadic regression -- WORKING ####
fit_dyadic <- bison_brm (
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff,
  fit_edge_short,
  counts_df_dyadic,
  num_draws=5, # Small sample size for demonstration purposes
  refresh=0
)
summary(fit_dyadic)

#### redo but with ages as posterior distribution -- WORKING ####
cdf_dyadic <- rbind(counts_df_short, counts_df_short, counts_df_short, counts_df_short, counts_df_short, # 10 draws per pair
                    counts_df_short, counts_df_short, counts_df_short, counts_df_short, counts_df_short)
cdf_dyadic <- rbind(cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic,                          # 100 draws
                    cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic)
cdf_dyadic <- rbind(cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic,                          # 1000 draws
                    cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic)
cdf_dyadic <- rbind(cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic,                                      # 8000 draws
                    cdf_dyadic, cdf_dyadic, cdf_dyadic, cdf_dyadic)

motnp_ages <- motnp_ages[motnp_ages$id %in% random_eles,]
motnp_ages$node <- as.factor(as.integer(as.factor(motnp_ages$id)))

cdf_dyadic$dyad <- paste0(cdf_dyadic$node_1_id, '_', cdf_dyadic$node_2_id)
cdf_dyadic$age_1 <- NA ; cdf_dyadic$age_2 <- NA
for(i in 1:length(unique(cdf_dyadic$dyad))){
  if(is.na(cdf_dyadic$age_1[i])){
    dyad <- cdf_dyadic[cdf_dyadic$dyad == cdf_dyadic$dyad[i],]
    cdf_dyadic <- anti_join(cdf_dyadic, dyad)
    id1 <- motnp_ages[motnp_ages$node == cdf_dyadic$node_1_id[i],]
    id2 <- motnp_ages[motnp_ages$node == cdf_dyadic$node_2_id[i],]
    dyad$age_1 <- id1$age
    dyad$age_2 <- id2$age
    cdf_dyadic <- rbind(cdf_dyadic, dyad)
  }
}
cdf_dyadic$age_diff <- cdf_dyadic$age_1 - cdf_dyadic$age_2

cdf_dyadic <- cdf_dyadic[,c('node_1_id','node_2_id','age_diff')]

str(counts_df_dyadic)
str(cdf_dyadic)

fit_dyadic <- bison_brm (
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff,
  fit_edge_short,
  counts_df_dyadic,
  num_draws=5, # Small sample size for demonstration purposes
  refresh=0
)
summary(fit_dyadic)

############### RERUN ALL BUT WITHOUT SUBSETTING DATA ######
rm(list = ls()[which(ls() != 'counts_df')])
#### edge weights ####
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

counts_df <- counts_df %>% 
  filter(id_1 %in% unique(motnp_ages$id)) %>% 
  filter(id_2 %in% unique(motnp_ages$id))

counts_df$node_1_id <- factor(counts_df$node_1_nogaps, levels = 1:(max(counts_df$node_1_nogaps)+1))
counts_df$node_2_id <- factor(counts_df$node_2_nogaps, levels = 1:(max(counts_df$node_1_nogaps)+1))
counts_df_model <- counts_df[, c('node_1_id','node_2_id','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

priors <- get_default_priors('binary')
prior_check(priors, 'binary')
fit_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)
plot_trace(fit_edge_weights, par_ids  =2)
plot_predictions(fit_edge_weights, num_draws = 20, type = "density")
plot_predictions(fit_edge_weights, num_draws = 20, type = "point")
summary(fit_edge_weights)
plot_network(fit_edge_weights, lwd = 5)

#### nodal regression ####
df_nodal <- distinct(counts_df[,c('node_1_id','id_1')])
colnames(df_nodal) <- c('node_2_id','id_2')
df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_id','id_2')])
colnames(df_nodal) <- c('node','id')

motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')

fit_nodal <- bison_brm(
  bison(node_eigen(node)) ~ age,
  fit_edge_weights,
  #chains = 4,
  motnp_ages
)
summary(fit_nodal)

#### dyadic regression ####
counts_df_dyadic <- counts_df[,c('dyad_id','node_1_id','node_2_id')]

counts_df_dyadic <- rbind(counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, # 10 draws each
                          counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic)
counts_df_dyadic <- rbind(counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, # 100 draws
                          counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic)
counts_df_dyadic <- rbind(counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, # 1000 draws
                          counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic)
counts_df_dyadic <- rbind(counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic,                   # 8000 draws
                          counts_df_dyadic, counts_df_dyadic, counts_df_dyadic, counts_df_dyadic)

counts_df_dyadic$age_1 <- NA ; counts_df_dyadic$age_2 <- NA
for(i in 1:length(unique(counts_df_dyadic$dyad_id))){
  if(is.na(counts_df_dyadic$age_1)){
    dyad <- counts_df_dyadic[counts_df_dyadic$dyad_id == counts_df_dyadic$dyad_id[i],]
    counts_df_dyadic <- anti_join(counts_df_dyadic, dyad)
    id1 <- motnp_ages[motnp_ages$node == counts_df_dyadic$node_1_id[i],]
    id2 <- motnp_ages[motnp_ages$node == counts_df_dyadic$node_2_id[i],]
    dyad$age_1 <- id1$age
    dyad$age_2 <- id2$age
    counts_df_dyadic <- rbind(counts_df_dyadic, dyad)
  }
}
counts_df_dyadic$age_diff <- counts_df_dyadic$age_1 - counts_df_dyadic$age_2

counts_df_dyadic <- counts_df_dyadic[,c('node_1_id','node_2_id','age_diff')]

fit_dyadic <- bison_brm (
  bison(edge_weight(node_1_id, node_2_id)) ~ age_diff,
  fit_edge_short,
  counts_df_dyadic,
  num_draws=5, # Small sample size for demonstration purposes
  refresh=0
)
summary(fit_dyadic)

