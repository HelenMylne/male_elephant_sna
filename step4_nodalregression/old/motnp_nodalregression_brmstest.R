#### Set up ####
# load packages
library(tidyverse)
library(dplyr)
library(rstan)
library(igraph)
library(dagitty)
library(cmdstanr)
library(bisonR)
library(brms)
library(mclust)

#library(tidyverse, lib.loc = '../packages/')
#library(dplyr, lib.loc = '../packages/')
#library(rstan, lib.loc = '../packages/')
#library(rethinking, lib.loc = '../packages/')
#library(igraph, lib.loc = '../packages/')
#library(dagitty, lib.loc = '../packages/')
#library(cmdstanr, lib.loc = '../packages/')
#library(bisonR, lib.loc = '../packages/')


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
priors$edge <- 'normal(-1,2)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
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

############### NETWORK COMPLEXITY ###############
# Edge mixture models can be used to detect different types of social connections by looking for clustering in edge weights. There are two key outputs from a bisonR edge mixture model: a posterior probability distribution over the number of components (clusters, interpreted associal connection types, from 1 to a maximum 𝐾, over the entire network), the probability that each dyad belongs to each possible component.

# simulate data from a binary edge model with 2 underlying components (dyads belong to either cluster with equal probability):
sim_data <- simulate_bison_model_mixture("binary", num_components = 2, component_weights = c(0.5, 0.5))
df <- sim_data$df_sim
head(df)
#> # A tibble: 6 × 5
#> # Groups:   node_1_id [1]
#>   node_1_id node_2_id event duration age_diff
#>   <fct>     <fct>     <dbl>    <dbl>    <dbl>
#> 1 1         2             1       14        1
#> 2 1         3             0       15        1
#> 3 1         4             8       10        1
#> 4 1         5            14       14        1
#> 5 1         6             0        1        1
#> 6 1         7            16       16        1
Since the data are already aggregated (e.g. a single row for each dyad), we can use the binary conjugate model for speed:
  
  fit_edge <- bison_model(
    (event | duration) ~ dyad(node_1_id, node_2_id), 
    data=df, 
    model_type="binary_conjugate"
  )
#> No priors set by user, using default priors instead. We recommend setting and checking priors explicitly for reliable inference.
We won’t conduct the usual checks for the sake of this tutorial, but at this point we’d recommend conducting diagnostic checks to ensure model fit. Once the model has fit well, we can use the bison_mixture() function to fit the edge mixture model:
  
  # verbose=FALSE for tutorial purposes
  fit_mixture <- bison_mixture(fit_edge, num_components=5, verbose=FALSE)
summary(fit_mixture)
#> === Fitted dyadic mixture model ===
#> Maximum number of components: 
#> Best model: 1 components
#> Probability of best model: 94.1%
#> === Component probabilities ===
#>            1     2     3 4 5
#> P(K=k) 0.941 0.052 0.007 0 0
#> === Component means for best model (K = 1) ===
#>              1
#> mean 0.5252286
#> === Edge component probabilities for best model (K = 1) ===
#>         1
#> 1 <-> 2 1
#> 1 <-> 3 1
#> 1 <-> 4 1
#> 1 <-> 5 1
#> 1 <-> 6 1
#> 1 <-> 7 1
#> ...
There’s quite a bit going on in the summary, so let’s break it down. The top section shows the model that fits best and the probability assigned to that model. The second section shows the probabilities of the models corresponding to each number of components (𝐾=1,2,3,...). You’ll notice the model with the highest probability is the same as the best fitting model from the section above. The final section shows the edge-level probabilities of component membership for the corresponding dyad. This can be interpreted as the probability a dyad belongs to a given component/cluster/social connection type.

The information shown in the summary above can also be accessed using the functions get_network_component_probabilities() and get_edge_component_probabilities(). These two functions make it possible to access the mixture model output programmatically for downstream analysis or plotting.

To get network-level probabilities of the number of components, we can use get_network_component_probabilities() function:
  
  get_network_component_probabilities(fit_mixture)
#>   num_components  probability
#> 1              1 9.411516e-01
#> 2              2 5.187969e-02
#> 3              3 6.532973e-03
#> 4              4 4.318312e-04
#> 5              5 3.869330e-06
To use probabilities over component membership for a given number of components, we can use the get_edge_component_probabilities() function, where 3 is the number of components we assume to exist in the network:
  
  get_edge_component_probabilities(fit_mixture, 3)
#>     node_1 node_2 P(K = 1) P(K = 2) P(K = 3)
#> 1        1      2    0.840    0.160    0.000
#> 2        1      3    0.875    0.115    0.010
#> 3        1      4    0.055    0.455    0.490
#> 4        1      5    0.010    0.245    0.745
#> 5        1      6    0.560    0.320    0.120
#> 6        1      7    0.005    0.225    0.770
#> 7        1      8    0.750    0.240    0.010
#> 8        1      9    0.905    0.090    0.005
#> 9        1     10    0.820    0.180    0.000
#> 10       1     11    0.780    0.205    0.015
#> 11       1     12    0.015    0.340    0.645
#> 12       1     13    0.015    0.205    0.780
#> 13       1     14    0.865    0.135    0.000
#> 14       1     15    0.010    0.425    0.565
#> 15       2      3    0.850    0.140    0.010
#> 16       2      4    0.740    0.245    0.015
#> 17       2      5    0.005    0.225    0.770
#> 18       2      6    0.000    0.330    0.670
#> 19       2      7    0.000    0.290    0.710
#> 20       2      8    0.840    0.155    0.005
#> 21       2      9    0.815    0.180    0.005
#> 22       2     10    0.105    0.445    0.450
#> 23       2     11    0.810    0.180    0.010
#> 24       2     12    0.005    0.255    0.740
#> 25       2     13    0.815    0.175    0.010
#> 26       2     14    0.760    0.230    0.010
#> 27       2     15    0.005    0.240    0.755
#> 28       3      4    0.805    0.185    0.010
#> 29       3      5    0.800    0.200    0.000
#> 30       3      6    0.025    0.395    0.580
#> 31       3      7    0.690    0.295    0.015
#> 32       3      8    0.120    0.470    0.410
#> 33       3      9    0.015    0.365    0.620
#> 34       3     10    0.010    0.260    0.730
#> 35       3     11    0.010    0.185    0.805
#> 36       3     12    0.005    0.195    0.800
#> 37       3     13    0.020    0.410    0.570
#> 38       3     14    0.845    0.150    0.005
#> 39       3     15    0.890    0.100    0.010
#> 40       4      5    0.020    0.265    0.715
#> 41       4      6    0.000    0.285    0.715
#> 42       4      7    0.015    0.220    0.765
#> 43       4      8    0.005    0.415    0.580
#> 44       4      9    0.755    0.235    0.010
#> 45       4     10    0.000    0.330    0.670
#> 46       4     11    0.000    0.395    0.605
#> 47       4     12    0.000    0.295    0.705
#> 48       4     13    0.140    0.410    0.450
#> 49       4     14    0.845    0.155    0.000
#> 50       4     15    0.015    0.345    0.640
#> 51       5      6    0.015    0.335    0.650
#> 52       5      7    0.885    0.115    0.000
#> 53       5      8    0.565    0.385    0.050
#> 54       5      9    0.060    0.450    0.490
#> 55       5     10    0.850    0.150    0.000
#> 56       5     11    0.010    0.220    0.770
#> 57       5     12    0.130    0.455    0.415
#> 58       5     13    0.370    0.475    0.155
#> 59       5     14    0.165    0.425    0.410
#> 60       5     15    0.010    0.185    0.805
#> 61       6      7    0.780    0.220    0.000
#> 62       6      8    0.810    0.185    0.005
#> 63       6      9    0.830    0.165    0.005
#> 64       6     10    0.000    0.330    0.670
#> 65       6     11    0.005    0.220    0.775
#> 66       6     12    0.850    0.150    0.000
#> 67       6     13    0.790    0.205    0.005
#> 68       6     14    0.010    0.275    0.715
#> 69       6     15    0.010    0.355    0.635
#> 70       7      8    0.775    0.215    0.010
#> 71       7      9    0.245    0.525    0.230
#> 72       7     10    0.890    0.105    0.005
#> 73       7     11    0.005    0.390    0.605
#> 74       7     12    0.690    0.280    0.030
#> 75       7     13    0.665    0.255    0.080
#> 76       7     14    0.565    0.325    0.110
#> 77       7     15    0.875    0.120    0.005
#> 78       8      9    0.005    0.340    0.655
#> 79       8     10    0.000    0.385    0.615
#> 80       8     11    0.840    0.150    0.010
#> 81       8     12    0.675    0.300    0.025
#> 82       8     13    0.830    0.165    0.005
#> 83       8     14    0.020    0.450    0.530
#> 84       8     15    0.005    0.295    0.700
#> 85       9     10    0.185    0.425    0.390
#> 86       9     11    0.875    0.120    0.005
#> 87       9     12    0.740    0.250    0.010
#> 88       9     13    0.840    0.140    0.020
#> 89       9     14    0.005    0.250    0.745
#> 90       9     15    0.875    0.120    0.005
#> 91      10     11    0.035    0.395    0.570
#> 92      10     12    0.220    0.510    0.270
#> 93      10     13    0.010    0.230    0.760
#> 94      10     14    0.145    0.405    0.450
#> 95      10     15    0.120    0.550    0.330
#> 96      11     12    0.015    0.335    0.650
#> 97      11     13    0.120    0.455    0.425
#> 98      11     14    0.615    0.330    0.055
#> 99      11     15    0.830    0.170    0.000
#> 100     12     13    0.000    0.330    0.670
#> 101     12     14    0.005    0.360    0.635
#> 102     12     15    0.835    0.160    0.005
#> 103     13     14    0.715    0.275    0.010
#> 104     13     15    0.865    0.135    0.000
#> 105     14     15    0.565    0.365    0.070
It’s often useful to know where the mixture components are location, so we can get their posterior means using the get_component_means() function for a given number of components. Here we’ll calculate the posterior means of the components for the mixture model with 3 components:
  
  get_component_means(fit_mixture, 3)
#>             50%         5%       95%
#> K = 1 0.1045326 0.01339178 0.1966635
#> K = 2 0.7640626 0.13868448 0.8916400
#> K = 3 0.9176529 0.84222298 0.9981651
Using components in downstream analysis
Now say that we’re interested in the number of strong vs weak partners. Specifically we want to know if the strength of connections only to strong partners predicts a node-level trait. This is a non-standard analysis so will require some manual data wrangling to fit the right kind of model.

For demonstration purposes, we’ll create a dataframe that might be used for this kind of analysis:
  
  df_nodal <- data.frame(node_id=as.factor(levels(df$node_1_id)), trait=rnorm(15))
df_nodal
#>    node_id       trait
#> 1        1 -1.35866457
#> 2        2  0.18107587
#> 3        3 -0.71568046
#> 4        4 -1.07174483
#> 5        5  0.88888720
#> 6        6 -0.57613074
#> 7        7 -1.62117078
#> 8        8  1.94552654
#> 9        9 -0.73863301
#> 10      10 -1.14172964
#> 11      11  1.33588630
#> 12      12  0.81998647
#> 13      13 -0.57319856
#> 14      14  0.01948893
#> 15      15  1.31491191
This dataframe will provide the information about traits that we need for the analysis. Next, we need to calculate strength using only individuals categorised as strong according to the mixture model. Since both component membership (strong or weak) and edge weights are probabilistic, we have to sample from both posteriors as we build a new posterior of mixture-based strength:
  
  # Set number of nodes and number of posterior samples
  num_nodes <- fit_edge$num_nodes
num_draws <- 5 # Keep short for demo purposes.

# Create a list of igraph networks from edgemodel to represent network posterior
nets <- bison_to_igraph(fit_edge, num_draws)

# Create an empty matrix to hold strengths of top mixture component
mix_strengths <- matrix(0, num_draws, num_nodes)

# Create an empty list for imputed versions of the dataframe
imputed_data <- list()

# Loop through each posterior sample and calculate strength of top mixture.
for (i in 1:num_draws) {
  # Calculate edge components (1 if strong, 0 if weak)
  edge_components <- 1 * (fit_mixture$edge_component_samples[i, 2, ] == 2)
  
  # If the edge is strong, don't change edge weight, but if it's weak then set to zero.
  E(nets[[i]])$weight <- E(nets[[i]])$weight * edge_components
  
  # Change the metric values of the imputed data to be mixture-based strength
  imputed_data[[i]] <- df_nodal
  imputed_data[[i]]$mix_strength <- strength(nets[[i]])
}
imputed_data[[2]]
#>    node_id       trait mix_strength
#> 1        1 -1.35866457     5.781460
#> 2        2  0.18107587     5.560490
#> 3        3 -0.71568046     6.031965
#> 4        4 -1.07174483     7.842593
#> 5        5  0.88888720     8.071741
#> 6        6 -0.57613074     7.619625
#> 7        7 -1.62117078     4.545718
#> 8        8  1.94552654     5.506187
#> 9        9 -0.73863301     3.732620
#> 10      10 -1.14172964     8.032430
#> 11      11  1.33588630     6.719255
#> 12      12  0.81998647     7.499526
#> 13      13 -0.57319856     4.383482
#> 14      14  0.01948893     6.334001
#> 15      15  1.31491191     6.358431
Now we have a set of imputed dataframes containing the posteriors of mixture-based node strength, we can fit the model using brm_multiple.

fit_brm <- brm_multiple(trait ~ mix_strength, imputed_data, silent=2, refresh=0)
summary(fit_brm)


############### RERUN ALL BUT WITHOUT SUBSETTING DATA ######
rm(list = ls()[which(ls() != 'counts_df')])
#### edge weights ####
motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')

counts_df <- counts_df %>% 
  filter(id_1 %in% unique(motnp_ages$id)) %>% 
  filter(id_2 %in% unique(motnp_ages$id))

#counts_df$node_1_id <- factor(counts_df$node_1_males, levels = 1:(max(counts_df$node_1_nogaps)+1))
#counts_df$node_2_id <- factor(counts_df$node_2_males, levels = 1:(max(counts_df$node_1_nogaps)+1))
#counts_df_model <- counts_df[, c('node_1_id','node_2_id','event_count','count_dyad')] %>% distinct()
#colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

counts_df$node_1_males <- as.integer(as.factor(counts_df$node_1_nogaps))
counts_df$node_2_males <- as.integer(as.factor(counts_df$node_2_nogaps))+1
counts_df_model <- counts_df[, c('node_1_males','node_2_males','event_count','count_dyad')] %>% distinct()
colnames(counts_df_model) <- c('node_1_id','node_2_id','event','duration')

priors <- get_default_priors('binary')
priors$edge <- 'normal(-1, 2)'    # slightly right skewed so that more density is towards 0 but still very weak and all values are possible
priors$fixed <- 'normal(0, 1.5)'  # slightly less wide that before, but zero centred so allows age to have a positive or negative impact on centrality and edge weight
prior_check(priors, 'binary')

fit_edge_weights <- bison_model(
  ( event | duration ) ~ dyad(node_1_id, node_2_id), 
  data = counts_df_model, 
  model_type = "binary",
  #mc_cores = 4,
  priors = priors
)
plot_trace(fit_edge_weights, par_ids = 2)
plot_predictions(fit_edge_weights, num_draws = 20, type = "density")
plot_predictions(fit_edge_weights, num_draws = 20, type = "point")

fit_edge_weights$chain         # columns = dyads, rows = chain position --> all draws of edge weight (untransformed)
fit_edge_weights$edge_samples  # columns = dyads, rows = chain position --> all draws of edge weight (appears identical to chain)
which(fit_edge_weights$chain != fit_edge_weights$edge_samples) # they're the same
fit_edge_weights$event_preds   # columns = dyads, rows = chain position --> all draws of number of times seen together (out of..?)

edgelist <- get_edgelist(fit_edge_weights, ci = 0.9, transform = TRUE)
summary(fit_edge_weights)  # THIS APPEARS TO BE CONFUSED ABOUT NUMBER OF NODES BUT IT'S NOT -- node_id values run 152-396 (151 females before males, 60ish unknowns after males)

# save workspace image for reloading at a later date that doesn't require running model again
save.image('motnp_bisonr_edgescalculated.RData')

#### plot network ####
load('motnp_bisonr_edgescalculated.RData')

plot_network(fit_edge_weights, lwd = 2, ci = 0.9) # this is all of the options for arguments in plot_network() -- currently illegible as nodes are too large. can this be saved as an igraph object?

plot_network
#function (obj, ci = 0.9, lwd = 2) 
#{
#  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
#  net <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), 
#                                     directed = obj$directed)
#  lb <- edgelist[, 3]
#  ub <- edgelist[, 5]
#  coords <- igraph::layout_nicely(net)
#  igraph::plot.igraph(net, edge.width = lb * lwd, layout = coords, 
#                      vertex.label.color = "white", vertex.color = bison_colors[1], 
#                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0)
#  igraph::plot.igraph(net, edge.width = ub * lwd, layout = coords, 
#                      vertex.label.color = "white", vertex.color = bison_colors[1], 
#                      edge.color = rgb(0, 0, 0, 0.3), edge.arrow.size = 0, 
#                      add = TRUE)
#  if (obj$directed) {
#    igraph::plot.igraph(net, edge.width = lb * 0, layout = coords, 
#                        vertex.label.color = "white", vertex.color = bison_colors[1], 
#                        edge.color = rgb(0.5, 0.5, 0.5), add = TRUE)
#  }
#}

plot_network_threshold <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                    vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                    vertex.color1 = 'transparent',
                                    vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                    vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                    vertex.color2 = 'seagreen1',
                                    vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  net <- igraph::graph_from_edgelist(as.matrix(edgelist[, 1:2]), 
                                     directed = obj$directed)
  md <- edgelist[, 3]
  ub <- edgelist[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      edge.width = ifelse(md < threshold, 0, md * lwd),
                      vertex.label.color = vertex.label.color1, vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0)
  igraph::plot.igraph(net, layout = coords,
                      edge.width = ifelse(md < threshold, 0, ub * lwd),
                      vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, 
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

plot_network_threshold(fit_edge_weights, ci = 0.9, threshold = 0.3,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black',
                       vertex.size2 = 7, edge.color2 = 'black')

nodes <- data.frame(bull = sort(unique(motnp_ages$id)),
                    age = NA,
                    sightings = NA)
for(i in 1:nrow(nodes)){
  x <- motnp_ages[motnp_ages$id == nodes$bull[i],]
  nodes$age[i] <- mean(x$age)
  if(nodes$bull[i] != 'M99'){
    y <- counts_df[counts_df$id_1 == nodes$bull[i], c('id_1','count_1')]
    nodes$sightings[i] <- y[1,2]
  }
}
nodes$sightings <- as.numeric(nodes$sightings)
str(nodes)

plot_network_threshold(fit_edge_weights)

plot_network_threshold(fit_edge_weights, lwd = 2, ci = 0.9, threshold = 0.1,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = 7, edge.color2 = 'black')

plot_network_threshold(fit_edge_weights, lwd = 2, ci = 0.9, threshold = 0.1,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

plot_network_threshold(fit_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')

plot_network_threshold(fit_edge_weights, lwd = 2, ci = 0.9, threshold = 0.2,
                       vertex.label.color1 = NA, edge.color1 = rgb(0,0,0,0.25),
                       vertex.label.color2 = 'black', vertex.color2 = nodes$age,
                       vertex.size2 = nodes$sightings, edge.color2 = 'black')


plot_network_threshold2 <- function (obj, ci = 0.9, lwd = 2, threshold = 0.3,
                                     vertex.label.color1 = 'transparent', vertex.label.font1 = NA, 
                                     vertex.color1 = 'transparent',
                                     vertex.size1 = 1, edge.color1 = rgb(0, 0, 0, 1),
                                     vertex.label.color2 = 'black', vertex.label.font2 = 'Helvetica', 
                                     vertex.color2 = 'seagreen1',
                                     vertex.size2 = 8, edge.color2 = rgb(0, 0, 0, 0.3))
{
  edgelist <- get_edgelist(obj, ci = ci, transform = TRUE)
  threshold_edges <- edgelist[edgelist$median >= threshold,]
  net <- igraph::graph_from_edgelist(as.matrix(threshold_edges[, 1:2]))
  md <- threshold_edges[, 3]
  ub <- threshold_edges[, 5]
  coords <- igraph::layout_nicely(net)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color1, label.family = vertex.label.font1, 
                      vertex.color = vertex.color1, vertex.size = vertex.size1,
                      edge.color = rgb(0, 0, 0, 1), edge.arrow.size = 0, edge.width = md * lwd)
  igraph::plot.igraph(net, layout = coords,
                      vertex.label.color = vertex.label.color2, label.family = vertex.label.font2,
                      vertex.color = vertex.color2, vertex.size = vertex.size2,
                      edge.color = edge.color1, edge.arrow.size = 0, edge.width = ub * lwd,
                      add = TRUE)
  if (obj$directed) {
    igraph::plot.igraph(net, edge.width = md * 0, layout = coords, 
                        vertex.label.color = vertex.label.color2, vertex.color = vertex.color2, 
                        edge.color = edge.color2, add = TRUE)
  }
}

plot_network_threshold2(obj = fit_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = 8)

plot_network_threshold2(obj = fit_edge_weights, threshold = 0.2,
                        vertex.color2 = nodes$age, vertex.size2 = nodes$sightings)

#### non-random edge weights ####
fit_edges_null <- bison_model(
  (event | duration) ~ 1, 
  data = counts_df_model, 
  model_type = "binary",
  priors = priors
)

model_comparison(list(non_random_model = fit_edge_weights, random_model = fit_edges_null)) # NETWORK IS RANDOM
# Method: stacking
#                  weight
# non_random_model 0.020
# random_model     0.980
# Warning message:  Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

# save workspace image
save.image('motnp_bisonr_randomnetwork.RData')

#### nodal regression -- mean age ####
df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
colnames(df_nodal) <- c('node_2_males','id_2')
df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
colnames(df_nodal) <- c('node','id')

motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')

motnp_ages_mean <- df_nodal
motnp_ages_mean$age <- NA
for(i in 1:nrow(motnp_ages_mean)){
  x <- motnp_ages[motnp_ages$id == motnp_ages_mean$id[i],]
  motnp_ages_mean$age[i] <- mean(x$age)
  rm(x)
}

fit_nodal_mean <- bison_brm(
  bison(node_eigen(node)) ~ age,
  fit_edge_weights,
  motnp_ages_mean,
  chains = 4
)
summary(fit_nodal_mean)

# save workspace image for reloading at a later date that doesn't require running model again
save.image('motnp_bisonr_nodalregression_meanage.RData')

# remove nodal mean age
rm(fit_nodal_mean)

#### nodal regression -- age distribution ####
fit_nodal <- bison_brm(
  bison(node_eigen(node)) ~ age,
  fit_edge_weights,
  motnp_ages,
  chains = 4
)
summary(fit_nodal)

save.image('motnp_bisonr_nodalregressionrun.RData')

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


