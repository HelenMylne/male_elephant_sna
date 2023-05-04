#### Information ####
# script takes data input from edge weight estimation for MOTNP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### Set up ####
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(car, lib.loc = '../packages/')       # library(car)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(bisonR, lib.loc = '../packages/')    # library(bisonR)
library(brms, lib.loc = '../packages/')      # library(brms)
#library(rstan, lib.loc = '../packages/')     # library(rstan)
#library(Rcpp, lib.loc = '../packages/')      # library(Rcpp)

#set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan/cmdstan-2.31.0')

load('motnp_bisonr_edgescalculated_strongprior.RData')
rm(counts_df_model, motnp_edges_null_strongpriors) ; gc()

#pdf('../outputs/motnp_dyadicregression_plots.pdf')

#### create data frame of male ages ####
## unique ids
#ids <- counts_df[,c('id_1','node_1_males')] %>% distinct()
#colnames(ids) <- c('id_2','node_2_males')
#ids <- rbind(ids, counts_df[nrow(counts_df),c('id_2','node_2_males')])
#colnames(ids) <- c('id','node_males')

## male ages
#motnp_ages <- as.data.frame(readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds'))
#motnp_ages <- motnp_ages[,colnames(motnp_ages) %in% ids$id]
#ids$mean_age <- NA
#for(i in 1:nrow(ids)){
#  ids$mean_age[i] <- mean(motnp_ages[,colnames(motnp_ages) == ids$id[i]])
#}
#motnp_ages <- motnp_ages %>%
#  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age') %>% 
#  filter(id %in% ids$id) %>% 
#  left_join(ids, by = 'id')
#motnp_ages$draw <- rep(1:8000, each = length(ids$id))

## create data frame of ages
#cdf_dyadic <- counts_df[,c('dyad_males','node_1_males','node_2_males')] # counts_df = data frame containing all information about dyad (inc. counts of sightings together/apart)
#colnames(motnp_ages)[3] <- 'node_1_males'
#cdf_dyadic <- left_join(cdf_dyadic,
#                              distinct(motnp_ages[,c('id','node_1_males','mean_age')]),
#                              by = 'node_1_males')
#colnames(cdf_dyadic)[4:5] <- c('id_1','age_1')
#colnames(motnp_ages)[3] <- 'node_2_males'
#rm(counts_df) ; gc()
#cdf_dyadic <- left_join(cdf_dyadic,
#                              distinct(motnp_ages[,c('id','node_2_males','mean_age')]),
#                              by = 'node_2_males')
#colnames(cdf_dyadic)[6:7] <- c('id_2','age_2')
#cdf_dyadic <- cdf_dyadic[,c(1:4,6,5,7)]
#colnames(cdf_dyadic)[2:3] <- c('node_1_id','node_2_id')

#cdf_dyadic$min_age <- NA
#cdf_dyadic$max_age <- NA
#for(i in 1:nrow(cdf_dyadic)){
#  cdf_dyadic$min_age[i] <- min(c(cdf_dyadic$age_1[i], cdf_dyadic$age_2[i]))
#  cdf_dyadic$max_age[i] <- max(c(cdf_dyadic$age_1[i], cdf_dyadic$age_2[i]))
#}
#length(which(cdf_dyadic$min_age >= cdf_dyadic$max_age))

#cdf_dyadic$age_difference <- cdf_dyadic$max_age - cdf_dyadic$min_age

#head(cdf_dyadic)

#write_csv(cdf_dyadic, '../data_processed/motnp_dyadicregression_modeldata.csv')
cdf_dyadic <- read_csv('../data_processed/motnp_dyadicregression_modeldata.csv')

#### prior predictive check ####
# prior predictive check
#priors <- get_default_priors('binary')
#priors$fixed
#prior_check(priors, 'binary')
#
#age_1 <- 10:60
#age_2 <- 10:60
#plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
#     xlab = 'age', ylab = 'eigenvector centrality')
#for(i in 1:100){
#  intercept <- rbeta(1,1,1)
#  beta <- rnorm(1, 0, 0.005)
#  lines(x = age_1, y = intercept + age_1*beta + age_2*beta, col = 'blue')
#}
#
#age_1 <- 10:60
#plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
#     xlab = 'age', ylab = 'eigenvector centrality')
#for(j in 1:50) {
#intercept <- rbeta(1,1,1)
#beta <- rnorm(1, 0, 0.005)
#  for(i in seq(10,60,by=10)){
#    age_2 <- i
#    lines(x = age_1, y = intercept + age_1*beta + age_2*beta, col = rgb(0,0,1,0.5))
#  }
#}
#
#prior <- bison_brm_get_prior(
#  bison(edge_weight(node_1_males, node_2_males)) ~ age_1 + age_2 + (1 | mm(node_1_males, node_2_males)),
#  list(motnp_edge_weights_strongpriors),
#  list(cdf_dyadic)
#)
#brms::prior_draws(prior)

#### fit model to mean age ####
mean_age_dyadic <- bison_brm (
  #bison(edge_weight(node_1_id, node_2_id)) ~ min_age + max_age + min_age:max_age + (1 | mm(node_1_id, node_2_id)),
  bison(edge_weight(node_1_id, node_2_id)) ~ age_difference + (1 | mm(node_1_id, node_2_id)),
  motnp_edge_weights_strongpriors,
  cdf_dyadic,
  #num_draws = 5, # Small sample size for demonstration purposes
  #refresh = 0,
  cores = 4, 
  chains = 4,
  #control = list(max_treedepth = 20),
  iter = 1000
)
summary(mean_age_dyadic)

hist(mean_age_dyadic$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

save.image('motnp_dyadicregression_meanage.RData')
dev.off()
