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
rm(counts_df_model, motnp_edges_null_strongpriors, edgelist, females_df, priors, motnp_ages, model_averaging) ; gc()
#rm(random_networks, gbi_males, m_mat) ; gc()

#pdf('../outputs/motnp_dyadicregression_plots.pdf')

#### create data frame of male ages ####
## unique ids
ids <- counts_df[,c('id_1','node_1_males')] %>% distinct()
colnames(ids) <- c('id_2','node_2_males')
ids <- rbind(ids, counts_df[nrow(counts_df),c('id_2','node_2_males')])
colnames(ids) <- c('id','node_males')

## male ages
#motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>%
#  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age') %>% 
#  filter(id %in% ids$id) %>% 
#  left_join(ids, by = 'id')
#motnp_ages$draw <- rep(1:8000, each = length(ids$id))
motnp_ages <- as.data.frame(readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds'))
motnp_ages <- motnp_ages[,colnames(motnp_ages) %in% ids$id]
ids$mean_age <- NA
for(i in 1:nrow(ids)){
  ids$mean_age[i] <- mean(motnp_ages[,colnames(motnp_ages) == ids$id[i]])
}

## create data frame of ages
#counts_df_dyadic <- counts_df[,c('dyad_males','node_1_males','node_2_males')]
#colnames(motnp_ages)[3] <- 'node_1_males'
#counts_df_dyadic <- left_join(counts_df_dyadic, motnp_ages, by = 'node_1_males', multiple = 'all')
#colnames(counts_df_dyadic)[4:5] <- c('id_1','age_1')
#colnames(motnp_ages)[3] <- 'node_2_males'
#rm(counts_df) ; gc()
#counts_df_dyadic <- left_join(counts_df_dyadic, motnp_ages, by = c('node_2_males','draw'))
#colnames(counts_df_dyadic)[7:8] <- c('id_2','age_2')
#counts_df_dyadic <- counts_df_dyadic[,c(1:3,6,4:5,7:8)]
#colnames(counts_df_dyadic)[2:3] <- c('node_1_id','node_2_id')
#head(counts_df_dyadic)
cdf_dyadic <- counts_df[,c('dyad_males','node_1_males','node_2_males')] # counts_df = data frame of dyads and their counts of sightings together and apart (data that went into edgeweight model)
colnames(ids) <- c('id_1','node_1_males','age_1')
cdf_dyadic <- left_join(cdf_dyadic, ids, by = 'node_1_males')
colnames(ids) <- c('id_2','node_2_males','age_2')
cdf_dyadic <- left_join(cdf_dyadic, ids, by = 'node_2_males')
length(which(is.na(cdf_dyadic$age_1) == TRUE))
length(which(is.na(cdf_dyadic$age_2) == TRUE))

write_csv(cdf_dyadic, '../data_processed/motnp_dyadicregression_modeldata.csv')
#cdf_dyadic <- read_csv('../data_processed/motnp_dyadicregression_modeldata.csv')

#### prior predictive check ####
# prior predictive check
priors <- get_default_priors('binary')
priors$fixed
prior_check(priors, 'binary')

age_1 <- 10:60
age_2 <- 10:60
plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  intercept <- rbeta(1,1,1)
  beta <- rnorm(1, 0, 0.005)
  lines(x = age_1, y = intercept + age_1*beta + age_2*beta, col = 'blue')
}

age_1 <- 10:60
plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(j in 1:50) {
intercept <- rbeta(1,1,1)
beta <- rnorm(1, 0, 0.005)
  for(i in seq(10,60,by=10)){
    age_2 <- i
    lines(x = age_1, y = intercept + age_1*beta + age_2*beta, col = rgb(0,0,1,0.5))
  }
}

prior <- bison_brm_get_prior(
  bison(edge_weight(node_1_males, node_2_males)) ~ age_1 + age_2 + (1 | mm(node_1_males, node_2_males)),
  list(motnp_edge_weights_strongpriors),
  list(cdf_dyadic)
)
brms::prior_draws(prior)


#### fit model to mean age ####
mean_age_dyadic <- bison_brm (
  bison(edge_weight(node_1_males, node_2_males)) ~ age_1 + age_2 + (1 | mm(node_1_males, node_2_males)),
  motnp_edge_weights_strongpriors,
  cdf_dyadic,
  #num_draws = 5, # Small sample size for demonstration purposes
  #refresh = 0,
  cores = 4, 
  chains = 4,
  iter = 1000,
  control = list(max_treedepth = 20)
)
summary(mean_age_dyadic)

dyad_data <- mean_age_dyadic$data
plot(dyad_data$ ~ dyad_data$, las = 1, pch = 19, col = rgb(0,0,1,0.2),
     xlab = '', ylab = '',
     main = '')

#mod_summary <- mean_age_dyadic$fit

hist(mean_age_dyadic$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

save.image('motnp_dyadicregression_meanage.RData')
dev.off()
