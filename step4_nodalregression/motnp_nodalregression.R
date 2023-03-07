#### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(brms, lib.loc = '../packages/')            # library(brms)
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
library(posterior, lib.loc = '../packages/')       # library(posterior)
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
#library(rstan, lib.loc = '../packages/')           # library(rstan)
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)

# load edge weight model and data frames
load('motnp_bisonr_edgescalculated_strongprior.RData') # currently all saved versions of this use prior N(0,1.5) for fixed, but running regression at the moment with normal(0,1) instead

#### nodal regression -- mean age estimate only, not full distribution ####
df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
colnames(df_nodal) <- c('node_2_males','id_2')
df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
colnames(df_nodal) <- c('node','id')

motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
  select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')
motnp_ages$draw <- rep(1:8000, length(unique(motnp_ages$id)))

#mean_motnp_ages <- df_nodal
#mean_motnp_ages$age <- NA
#for(i in 1:nrow(mean_motnp_ages)){
#  x <- motnp_ages[motnp_ages$id == mean_motnp_ages$id[i],]
#  mean_motnp_ages$age[i] <- mean(x$age)
#  rm(x)
#}

#rm(list= ls()[!(ls() %in% c('mean_motnp_ages','motnp_edge_weights_strongpriors','motnp_ages')) ]) ; gc()

# define PDF output
#pdf('../outputs/motnp_nodalregression_plots_meanage.pdf')

## eigenvector only ####
# prior predictive check
priors$fixed
prior_check(priors, 'binary')

age <- 1:60
plot(NULL, xlim = c(0,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  intercept <- rbeta(1,1,1)
  beta <- rnorm(1, 0, 0.005)
  lines(x = age, y = intercept + age*beta, col = 'blue')
  lines(x = age, y = intercept + age*plogis(beta), col = 'red') # can't ever be a negative effect so that can't be right
}

mean_age <- mean(motnp_ages$age)
plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
     xlab = 'age', ylab = 'eigenvector centrality')
for(i in 1:100){
  intercept <- rbeta(1,2,2)   # this isn't right but I don't think there is a prior for the intercept?? I've gone for a symmetrical one here that in itself explores most of the parameter space and allows it to see whether some of the lines are steep enough to go from top to bottom, but on the assumption that when combined, they will explore only the space relevant to their starting position
  beta <- rnorm(1, 0, 0.005)
  lines(x = age, y = intercept + (age - mean_age)*beta, col = rgb(0,0,1,0.5)) # vast majority come out somewhere sensible, and those that don't would if they started at a different value for age 10 so that comes down to my ability to work out what the intercept prior should actually be -- think this is a good prior for the slope
}

#prior_eigen <- bison_brm_get_prior(
#  bison(node_eigen(node)) ~ age,
#  list(motnp_edge_weights_strongpriors),
#  counts_df_model
#)
#prior_eigen$fixed <- 'normal(0,1)'

#mean_motnp_eigen <- bison_brm(
#  bison(node_eigen(node)) ~ age,
#  motnp_edge_weights_strongpriors,
#  mean_motnp_ages,
#  chains = 4,
#  iter = 10000,
#  thin = 2
#)
#summary(mean_motnp_eigen) # FIT 100 IMPUTED MODELS, 4 CHAINS PER MODEL, EACH 1000 DRAWS LONG (+1000 DRAWS WARM UP). WARNING AT END OF MODEL RUN THAT CHAINS <3 DRAWS LONG AS ACTUAL CHAIN IS ONLY 1 WARMUP AND 1 SAMPLE. ONLY IMPUTED CHAINS FOLLOW THE SPECIFIED ITERATIONS AND THINNING

#mean_eigen_values <- mean_motnp_eigen$data
#plot(mean_eigen_values$bison_node_eigen ~ mean_eigen_values$age, 
#     las = 1, pch = 19, col = rgb(0,0,1,0.2),
#     xlab = 'mean age estimate', ylab = 'eigenvector centrality',
#     main = 'no effect of age on eigenvector centrality')

#mean_eigen_summary <- mean_motnp_eigen$fit

#hist(mean_motnp_eigen$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

#post_eigen <- as_draws_df(mean_motnp_eigen)
#hist(post_eigen$b_age) # what scale is this on? should it go through plogis() (aka inverse logit) or not -- what scale does eigenvector read as?
#plot(data = post_eigen[post_eigen$.chain == 1,], b_age ~ .draw, type = 'l')

## save output
#save.image('motnp_nodalregression_meanage.RData')

## strength only ####
#prior_strength <- bison_brm_get_prior(
#  bison(node_eigen(node)) ~ age,
#  list(motnp_edge_weights_strongpriors),
#  counts_df_model
#)
#prior_strength$fixed <- 'normal(0,1)'

#mean_motnp_strength <- bison_brm(
#  bison(node_strength(node)) ~ age,
#  motnp_edge_weights_strongpriors,
#  mean_motnp_ages,
#  chains = 4,
#  iter = 10000,
#  thin = 2
#)

#summary(mean_motnp_strength) # FIT 100 IMPUTED MODELS, 4 CHAINS PER MODEL, EACH 1000 DRAWS LONG (+1000 DRAWS WARM UP). WARNING AT END OF MODEL RUN THAT CHAINS <3 DRAWS LONG AS ACTUAL CHAIN IS ONLY 1 WARMUP AND 1 SAMPLE. ONLY IMPUTED CHAINS FOLLOW THE SPECIFIED ITERATIONS AND THINNING

#mean_strength_values <- mean_motnp_strength$data
#plot(mean_strength_values$bison_node_strength ~ mean_strength_values$age, 
#     las = 1, pch = 19, col = rgb(0,0,1,0.2),
#     xlab = 'mean age estimate', ylab = 'strength',
#     main = 'no effect of age on node strength')

#mean_strength_summary <- mean_motnp_strength$fit

#hist(mean_motnp_strength$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

#post_strength <- as_draws_df(mean_motnp_strength)
#hist(post_strength$b_age) # what scale is this on? should it go through plogis() (aka inverse logit) or not -- what scale does strength read as?
#plot(data = post_strength[post_strength$.chain == 1,], b_age ~ .draw, type = 'l')

## save output
#save.image('motnp_nodalregression_meanage.RData')

#### full age distribution ####
## eigenvector only ####
pdf('../outputs/motnp_nodalregression_eigen_agedistribution.pdf')

## define priors
#prior <- bison_brm_get_prior(
#  bison(node_eigen(node)) ~ age,
#  list(motnp_edge_weights_strongpriors),
#  motnp_ages
#)
#prior$prior[1] <- "normal(0,10)"
#bison_brm(
#  bison(global_cv(bison_network)) ~ bison_network,
#  list(fit_edge_1, fit_edge_2),
#  df_global,
#  prior=prior
#)

## define priors
motnp_edge_weights_strongpriors$model_data$prior_fixed_sigma <- 0.005

## reduce age data to manageable number of draws
test_ages <- motnp_ages[motnp_ages$draw %in% sample(motnp_ages$draw, 2000, replace = F),]

## run model
motnp_eigen <- bison_brm(
  bison(node_eigen(node)) ~ age,
  motnp_edge_weights_strongpriors,
  #motnp_ages,
  test_ages,
  chains = 4,
  iter = 1000,
  cores = 4
)
summary(motnp_eigen)

eigen_values <- motnp_eigen$data
plot(eigen_values$bison_node_eigen ~ eigen_values$age, las = 1, pch = 19, col = rgb(0,0,1,0.2),
     xlab = 'mean age estimate', ylab = 'eigenvector centrality',
     main = 'effect of age on eigenvector centrality')

mod_summary <- motnp_eigen$fit

hist(motnp_eigen$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

# compare to null model
motnp_eigen_null <- bison_brm(
  bison(node_eigen(node)) ~ 1,
  motnp_edge_weights_strongpriors,
  #motnp_ages,
  test_ages,
  chains = 4,
  cores = 4
)
#model_comparison(list(non_random_model = motnp_eigen, random_model = motnp_eigen_null))
#model_averaging(list(non_random_model = motnp_eigen, random_model = motnp_eigen_null))

save.image('motnp_nodalregression_eigenvector.RData')

rm(motnp_eigen, eigen_values, mod_summary, motnp_eigen_null) ; gc()
dev.off()

## strength only ####
#pdf('../outputs/motnp_nodalregression_strength_agedistribution.pdf')

#motnp_strength <- bison_brm(
#  bison(node_strength(node)) ~ age,
#  motnp_edge_weights_strongpriors,
#  motnp_ages,
#  chains = 4,
#  iter = 10000
#)
#summary(motnp_strength)

#strength_values <- motnp_strength$data
#plot(strength_values$bison_node_strength ~ strength_values$age, las = 1, pch = 19, col = rgb(0,0,1,0.2),
#     xlab = 'mean age estimate', ylab = 'strength centrality',
#     main = 'effect of age on strength centrality')

#mod_summary <- motnp_strength$fit

#hist(motnp_strength$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

# compare to null model
#motnp_strength_null <- bison_brm(
#  bison(node_strength(node)) ~ 1,
#  motnp_edge_weights,
#  motnp_ages,
#  chains = 4
#)
#model_comparison(list(non_random_model = motnp_strength, random_model = motnp_strength_null))

#save.image('motnp_nodalregression_strength.RData')

#dev.off()
