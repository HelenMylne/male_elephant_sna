#### information #####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

library(cmdstanr, lib.loc = '../packages/')        # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/')       # library(tidyverse)
library(brms, lib.loc = '../packages/')            # library(brms)
library(Rcpp, lib.loc = '../packages/')            # library(Rcpp)
library(ggdist, lib.loc = '../packages/')          # library(ggdist)
library(posterior, lib.loc = '../packages/')       # library(posterior)
library(bayesplot, lib.loc = '../packages/')       # library(bayesplot)
#library(rstan, lib.loc = '../packages/')           # library(rstan)
library(igraph, lib.loc = '../packages/')          # library(igraph)
library(LaplacesDemon, lib.loc = '../packages/')   # library(LaplacesDemon)
library(bisonR, lib.loc = '../packages/')          # library(bisonR)
library(extraDistr, lib.loc = '../packages/')      # library(extraDistr)

# load edge weight model and data frames
load('motnp_bisonr_edgescalculated_strongprior.RData') # currently all saved versions of this use prior N(0,1.5) for fixed, but running regression at the moment with normal(0,1) instead
rm(motnp_edges_null_strongpriors) ; gc()

#### create nodal regression data ####
df_nodal <- distinct(counts_df[,c('node_1_males','id_1')])
colnames(df_nodal) <- c('node_2_males','id_2')
df_nodal <- rbind(df_nodal, counts_df[nrow(counts_df),c('node_2_males','id_2')])
colnames(df_nodal) <- c('node','id')

mean_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
  select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  apply(MARGIN = 2, FUN = mean) %>% 
  as.data.frame()
mean_ages$mean_age <- mean_ages$.
mean_ages <- mean_ages %>% 
  mutate(id = sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  select(id, mean_age) %>% 
  left_join(df_nodal, by = 'id')

motnp_ages <- readRDS('../data_processed/motnp_ageestimates_mcmcoutput.rds') %>% 
  select(sort(unique(c(counts_df$id_1, counts_df$id_2)))) %>% 
  pivot_longer(cols = everything(), names_to = 'id', values_to = 'age')
motnp_ages <- left_join(motnp_ages, df_nodal, by = 'id')
motnp_ages$draw <- rep(1:8000, length(unique(motnp_ages$id)))

#### set up priors ####
pdf('../outputs/motnp_nodalregression_strength_agedistribution.pdf')
strengths <- extract_metric(motnp_edge_weights_strongpriors, 'node_strength')
hist(strengths, breaks = 100)       # nearly normal? maybe Weibull distributed
mean_strength <- apply(X = strengths, MARGIN = 2, FUN = mean)
hist(mean_strength, breaks = 20)    # Weibull distributed

#N <- 1000
#plot(NULL, las = 1, xlim = c(0,250), ylim = c(0,0.1),
#     xlab = 'strength', ylab = 'density')
#lines(density(rweibull(N, shape = 5, scale = 150)), col = 'red')
#lines(density(rweibull(N, shape = 3, scale = 150)), col = 'blue')
#lines(density(rweibull(N, shape = 3, scale = 100)), col = 'purple')
#lines(density(rweibull(N, shape = 7, scale = 150)), col = 'green')
#plot(NULL, las = 1, xlim = c(0,25), ylim = c(0,0.2),
#     xlab = 'strength', ylab = 'density')
#lines(density(rweibull(N, shape = 7, scale = 15)), col = 'red', lty = 2)
#lines(density(rweibull(N, shape = 5, scale = 15)), col = 'purple', lty = 2)
#lines(density(rweibull(N, shape = 7, scale = 20)), col = 'blue', lty = 2)
#lines(density(rweibull(N, shape = 5, scale = 20)), col = 'turquoise', lty = 2)

### prior predictive check
motnp_edge_weights_strongpriors$model_data$prior_fixed_mu
motnp_edge_weights_strongpriors$model_data$prior_fixed_sigma

## plot raw data (mean values only)
mean_ages$mean_strength <- mean_strength
plot(mean_ages$mean_strength ~ mean_ages$mean_age, las = 1,
     pch = 16, col = rgb(0,0,1,0.3),
     xlab = 'mean age estimate', ylab = 'mean strength estimate')
motnp_ages <- left_join(motnp_ages, mean_ages[,1:3], by = c('id','node'))
length(which(is.na(motnp_ages$mean_age) == TRUE))

## set up prior predictive check
N <- 100
x_bar <- mean(mean_ages$mean_age)
a <- rweibull(N, shape = 5, scale = 12)
b_mu <- 0
b_sigma <- 0.15
b <- rnorm(N, b_mu, b_sigma)

## plot prior predictive
plot(NULL, las = 1, xlim = c(10,60), ylim = c(0,25),
     xlab = 'age', ylab = 'strength')
abline(h = c(5,18), lty = 2)
for(i in 1:N){
  curve(a[i] + b[i]*(x - x_bar),
        from = 10, to = 60, add = T, col = rgb(0,0,1,0.4))
}

motnp_edge_weights_strongpriors$model_data$prior_fixed_mu <- b_mu
motnp_edge_weights_strongpriors$model_data$prior_fixed_sigma <- b_sigma

#### fit strength model ####
motnp_strength <- bison_brm(
  #bison(node_strength(node)) ~ age,
  bison(node_strength(node)) ~ mean_age,
  motnp_edge_weights_strongpriors,
  motnp_ages,
  chains = 4,
  iter = 1000
)
summary(motnp_strength)

hist(motnp_strength$rhats[,2], las = 1, main = 'Rhat values for 100 imputed model runs', xlab = 'Rhat')

# compare to null model
#motnp_strength_null <- bison_brm(
#  bison(node_strength(node)) ~ 1,
#  motnp_edge_weights,
#  motnp_ages,
#  chains = 4
#)
#model_comparison(list(non_random_model = motnp_strength, random_model = motnp_strength_null))

save.image('motnp_nodalregression_strength.RData')

dev.off()
