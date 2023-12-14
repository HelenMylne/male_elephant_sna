#### information ####
# script takes data input from edge weight estimation for ANP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### set up ####
#library(tidyverse) ; library(LaplacesDemon) ; library(car) ; library(cmdstanr) ; library(bisonR) ; library(brms)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(car, lib.loc = '../packages/')       # library(car)
#library(bisonR, lib.loc = '../packages/')    # library(bisonR)
#library(brms, lib.loc = '../packages/')      # library(brms)
library(LaplacesDemon, lib.loc = '../packages/')

set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan2/cmdstan-2.33.1/')

load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
# rm(edgelist, edge_binary, nodes) ; gc()

pdf('../outputs/anp1_dyadicregression_plots.pdf')

#### prior predictive check ####
age_min <- 5:50
age_max <- seq(10, 50, length.out = 12)
par(mfrow = c(4,3))
for(age in age_max){
  plot(NULL, xlim = c(5,50), ylim = c(0,1), las = 1,
       xlab = 'age of younger', ylab = 'edge weight', main = paste0('oldest = ',round(age)))
  for(i in 1:100){
    beta_min <- rnorm(1, 0, 0.1)
    beta_max <- rnorm(1, 0, 0.1)
    #beta_int <- rnorm(1, 0, 0.005)
    min_new <- age_min[age > age_min]
    lines(x = min_new, y = plogis(min_new*beta_min +
                                    #min_new*age*beta_int +
                                    age*beta_max),
          col = rgb(0,0,1,0.2))
  }
}
rm(age, age_min, age_max, beta_max, beta_min, i, min_new) ; gc()
par(mfrow = c(1,1))

# add time marker
print(paste0('prior predictive check completed at ', Sys.time()))

#### fit multivariate Gaussian distribution to output of edge weight model ####
# To parameterise the multivariate normal approximation, we use the sample mean and covariance matrix, calculated from the posterior edge weight samples using the following code:
### get the weights on the logit scale (they are not currently because we used a beta prior and identity link here rather than logistic link)
logit_weights <- apply(edge_samples, 2, qlogis)

### fit a multivariate normal dist to the edges -- quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
logit_edge_draws_mu <- apply(logit_weights, 2, mean)
logit_edge_draws_cov <- cov(logit_weights)

#### plot to see how well the approximation is working ####
### Randomly selecting samples to examine
num_check <- 20
selected_samples <- sample(1:n_dyads, num_check, replace = FALSE)

### Setting grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

### plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values_logit <- rnorm(1e5, mean = mu, sd = sd)
  fitted_values_original <- plogis(fitted_values_logit)
  
  hist(unlist(edge_samples[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values_original), col="blue", lw=1.5)
}

for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values_logit <- rnorm(1e5, mean = mu, sd = sd)
  fitted_values_original <- plogis(fitted_values_logit)
  
  plot(unlist(edge_samples[,i]), unlist(edge_samples[,i+1]), col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

### reset plot window
par(mfrow=c(1,1))

# save image so far
save.image('anpshort1_dyadicregression_minmax.RData')

# add time marker
print(paste0('multivariate Gaussian approximation fitted at ', Sys.time()))

#### fit dyadic regression ####
### identify older and younger of dyad
cdf_1$age_min <- NA ; cdf_1$age_max <- NA
for(i in 1:nrow(cdf_1)){
  x <- c(cdf_1$age_start_1[i],cdf_1$age_start_2[i])
  cdf_1$age_min[i] <- min(x)
  cdf_1$age_max[i] <- max(x)
}

## convert node IDs to values 1:52 not numbers based on casename
all_node_IDs <- unique(c(cdf_1$id_1, cdf_1$id_2))
n_nodes <- length(all_node_IDs)
cdf_1 <- cdf_1 %>% 
  rename(node_1_original = node_1,
         node_2_original = node_2) %>% 
  mutate(node_1_period = as.integer(factor(id_1, levels = all_node_IDs)),
         node_2_period = as.integer(factor(id_2, levels = all_node_IDs)),)

## create data list
dyad_data <- list(
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
  age_min = cdf_1$age_min,                  # age of younger dyad member
  age_max = cdf_1$age_max,                  # age of  older  dyad member
  #age_diff = cdf_1$age_diff,                # age difference between dyad members
  node_1 = cdf_1$node_1_period,             # node IDs for multimembership effects
  node_2 = cdf_1$node_2_period              # node IDs for multimembership effects
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression.stan')

# add time marker
print(paste0('start model run at ', Sys.time()))

## fit dyadic regression
fit_dyadreg_anp1 <- dyadic_regression$sample(
  data = dyad_data,
  chains = n_chains,
  parallel_chains = n_chains)

# save image so far
save.image('anpshort1_dyadicregression_minmax.RData')

# add time marker
print(paste0('finish model run at ', Sys.time()))

