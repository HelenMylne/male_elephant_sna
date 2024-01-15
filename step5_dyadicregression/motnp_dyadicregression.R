#### information ####
# script takes data input from edge weight estimation for MOTNP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### set up ####
# library(cmdstanr) ; library(tidyverse) ; library(car) ; library(LaplacesDemon)
library(cmdstanr, lib.loc = '../packages/')      # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/')     # library(tidyverse)
library(car, lib.loc = '../packages/')           # library(car)
library(LaplacesDemon, lib.loc = '../packages/') # library(LaplacesDemon)
#library(bisonR, lib.loc = '../packages/')        # library(bisonR)
#library(brms, lib.loc = '../packages/')          # library(brms)
#library(rstan, lib.loc = '../packages/')         # library(rstan)
#library(Rcpp, lib.loc = '../packages/')          # library(Rcpp)

set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0')

load('motnp_edgeweights_conditionalprior.RData')
rm(edge_binary, fit_edges_motnp, edgelist, edges, motnp_ages, x, summary, i, make_edgelist, plot_network_threshold) ; gc()

pdf('../outputs/motnp_dyadicregression.pdf')
theme_set(theme_classic())

#### identify younger and older ####
## categorise ages
counts_df$age_cat_num_1 <- as.numeric(counts_df$age_cat_id_1)-2
table(counts_df$age_cat_num_1, counts_df$age_category_1)
counts_df$age_cat_num_2 <- as.numeric(counts_df$age_cat_id_2)-2
table(counts_df$age_cat_num_2, counts_df$age_category_2)

## identify younger and older
counts_df$age_min <- NA ; counts_df$age_max <- NA
for(i in 1:nrow(counts_df)){
  counts_df$age_min[i] <- min(counts_df$age_cat_num_1[i], counts_df$age_cat_num_2[i])
  counts_df$age_max[i] <- max(counts_df$age_cat_num_1[i], counts_df$age_cat_num_2[i])
}

## clean up counts_df a bit
counts_df <- counts_df %>% 
  dplyr::select(dyad_id, id_1, id_2, node_1, node_2, node_1_males, node_2_males, age_category_1, age_category_2, age_cat_num_1, age_cat_num_2, age_min, age_max)

#### summarise edges ####
counts_df$mean_edge <- apply(edge_samples, 2, mean)

#### plot raw data ####
ggplot()+
  geom_jitter(data = counts_df,
              aes(x = age_min, y = mean_edge,
                  colour = as.factor(age_max)),
              shape = 19, alpha = 0.2)+
  scale_x_continuous('age category of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  labs(colour = 'age category \nof older \ndyad member')

edges <- edge_samples %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'edge_draw') %>% 
  mutate(dyad_id = as.numeric(dyad_id)) %>% 
  left_join(counts_df, by = 'dyad_id')
ggplot()+
  geom_point(data = edges,
             aes(x = age_min, y = edge_draw,
                 colour = as.factor(age_max)),
             shape = 19, alpha = 0.05)+
  geom_point(data = counts_df, aes(x = age_min, y = mean_edge),
             shape = 19, colour = 'white', size = 1)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

## heat map: x = min age, y = max age, colour = edge weight
ggplot(counts_df)+
  geom_tile(aes(x = age_min, y = age_max,
                fill = mean_edge))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age of older dyad member')+
  labs(fill = 'logit(mean edge weight)')

#### fit multivariate Gaussian distribution to output of age model ####
### quantities will be given to Stan model as data to model joint posteriors of edge weight in the regression
logit_edge_draws <- logit(edge_samples)
logit_edge_draws_mu <- apply(logit_edge_draws, 2, mean)
logit_edge_draws_cov <- cov(logit_edge_draws)

#### plot to check how well multivariate approximation is working ####
### Randomly selecting samples to examine
num_check <- 50
selected_samples <- sample(1:(n_dyads-1), num_check, replace = FALSE)

### Setting grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

### plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values <- rnorm(1e5, mean = mu, sd = sd)
  
  hist(unlist(logit_edge_draws[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="blue", lw=1.5)
}

for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values <- rnorm(1e5, mean = mu, sd = sd)
  
  plot(unlist(logit_edge_draws[,i]), unlist(logit_edge_draws[,i+1]),
       col = rgb(0,0,1,0.05), las = 1, main = paste("cov ", i ,"&",i+1))
}

### reset plot window and clean up
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))
rm(cols, fitted_values, i, j, num_check, rows, sd, selected_samples) ; gc()

#### prior predictive check ####
n <- 100
n_age_cat <- length(unique(counts_df$age_max))
beta_age_min <- rnorm(n, 0, 1.5)
beta_age_max <- rnorm(n, 0, 1.5)
intercept <- rnorm(n, 0, 1)
min_dirichlet <- rdirichlet(n, rep(1, n_age_cat))
max_dirichlet <- rdirichlet(n, rep(1, n_age_cat))
plot(NULL, las = 1, xlab = 'minimum age (standardised)', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5), xlim = c(1,n_age_cat))
abline(h = min(logit_edge_draws_mu), lty = 2) ; abline(h = max(logit_edge_draws_mu), lty = 2)
x <- 1:n_age_cat
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age_min[i]*sum(min_dirichlet[i,][1:x[j]]) + beta_age_max[i]*sum(max_dirichlet[i,][1:x[j]])
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age_max, beta_age_min, intercept, min_dirichlet, max_dirichlet, i, y, j, x) ; gc()

#### fit dyadic regression ####
## create data list
dyad_data <- list(
  num_dyads = n_dyads,                      # number of dyads
  num_nodes = n_nodes,                      # number of nodes
  num_age_cat = n_age_cat,                  # number of unique age categories
  length_dirichlet = n_age_cat + 1,         # number of unique age categories + 1
  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
  age_min_cat = sim_dyads$age_min,          # age of younger dyad member
  age_max_cat = sim_dyads$age_max,          # age of  older  dyad member
  node_1 = sim_dyads$node_1,                # node IDs for multimembership effects
  node_2 = sim_dyads$node_2,                # node IDs for multimembership effects
  prior_min = rep(1, n_age_cat),            # prior for minimum age slope
  prior_max = rep(1, n_age_cat)             # prior for maximum age slope
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression_motnp.stan')
n_chains <- 4
n_samples <- 1000

## fit dyadic regression
fit_dyadreg_motnp <- dyadic_regression$sample(
  data = dyad_data,
  iter_warmup = n_samples,
  iter_sampling = n_samples,
  chains = n_chains,
  parallel_chains = n_chains)

## save output
save.image('motnp_dyadicregression.RData')
dev.off()

#### check outputs ####

#### plot predictions ####

