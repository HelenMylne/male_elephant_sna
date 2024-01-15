#### information ####
# script takes data input from edge weight estimation for MPNP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
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

load('mpnp_edgecalculations/mpnplong_edgeweights_conditionalprior.RData')
rm(edge_binary, fit_edges_mpnp1, edgelist, eles, missing_age, mean_ages, summary, make_edgelist, plot_network_threshold_mpnp) ; gc()

pdf('../outputs/mpnp_dyadicregression.pdf')
theme_set(theme_classic())

### filter down to only elephants with known age categories ####
## ages
hist(nodes$age, breaks = 50)
nodes$age_cat <- ifelse(nodes$age < 10, 1,
                        ifelse(nodes$age < 16, 2,
                               ifelse(nodes$age < 21, 3,
                                      ifelse(nodes$age < 26, 4,
                                             ifelse(nodes$age < 36, 5, 6)))))
## nodes data frame
nodes_to_include <- nodes %>% 
  filter(! is.na(age_cat))
n_age_cat <- length(unique(nodes_to_include$age_cat))
ele_ids <- nodes_to_include$id
n_nodes <- length(ele_ids)

## dyads data frame
dyads_to_include <- counts_df %>% 
  filter(id_1 %in% ele_ids) %>% 
  filter(id_2 %in% ele_ids)
n_dyads <- nrow(dyads_to_include)

#### identify younger and older ####
## categorise ages
node_join <- nodes_to_include %>% 
  mutate(node_rank = as.integer(as.factor(node))) %>% 
  dplyr::select(node_rank, node, id, age_cat) %>% 
  rename(node_1 = node, id_1 = id, node_rank_1 = node_rank, age_cat_1 = age_cat) %>% 
  mutate(node_2 = node_1, id_2 = id_1, node_rank_2 = node_rank_1, age_cat_2 = age_cat_1)
dyads_to_include <- dyads_to_include %>% 
  left_join(node_join[,c('node_rank_1','node_1','id_1','age_cat_1')], by = c('node_1', 'id_1')) %>% 
  left_join(node_join[,c('node_rank_2','node_2','id_2','age_cat_2')], by = c('node_2', 'id_2'))

## identify younger and older
dyads_to_include$age_min <- NA ; dyads_to_include$age_max <- NA
for(i in 1:nrow(dyads_to_include)){
  dyads_to_include$age_min[i] <- min(dyads_to_include$age_cat_1[i], dyads_to_include$age_cat_2[i])
  dyads_to_include$age_max[i] <- max(dyads_to_include$age_cat_1[i], dyads_to_include$age_cat_2[i])
}

## clean up dyads_to_include a bit
dyads_to_include <- dyads_to_include %>% 
  dplyr::select(dyad_id, id_1, id_2, node_1, node_2, node_rank_1, node_rank_2, age_cat_1, age_cat_2, age_min, age_max)

#### summarise edges ####
edge_samples_all <- edge_samples
edge_samples <- edge_samples[, dyads_to_include$dyad_id]
dyads_to_include$mean_edge <- apply(edge_samples, 2, mean)

#### plot raw data ####
ggplot()+
  geom_jitter(data = dyads_to_include,
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
  left_join(dyads_to_include, by = 'dyad_id')
ggplot()+
  geom_point(data = edges,
             aes(x = age_min, y = edge_draw,
                 colour = as.factor(age_max)),
             shape = 19, alpha = 0.05)+
  geom_point(data = dyads_to_include, aes(x = age_min, y = mean_edge),
             shape = 19, colour = 'white', size = 1)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

## heat map: x = min age, y = max age, colour = edge weight
ggplot(dyads_to_include)+
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
n_age_cat <- length(unique(dyads_to_include$age_max))
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
  age_min_cat = dyads_to_include$age_min,   # age of younger dyad member
  age_max_cat = dyads_to_include$age_max,   # age of  older  dyad member
  node_1 = dyads_to_include$node_1,         # node IDs for multimembership effects
  node_2 = dyads_to_include$node_2,         # node IDs for multimembership effects
  prior_min = rep(1, n_age_cat),            # prior for minimum age slope
  prior_max = rep(1, n_age_cat)             # prior for maximum age slope
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression_motnp.stan')
n_chains <- 4
n_samples <- 1000

## fit dyadic regression
fit_dyadreg_mpnp <- dyadic_regression$sample(
  data = dyad_data,
  iter_warmup = n_samples,
  iter_sampling = n_samples,
  chains = n_chains,
  parallel_chains = n_chains)

## save output
save.image('mpnp_dyadicregression.RData')
dev.off()

#### check outputs ####

#### plot predictions ####

