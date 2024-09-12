## script for working out the problems with cov() in dyadic regressions ##
#### set up ####
library(tidyverse)
library(rstan)
library(LaplacesDemon)

#### MOTNP: import data ####
# load('motnp_edgeweights_conditionalprior.RData') ; write_csv(counts_df, '../data_processed/step3_edgeweightestimation/motnp_countsdf.csv') ; saveRDS(edge_samples,'../data_processed/step3_edgeweightestimation/motnp_edgeweightmatrix_conditionalprior.RDS')
dyads <- read_csv('../data_processed/step3_edgeweightestimation/motnp_countsdf.csv')
nodes <- read_csv('../data_processed/step1_dataprocessing/motnp_elenodes.csv') %>% 
  filter(sex == 'M' & age_category %in% c('10-15','15-19','20-25','25-40','40+')) %>% 
  filter(id_no != 'M0051') %>% filter(id_no != 'M0113')
edges <- readRDS('../data_processed/step3_edgeweightestimation/motnp_edgeweightmatrix_conditionalprior.RDS')

## identify younger and older member of each dyad
nodes$age_cat_num <- ifelse(nodes$age_category == '10-15', 1,
                            ifelse(nodes$age_category == '15-19', 2,
                                   ifelse(nodes$age_category == '20-25', 3,
                                          ifelse(nodes$age_category == '25-40', 4,
                                                 5))))
dyads <- dyads %>% mutate(age_cat_min = NA, age_cat_max = NA)
for(i in 1:nrow(dyads)){
  dyads$age_cat_min[i] <- min(dyads$age_cat_id_1[i],dyads$age_cat_id_2[i])
  dyads$age_cat_max[i] <- max(dyads$age_cat_id_1[i],dyads$age_cat_id_2[i])
}

#### MOTNP: calculate multivariate Gaussian distribution ####
logit_edges <- logit(edges)
logit_edge_mu <- apply(logit_edges, 2, mean)
logit_edge_cov <- cov(logit_edges)

## randomly select samples to examine
num_check <- 25
selected_samples <- sample(1:(n_dyads-1), num_check, replace = FALSE)

## set grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

## plot
for (i in selected_samples) {
  mu <- logit_edge_mu[i]
  sd <- sqrt(logit_edge_cov[i,i])
  
  fitted_values <- rnorm(1e5, mean = mu, sd = sd)
  
  hist(unlist(sim_dat[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="blue", lw=1.5)
}

for (i in selected_samples) {
  mu <- logit_edge_mu[i]
  sd <- sqrt(logit_edge_cov[i,i])
  
  fitted_values <- rnorm(1e5, mean = mu, sd = sd)
  
  plot(unlist(sim_dat[,i]), unlist(sim_dat[,i+1]), col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

## reset plot window
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))

#### MOTNP: fit model ####
## create data list
n_age_cat <- length(unique(nodes$age_category))
dyad_data <- list(
  num_dyads = nrow(dyads),                            # number of dyads
  num_nodes = nrow(nodes),                            # number of nodes
  num_age_cat = n_age_cat,            # number of unique age categories
  length_dirichlet = n_age_cat + 1,   # number of unique age categories + 1
  logit_edge_mu = logit_edge_mu,      # means of the logit edge weights
  logit_edge_cov = logit_edge_cov,    # standard deviation of logit edge weights
  age_min_cat = dyads$age_cat_min,    # age of younger dyad member
  age_max_cat = dyads$age_cat_max,    # age of younger dyad member
  node_1 = dyads$node_1,              # node IDs for multimembership effects
  node_2 = dyads$node_2,              # node IDs for multimembership effects
  prior_min = rep(1, n_age_cat),      # prior for minimum age slope
  prior_max = rep(1, n_age_cat)       # prior for maximum age slope
)

## load dyadic regression model
dyadic_regression <- stan_model('models/old_regression_models/dyadic_regression_motnp.stan')

#### ANP: import data ####
load('anp_edgecalculations/excluding_under10s/') ; write_csv(counts_df, '../data_processed/step3_edgeweightestimation/') ; saveRDS(edge_samples,'../data_processed/step3_edgeweightestimation/')
dyads <- read_csv('../data_processed/step3_edgeweightestimation/')
nodes <- read_csv('../data_processed/step1_dataprocessing/')
edges <- readRDS('../data_processed/step3_edgeweightestimation/')

## identify younger and older member of each dyad
nodes$age_cat_num <- ifelse(nodes$age_category == '10-15', 1,
                            ifelse(nodes$age_category == '15-19', 2,
                                   ifelse(nodes$age_category == '20-25', 3,
                                          ifelse(nodes$age_category == '25-40', 4,
                                                 5))))
dyads <- dyads %>% mutate(age_cat_min = NA, age_cat_max = NA)
for(i in 1:nrow(dyads)){
  dyads$age_cat_min[i] <- min(dyads$age_cat_id_1[i],dyads$age_cat_id_2[i])
  dyads$age_cat_max[i] <- max(dyads$age_cat_id_1[i],dyads$age_cat_id_2[i])
}



