#### information ####
# script to test dyadic regression with age distribution

#### set up ####
#library(tidyverse) ; library(car) ; library(cmdstanr) ; library(bisonR) ; library(brms)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(car, lib.loc = '../packages/')       # library(car)
#library(bisonR, lib.loc = '../packages/')    # library(bisonR)
#library(brms, lib.loc = '../packages/')      # library(brms)

#set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan/cmdstan-2.31.0')

load('motnp_edgeweights_conditionalprior.RData')
# rm(edgelist, edge_binary, nodes) ; gc()

pdf('../outputs/motnp_dyadicregression_plots.pdf')

#### prior predictive check ####
age_min <- 5:50
age_diff <- c(0,2,5,10,15,20,25,30)
par(mfrow = c(4,2), mai = c(0.3,0.3,0.3,0.1))
for(age in age_diff){
  plot(NULL, xlim = range(age_min), ylim = c(0,1), las = 1,
       xlab = 'age of younger', ylab = 'edge weight',
       main = paste0('diff = ',age))
  for(i in 1:100){
    beta_min <- rnorm(1, 0, 0.1)
    beta_diff <- rnorm(1, 0, 0.1)
    lines(x = age_min, y = plogis(age_min*beta_min + 
                                    age*beta_diff), 
          col = rgb(0,0,1,0.2))
    lines(x = age_min, y = age_min*beta_min + age*beta_diff, 
          col = rgb(1,0,0,0.2))
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
save.image('motnp_dyadicregression.RData')

# add time marker
print(paste0('multivariate Gaussian approximation fitted at ', Sys.time()))

#### identify older and younger of dyad ####
# load('motnp_dyadicregression.RData')
load('motnp_ageestimation.RData')

# import age data
ages <- readRDS('../data_processed/step2_ageestimation/motnp_ageestimates_mcmcoutput.rds') %>% 
  select(unique(c(counts_df$id_1, counts_df$id_2)))
age_mat <- as.matrix(ages)

# create matrix of age contrasts
all_node_IDs <- colnames(ages)
contrasts <- array(NA, dim = c(ncol(ages), ncol(ages), nrow(ages)),
                   dimnames = list(all_node_IDs, all_node_IDs, NULL))
for(id1 in 1:length(all_node_IDs)) {
  for(id2 in 1:length(all_node_IDs)) {
    for(draw in 1:nrow(ages)) {
      if(id1 == id2){
        contrasts[id1, id2, draw] <- 0
      }
      if(id1 > id2){
        contrasts[id1, id2, draw] <- contrasts[id2, id1, draw]
      }
      contrasts[id1, id2, draw] <- as.numeric(age_mat[draw,id1] - age_mat[draw,id2])
    }
  }
}

# save image so far
counts_df <- summary
rm(age_mat, edgelist, motnp_ages, summary, x, cols, draw, i, id1, id2, mu, num_check, rows, sd, selected_samples, make_edgelist, plot_network_threshold)
save.image('motnp_dyadicregression.RData')

# calculate age distributions -- test quality of normal approximation
par(mfrow = c(5,3), mai = c(0.3,0.3,0.3,0.3))
#for( i in sample(colnames(ages), 15, replace = F)){
for( i in colnames(ages)){
  age_i <- as.matrix(ages[,i])
  hist(age_i, main = i, probability = TRUE, ylim = c(0,0.2))
  lines(density(rnorm(n = length(age_i),
              mean = mean(age_i), sd = sd(age_i))),
        col = 'blue')
}

par(mfrow = c(2,3))
example_nodes <- c('M103','M102','M10','M1','M11')
example_node_ages <- c('10-15','15-20','20-25','25-40','40+')
for( i in 1:length(example_nodes)){
  age_i <- as.matrix(ages[,example_nodes[i]])
  hist(age_i, main = example_node_ages[i], breaks = 20,
       probability = TRUE, ylim = c(0,0.2))
  lines(density(rnorm(n = length(age_i),
                      mean = mean(age_i), sd = sd(age_i))),
        col = 'blue')
}
# 10-15, 15-20 and 20-25 = normal approximation is bang on
# 25-40 and 40+ = normal approximation is a bit too tall in the middle but width is good

# calculate age distributions -- test quality of gompertz bathtub approximation
par(mfrow = c(5,3), mai = c(0.3,0.3,0.3,0.3))
#for( i in sample(colnames(ages), 15, replace = F)){
for( i in colnames(ages)){
  age_i <- as.matrix(ages[,i])
  hist(age_i, main = i, probability = TRUE, ylim = c(0,0.2))
  lines(density(rnorm(n = length(age_i),
                      mean = mean(age_i), sd = sd(age_i))),
        col = 'blue')
}
# 10-15, 15-20 and 20-25 = normal approximation is bang on
# 25-40 and 40+ = normal approximation is a bit too tall in the middle but width is good





# set values for normal approximation of each node's age
nodes$age_mean <- NA ; nodes$age_stdv <- NA
for( i in colnames(ages)){
  age_i <- as.matrix(ages[,i])
  nodes$age_mean[nodes$id == i] <- mean(age_i)
  nodes$age_stdv[nodes$id == i] <- sd(age_i)
}
rm(age_i) ; gc()

# set values for normal approximation of minimum age for each dyad
counts_df$age_min_mu <- NA ; counts_df$age_min_sd <- NA
for(i in 1:nrow(counts_df)){
  dyad <- nodes %>% filter(id %in% c(counts_df$id_1[i], counts_df$id_2[i]))
  counts_df$age_min_mu[i] <- min(dyad$age_mean)
  counts_df$age_min_sd[i] <- dyad$age_stdv[which(dyad$age_mean == min(dyad$age_mean))]
}
hist(counts_df$age_min_mu)

# calculate age difference distributions -- test quality of normal approximation
for( i in sample(counts_df$dyad_id, 150, replace = F) ){
  age_i <- contrasts[counts_df$id_1[counts_df$dyad_id == i],
                     counts_df$id_2[counts_df$dyad_id == i], ]
  hist(age_i, main = i, probability = TRUE)
  lines(density(rnorm(n = length(age_i),
                      mean = mean(age_i), sd = sd(age_i))),
        col = 'blue')
}
# often approximation puts slightly too much probability density over the central ages but generally a very good fit to the contrasts data
par(mfrow = c(1,1))

# set values for normal approximation of age difference for each dyad
counts_df$age_diff_mu <- NA ; counts_df$age_diff_sd <- NA
for(i in 1:nrow(counts_df)){
  dyad <- contrasts[counts_df$id_1[i], counts_df$id_2[i], ] # always samples from top triangle so no issue with taking inverse values half the time, but could be in either direction -- need to make sure differences are positive, so it will be the max-min one
  counts_df$age_diff_mu[i] <- abs(mean(dyad))  # absolute value to ensure right way around
  counts_df$age_diff_sd[i] <- sd(dyad)
}
hist(counts_df$age_diff_mu)

#### fit dyadic regression ####
load('motnp_dyadicregression.RData')

## create data list
#n_nodes <- length(all_node_IDs)
#dyad_data <- list(
#  num_dyads = n_dyads,                      # number of dyads
#  num_nodes = n_nodes,                      # number of nodes
#  logit_edge_mu = logit_edge_draws_mu,      # sample means of the logit edge weights
#  logit_edge_cov = logit_edge_draws_cov,    # sample covariance of logit edge weights
#  age_min_mu = counts_df$age_min_mu,        # age of younger dyad member
#  age_diff_mu = counts_df$age_diff_mu,      # age difference between dyad members
#  age_min_sd = counts_df$age_min_sd,        # age of younger dyad member
#  age_diff_sd = counts_df$age_diff_sd,      # age difference between dyad members
#  node_1 = counts_df$node_1_males,          # node IDs for multimembership effects
#  node_2 = counts_df$node_2_males           # node IDs for multimembership effects
#)

samples_nodes <- sample(all_node_IDs, 20, replace = F)
test_df <- counts_df %>%
  filter(id_1 %in% samples_nodes) %>% 
  filter(id_2 %in% samples_nodes) %>% 
  mutate(node_1_test = as.integer(as.factor(id_1)),
         node_2_test = (as.integer(as.factor(id_2))+1))
logit_test_mu <- logit_edge_draws_mu[names(logit_edge_draws_mu) %in% test_df$dyad_id]
logit_test_cov <- logit_edge_draws_cov[rownames(logit_edge_draws_cov) %in% test_df$dyad_id, 
                                       colnames(logit_edge_draws_cov) %in% test_df$dyad_id]
test_data <- list(
  num_dyads = nrow(test_df),          # number of dyads
  num_nodes = length(samples_nodes),  # number of nodes
  logit_edge_mu = logit_test_mu,      # sample means of the logit edge weights
  logit_edge_cov = logit_test_cov,    # sample covariance of logit edge weights
  age_min_mu = test_df$age_min_mu,    # age of younger dyad member
  age_min_sd = test_df$age_min_sd,    # age of younger dyad member
  age_diff_mu = test_df$age_diff_mu,  # age difference between dyad members
  age_diff_sd = test_df$age_diff_sd,  # age difference between dyad members
  node_1 = test_df$node_1_test,       # node IDs for multimembership effects
  node_2 = test_df$node_2_test        # node IDs for multimembership effects
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression_agedistributions.stan')

# add time marker
print(paste0('start model run at ', Sys.time()))

## fit dyadic regression
fit_dyadreg_motnp <- dyadic_regression$sample(
  #data = dyad_data,
  data = test_data,
  #chains = n_chains,
  #parallel_chains = n_chains)
  chains = 1)

# save image so far
save.image('motnp_dyadicregression.RData')

# add time marker
print(paste0('finish model run at ', Sys.time()))

#### check outputs ####
#load('motnp_dyadicregression.RData')
rm(counts_df, counts_df, dyad_data, edge_binary, edgelist, edges, i, x, make_edgelist, plot_network_threshold_motnp) ; gc()

# obtain summary
fit_dyadreg_motnp$summary()

## extract draws
draws <- fit_dyadreg_motnp$draws(format = 'df')

## extract dyadic regression slopes
b_diff <- draws[,'beta_age_diff']
b_min <- draws[,'beta_age_min']
sigma <- draws[,'sigma']

## extract overall age effects
dimnames(draws)[[2]]
draws <- draws[,4:length(draws[1,])]

## extract minimum ages
age_min <- draws[,1:test_data$num_dyads]
age_diff <- draws[,(1+test_data$num_dyads):(2*test_data$num_dyads)]
draws <- draws[,((2*test_data$num_dyads)+1):length(draws[1,])]

## extract multi-membership samples -- CHECK THIS IS IN THE RIGHT ORDER
mm_matrix <- draws[,1:test_data$num_nodes]
sigma <- draws[,(test_data$num_nodes+1)]
sigma_mm <- draws[,(test_data$num_nodes+2)]
draws <- draws[,(test_data$num_nodes+3):length(draws[1,])]

## extract covariance values
l_cov <- draws[,1:(test_data$num_dyads^2)]
draws <- draws[,(1+(test_data$num_dyads^2)):length(draws[1,])]

## extract predictor values
predictor <- draws[,1:test_data$num_nodes]
rm(draws) ; gc()

# add time marker
print(paste0('parameters extracted at ', Sys.time()))

## save parameter values
#extract_slopes <- function(draws, n_samples = 1000, n_chains = 4){
#  draws <- as.data.frame(draws) %>% 
#    pivot_longer(everything(), names_to = 'chain', values_to = 'slope_draw') %>% 
#    separate(chain, sep = c(1,2), remove = T, into = c('chain','.', 'parameter')) %>% 
#    select(-.)
#  draws$chain <- as.numeric(draws$chain)
#  draws$position <- rep(1:n_samples, each = n_chains)
#  draws$draw_id <- draws$position + (draws$chain-1)*n_samples
#  return(draws)
#}
#b_diff <- extract_slopes(b_diff)
#b_min <- extract_slopes(b_min)
#sigma <- extract_slopes(sigma)

make_slopes_df <- function(draws_df, n_chains = 4){
  n_samples <- nrow(draws_df)/n_chains
  n_cols <- ncol(draws_df)
  if(n_cols > 1){
    draws_df <- draws_df %>% 
      pivot_longer(cols = everything(),
                   names_to = 'parameter') %>% 
      mutate(type = 'FILL THIS IN MANUALLY')
  } else{
    parameter <- colnames(draws_df)
    colnames(draws_df) <- 'value'
    draws_df <- draws_df %>% 
      mutate(parameter = parameter,
             type = parameter) %>% 
      relocate(parameter)
  }
  draws_df <- draws_df %>% 
    mutate(chain = rep(rep(1:n_chains, each = n_samples),n_cols),
           position = rep(rep(1:n_samples, n_chains), each = n_cols),
           draw_id = position + (chain-1)*n_samples)
}
b_diff <- make_slopes_df(b_diff)
b_min <- make_slopes_df(b_min)
sigma <- make_slopes_df(sigma)
sigma_mm <- make_slopes_df(sigma_mm)
age_diff <- make_slopes_df(age_diff) %>% 
  mutate(type = 'age_diff')
age_min <- make_slopes_df(age_min) %>% 
  mutate(type = 'age_min')
parameters <- rbind(b_min, b_diff,
                    sigma, sigma_mm,
                    age_min, age_diff)

## save data 
saveRDS(parameters, '../data_processed/motnp_dyadicregression_slopeparameters.RDS')

# add time marker
print(paste0('parameters saved to file at ', Sys.time()))

## traceplots
ggplot(data = parameters[! parameters$type %in% c('age_diff', 'age_min'),])+
  geom_line(aes(x = position, y = value, colour = as.factor(chain)))+
  theme_classic()+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap( . ~ parameter , scales = 'free_y')

ggplot(data = parameters[parameters$parameter %in% c('age_diff[1]','age_diff[2]','age_diff[3]','age_diff[4]','age_diff[5]','age_diff[6]','age_diff[7]','age_diff[8]','age_diff[9]','age_diff[10]','age_diff[11]','age_diff[12]'),])+
  geom_line(aes(x = position, y = value, colour = as.factor(chain)))+
  theme_classic()+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap( . ~ parameter , scales = 'free_y')

ggplot(data = parameters[parameters$parameter %in% c('age_min[1]','age_min[2]','age_min[3]','age_min[4]','age_min[5]','age_min[6]','age_min[7]','age_min[8]','age_min[9]','age_min[10]','age_min[11]','age_min[12]'),])+
  geom_line(aes(x = position, y = value, colour = as.factor(chain)))+
  theme_classic()+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap( . ~ parameter , scales = 'free_y')

# add time marker
print(paste0('traceplots run at ', Sys.time()))

# clean workspace
rm(age_diff, age_i, age_min, ages, all_node_IDs, b_diff, b_min, contrasts, dyad, dyadic_regression, edge_samples, fit_edges_motnp, fitted_values_logit, fitted_values_original, i, l_cov,logit_weights, mm_matrix, samples_nodes, sigma, sigma_mm) ; gc()

#### plot edges against age values ####
#edges <- readRDS('../data_processed/motnp_edgedistributions_conditionalprior.RDS')

### create mean data frame to plot average values over full distribution
mean_edges <- data.frame(dyad = test_df$dyad_id,
                         node_1 = test_df$node_1_test,
                         node_2 = test_df$node_2_test,
                         together = test_df$event_count,
                         count_dyad = test_df$count_dyad,
                         age_1 = test_df$age_category_1,
                         age_2 = test_df$age_category_2,
                         age_min = test_df$age_min_mu,
                         age_diff = test_df$age_diff_mu,
                         edge_mean = NA)
for(i in 1:nrow(mean_edges)) {
  x <- edges[edges$dyad == mean_edges$dyad[i],]
  mean_edges$edge_mean[i] <- mean(x$edge_draw)
}

### plot raw data
theme_set(theme_classic())
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean),
             shape = 19, colour = rgb(0,0,1,0.1)
  )+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean,
                                    colour = age_diff),
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

edges <- left_join(edges, mean_edges, by = 'dyad')
edges$age_cat_min <- ifelse(edges$age_min < 10, '<10',
                     ifelse(edges$age_min < 15, '10-15',
                     ifelse(edges$age_min < 20, '15-20',
                     ifelse(edges$age_min < 25, '20-25',
                     ifelse(edges$age_min < 40, '25-40','40+')))))
edges$age_cat_min <- factor(edges$age_cat_min, levels = c('<10','10-15','15-20','20-25','25-40','40+'))
mean_edges$age_cat_min <- ifelse(mean_edges$age_min < 10, '<10',
                          ifelse(mean_edges$age_min < 15, '10-15',
                          ifelse(mean_edges$age_min < 20, '15-20',
                          ifelse(mean_edges$age_min < 25, '20-25',
                          ifelse(mean_edges$age_min < 40, '25-40','40+')))))
mean_edges$age_cat_min <- factor(mean_edges$age_cat_min, levels = c('<10','10-15','15-20','20-25','25-40','40+'))
ggplot()+
  geom_violin(data = edges, aes(x = age_cat_min, y = edge_draw))+
  geom_jitter(data = mean_edges, aes(x = age_cat_min, y = edge_mean, size = count_dyad),
              shape = 1, width = 0.3, colour = rgb(0,0,1,0.5))+
  scale_x_discrete('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  coord_cartesian(ylim = c(0,0.2))

ggplot()+
  geom_point(data = edges, aes(x = count_dyad, y = edge_draw), colour = rgb(0,0,1,0.01))+
  geom_point(data = mean_edges, aes(x = count_dyad, y = edge_mean), colour = 'red')+
  scale_x_continuous('number of sightings of dyad')+
  scale_y_continuous('edge weight')

## heat map: x = min age, y = diff age, colour = edge weight
ggplot()+
  geom_tile(data = mean_edges,
            mapping = aes(x = round(age_min,0), y = round(age_diff,0), fill = edge_mean))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age difference')

# add time marker
print(paste0('age vs edge predictions completed at ', Sys.time()))

#### plot predictions ####
## posterior predictive check -- CURRENTLY HAVE NOTHING HERE

### end pdf
dev.off()

# save image
save.image('motnp_dyadicregression.RData')

### clear workspace
rm(counts_df, n_dyads, counts_ls, n_eles, fit_edges_motnp, edges, plot_dyads, plot_edges, nodes, edge_weights_matrix, edge_samples)

# add time marker
print(paste0('completed at ', Sys.time()))
