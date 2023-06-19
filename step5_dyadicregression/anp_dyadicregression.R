#### information ####
# script takes data input from edge weight estimation for ANP population (input = 1000 draws for each of 4 chains per dyad posterior distribution)
# runs through dyadic regression as specified by Jordan Hart in BISoN examples (https://github.com/JHart96/bison_examples/blob/main/examples/dyadic_regression_stan.md)

#### set up ####
#library(tidyverse) ; library(car) ; library(cmdstanr) ; library(bisonR) ; library(brms)
library(tidyverse, lib.loc = '../packages/') # library(tidyverse)
library(car, lib.loc = '../packages/')       # library(car)
library(cmdstanr, lib.loc = '../packages/')  # library(cmdstanr)
library(bisonR, lib.loc = '../packages/')    # library(bisonR)
library(brms, lib.loc = '../packages/')      # library(brms)
#library(rstan, lib.loc = '../packages/')     # library(rstan)
#library(Rcpp, lib.loc = '../packages/')      # library(Rcpp)

#set_cmdstan_path('R:/rsrch/df525/phd/hkm513/packages/.cmdstan/cmdstan-2.31.0')

load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
rm(edgelist, counts_df, edge_binary, nodes) ; gc()

pdf('../outputs/anp1_dyadicregression_withmm_plots.pdf')

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

#### fit multivariate Gaussian distribution to output of edge weight model ####
# edge weights from beta() model stored in edge_samples, where column is dyad and row is edge draw (1:4000)
# note in bison papers, Jordan calculated edges under a logit-normal model
trfm_edges <- (edge_samples)^(1/3)                   # transform edges to be more normal: log transform takes too left skewed, plogis doesn't change shape at all
plot(density(edge_samples[, 1]), lwd = 2, las = 1,
     main = "edge transformation: black = edges, blue = cube root edges", xlab = "edge weight", xlim = c(0,1))
lines(density(trfm_edges[, 1]), col = rgb(0,0,1,0.5), lwd = 2)

# calculate mean of logit edge samples
edge_mu <- apply(trfm_edges, 2, mean)

# calculate covariance of logit edge samples
edge_cov <- cov(trfm_edges)

# compute multivariate Gaussian approximation
mv_edge_samples <- MASS::mvrnorm(1e5, edge_mu, edge_cov)

# plot edge weight vs normal approximation
par(mai = c(1,1,1.2,0.4))
plot(density(edge_samples[, 1]), lwd = 2, las = 1, xlim = c(0,1), ylim = c(0,5), xlab = "edge weight",
     main = "normal approximation of edges:\nblack = raw edge weights\nred = cube root edge weights\nblue = mvnorm approximation")
lines(density(trfm_edges[, 1]), col = 'red', lwd = 2)
lines(density(mv_edge_samples[, 1]), col = 'blue', lwd = 2)

# plot covariance of edges
plot(mv_edge_samples[, 1], mv_edge_samples[, 2], col = rgb(0,0,1,0.05), las = 1,
     main = "Covariance between edges 1 & 2", xlab = "Edge 1 samples", ylab = "Edge 2 samples")

# save image so far
save.image('anpshort1_dyadicregression_withmm.RData')

#### fit dyadic regression ####
### identify older and younger of dyad
cdf_1$age_min <- NA ; cdf_1$age_max <- NA
for(i in 1:nrow(cdf_1)){
  x <- c(cdf_1$age_start_1[i],cdf_1$age_start_2[i])
  cdf_1$age_min[i] <- min(x)
  cdf_1$age_max[i] <- max(x)
}

## create data list
dyad_data <- list(
  num_dyads = n_dyads,                                    # number of dyads
  num_nodes = length(unique(c(cdf_1$id_1, cdf_1$id_2))),  # number of nodes
  edge_mu = edge_mu,                                      # sample means of logit edge weights
  edge_cov = edge_cov,                                    # sample covariance of logit edge weights
  age_min = cdf_1$age_min,                                # age of younger dyad member
  age_max = cdf_1$age_max,                                # age of  older  dyad member
  node_1 = cdf_1$node_1,                                  # node IDs for multimembership effects
  node_2 = cdf_1$node_2                                   # node IDs for multimembership effects
)

## load dyadic regression model
dyadic_regression <- cmdstan_model('models/dyadic_regression.stan')

## fit dyadic regression
fit_dyadreg_anp1 <- dyadic_regression$sample(
  data = dyad_data,
  chains = n_chains,
  parallel_chains = n_chains)

## save image
save.image('anpshort1_dyadicregression_conditionaledge.RData')

#### check outputs ####
#load('anpshort1_dyadicregression_conditionaledge')
fit_dyadreg_anp1$summary()

## extract draws
draws <- fit_dyadreg_anp1$draws()

## extract dyadic regression slopes
b_min <- draws[,,'b_min']
b_max <- draws[,,'b_max']
#b_int <- draws[,,'b_int']
sigma <- draws[,,'sigma']

## extract overall age effects
draws <- draws[,,5:length(draws[1,1,])]   # ONCE ADD INTERACTION, CHANGE TO: draws <- draws[,,6:length(draws[1,1,])]
age_effects <- draws[,,1:n_dyads]

## extract multimembership samples -- CHECK THIS IS IN THE RIGHT ORDER
#draws <- draws[,,(n_dyads+1):length(draws[1,1,])]
#mm_matrix <- draws[,,1:n_dyads]
#sigma_mm <- draws[,,(n_dyads+1)]
#length(draws[1,1,]) == length(mm_matrix[1,1,]) + length(sigma_mm[1,1,])
rm(draws) ; gc()

## save parameter values
extract_slopes <- function(draws, n_samples = 1000, n_chains = 4){
  draws <- as.data.frame(draws) %>% 
    pivot_longer(everything(), names_to = 'chain', values_to = 'slope_draw') %>% 
    separate(chain, sep = c(1,2), remove = T, into = c('chain','.', 'parameter')) %>% 
    select(-.)
  draws$chain <- as.numeric(draws$chain)
  draws$position <- rep(1:n_samples, each = n_chains)
  draws$draw_id <- draws$position + (draws$chain-1)*n_samples
  return(draws)
}
b_min <- extract_slopes(b_min)
b_max <- extract_slopes(b_max)
sigma <- extract_slopes(sigma)
#b_int <- extract_slopes(b_int)
parameters <- rbind(b_min, b_max, #b_int,
                    sigma)

## save data 
saveRDS(parameters, '../data_processed/anp1_dyadicregression_slopeparameters.RDS')

## traceplots
ggplot(data = parameters)+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme_classic()+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap( . ~ parameter )

#### plot edges against age values ####
#edges <- readRDS('../data_processed/anpshort1_edgedistributions_conditionalprior.RDS')

### create mean data frame to plot average values over full distribution
mean_edges <- data.frame(dyad = cdf_1$dyad_id,
                         node_1 = cdf_1$node_1, node_2 = cdf_1$node_2,
                         together = cdf_1$event_count, count_dyad = cdf_1$period_count_dyad,
                         age_1 = cdf_1$age_start_1, age_2 = cdf_1$age_start_2,
                         age_min = cdf_1$age_min, age_max = cdf_1$age_max, age_diff = cdf_1$age_diff,
                         edge_mean = NA)
for(i in 1:nrow(mean_edges)) {
  x <- edges[edges$dyad == mean_edges$dyad[i],]
  mean_edges$edge_mean[i] <- mean(x$edge_draw)
}

### add categories for box plots
mean_edges$age_cat_min <- ifelse(mean_edges$age_min < 15, '10-15',
                                 ifelse(mean_edges$age_min < 20, '15-20',
                                        ifelse(mean_edges$age_min < 25, '20-25',
                                               ifelse(mean_edges$age_min < 40, '25-40', '40+'))))
mean_edges$age_cat_max <- ifelse(mean_edges$age_max < 15, '10-15',
                                 ifelse(mean_edges$age_max < 20, '15-20',
                                        ifelse(mean_edges$age_max < 25, '20-25',
                                               ifelse(mean_edges$age_max < 40, '25-40', '40+'))))

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
                                    colour = age_max), # this makes it look more interesting than it is -- there is a definite pattern at first glance, but actually appears to be nothing more than that the older the youngest male is, the older the oldest must be to be older than the youngest, whereas when the youngest is very young the oldest doesn't have to be very old to still be the oldest
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean, colour = age_diff), # much less exciting now
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

edges <- left_join(edges, mean_edges, by = 'dyad')
ggplot()+
  geom_point(data = edges, aes(x = age_min, y = edge_draw),
             shape = 19, colour = rgb(0,0,1,0.01))+
  geom_point(data = mean_edges, aes(x = age_min, y = edge_mean),
             shape = 19, colour = 'red')+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')
ggplot()+
  geom_violin(data = edges, aes(x = age_cat_min, y = edge_draw))+
  geom_jitter(data = mean_edges, aes(x = age_cat_min, y = edge_mean, size = count_dyad),
              shape = 1, width = 0.3, colour = rgb(0,0,1,0.1))+
  scale_x_discrete('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')

ggplot()+
  geom_point(data = edges, aes(x = count_dyad, y = edge_draw), colour = rgb(0,0,1,0.01))+
  geom_point(data = mean_edges, aes(x = count_dyad, y = edge_mean), colour = 'red')+
  scale_x_continuous('number of sightings of dyad')+
  scale_y_continuous('edge weight')

ggplot()+
  geom_point(data = mean_edges, aes(x = age_min, y = age_max, colour = edge_mean))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age of older dyad member')

## heat map: x = min age, y = max age, colour = edge weight
ggplot()+
  geom_tile(data = mean_edges, mapping = aes(x = age_min, y = age_max, fill = edge_mean))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('age of older dyad member')

#### plot predictions ####
## posterior predictive check
plot(density(edge_samples[1, ]), main = "Posterior predictive density of responses (edge weights)",
     ylim = c(0, 12), col = rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:1000, 1)
  lines(density(edge_samples[j, ]), col = rgb(0, 0, 0, 0.25))
  mu_plot <- b_min$slope_draw[j]*dyad_data$age_min + b_max$slope_draw[j]*dyad_data$age_max
  sigma_plot <- dyad_data$edge_cov + diag(rep(sigma$slope_draw[j], n_dyads))
  mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)
  lines(density(mv_norm^3), col = rgb(1, 0, 0, 0.25))
}

# compute posterior means
age_min <- seq(5,46,1)
age_max <- seq(5,46,1)
mu <- array(NA, dim = c(n_samples, length(age_min), length(age_max)), dimnames = list(1:n_samples, age_min, age_max))
for(sample in 1:n_samples){
  for(min in 1:length(age_min)){
    for(max in 1:length(age_max)){
      if( min <= max ){
        mu[sample,min,max] <- b_min$slope_draw[sample]*age_min[min] + b_max$slope_draw[sample]*age_max[max]# + b_int$slope_draw[sample]*age_min[min]*age_max[max]
      }
    }
  }
}
#mu <- plogis(mu)     # ONLY NECESSARY IF PRODUCING VALUES ON LOGIT NORMAL SCALE RATHER THAN BETA??
mu <- mu^3           # reverse data transformation

mu_max_age <- function(edge_prediction_array, age_threshold){
  minimum_ages <- as.numeric(colnames(edge_prediction_array[1,,]))
  maximum_ages <- as.numeric(names(edge_prediction_array[1,1,]))
  predictions <- edge_prediction_array[,which(minimum_ages <= age_threshold),which(maximum_ages == age_threshold)]
  if(is.null(dim(predictions)) == TRUE) { mean_edge <- mean(predictions) } else { mean_edge <- apply(predictions, 2, mean) }
  return(mean_edge)
}

pi_max_age <- function(edge_prediction_array, age_threshold){
  minimum_ages <- as.numeric(colnames(edge_prediction_array[1,,]))
  maximum_ages <- as.numeric(names(edge_prediction_array[1,1,]))
  predictions <- edge_prediction_array[,which(minimum_ages <= age_threshold),which(maximum_ages == age_threshold)]
  if(is.null(dim(predictions)) == TRUE) { 
    mean_edge <- rethinking::HPDI(predictions) 
  } else { 
    mean_edge <- apply(predictions, 2, rethinking::HPDI) 
  }
  return(mean_edge)
}

ages <- c(20,25,30,35,41,46)

mu_mean_max20 <- mu_max_age(mu, 20)
mu_mean_max25 <- mu_max_age(mu, 25)
mu_mean_max30 <- mu_max_age(mu, 30)
mu_mean_max35 <- mu_max_age(mu, 35)
mu_mean_max41 <- mu_max_age(mu, 41)
mu_mean_max46 <- mu_max_age(mu, 46)

pi_mean_max20 <- pi_max_age(mu, 20)
pi_mean_max25 <- pi_max_age(mu, 25)
pi_mean_max30 <- pi_max_age(mu, 30)
pi_mean_max35 <- pi_max_age(mu, 35)
pi_mean_max41 <- pi_max_age(mu, 41)
pi_mean_max46 <- pi_max_age(mu, 46)

## simulate from posterior -- NOT YET WORKED THIS BIT OUT
sim_edges <- matrix(nrow = n_samples, ncol = length(age_min))
for(i in 1:nrow(sim_edges)){
  for(j in 1:ncol(sim_edges)){
    mean_sim <- b_min$slope_draw[i]*age_min[j] + b_max$slope_draw[i]*age_max[j]
    sigma_sim <- sigma$slope_draw[i]
    x <- MASS::mvrnorm(n = 1, mu = (mean_sim^3), Sigma = sigma_sim)
    sim_edges[i,j] <- x#^3
  }
}
hist(sim_edges)  # almost entirely equal
sim_pi <- apply(sim_edges, 2, rethinking::HPDI, prob = 0.95)

# plot
par(mfrow = c(3,2))
mu_list <- list(mu_mean_max20, mu_mean_max25, mu_mean_max30, mu_mean_max35, mu_mean_max41, mu_mean_max46)
pi_list <- list(pi_mean_max20, pi_mean_max25, pi_mean_max30, pi_mean_max35, pi_mean_max41, pi_mean_max46)
par(mai = c(0.6,0.6,0.6,0.4))
for(i in 1:6){
  plot(edge_mean ~ age_min, data = mean_edges[mean_edges$age_max == ages[i],], col = rgb(0,0,1,0.5), # raw sightings
       pch = 19, las = 1, xlim = c(5,45), ylim = c(0,0.5),
       xlab = 'age of younger dyad member', ylab = 'edge weight',
       main = paste0('older = ',ages[i],' years old'))
  rethinking::shade(pi_list[[i]], age_min[which(age_min <= ages[i])], lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean line
  lines(age_min[which(age_min <= ages[i])], mu_list[[i]], lwd = 2, col = 'purple')                          # add mean shading
  rethinking::shade(sim_pi, age_min)                                                                   # add predictions
}
par(mfrow = c(1,1))

## heat map: x = min age, y = max age, colour = edge weight
mean_predict <- matrix(NA, nrow = length(age_max), ncol = length(age_min))
for(i in 1:nrow(mean_predict)){
  for(j in 1:ncol(mean_predict)){
    if(i <= j){
      x <- mu_max_age(mu, age_max[j])
      mean_predict[,j] <- c(x ,rep(NA, length(mean_predict[,j]) - length(x)))
    }
  }
}
colnames(mean_predict) <- paste0('max.',age_max)
rownames(mean_predict) <- paste0('min.',age_min)
mean_predict_long <- as.data.frame(mean_predict) %>% 
  pivot_longer(everything(), names_to = 'age_max', values_to = 'mean_edge_weight') %>% 
  separate(age_max, into = c('max.min','age_max'), sep = 4) %>% 
  select(-max.min) %>% 
  mutate(age_min = rep(rownames(mean_predict), each = ncol(mean_predict))) %>% 
  separate(age_min, into = c('max.min','age_min'), sep = 4) %>% 
  select(-max.min) %>% 
  mutate(age_diff = as.numeric(age_max) - as.numeric(age_min))
#filter(!is.na(mean_edge_weight)) %>% 
mean_predict_long$age_max <- as.numeric(mean_predict_long$age_max)
mean_predict_long$age_min <- as.numeric(mean_predict_long$age_min)
ggplot()+
  geom_tile(data = mean_predict_long, mapping = aes(x = age_min, y = age_max, fill = mean_edge_weight))
ggplot()+
  geom_line(data = mean_predict_long, mapping = aes(x = age_diff, y = mean_edge_weight, colour = age_min, group = age_min))+
  scale_x_continuous(limits = c(0,40))

### end pdf
dev.off()

### clear workspace
rm(cdf_1, n_dyads, counts_ls, n_eles, fit_edges_anp1, edges, plot_dyads, plot_edges, nodes, edge_weights_matrix, edge_samples)

