#### Information ####
# 1) to estimate the ages of individuals in other populations, first came up with a survival curve for the Amboseli elephants which have good birth and death data.
# 2) selected the Gompertz bathtub distribution
# 3) use distribution to convert MOTNP ages in categories to probability distributions of true age (first part of this script)
# 4) extract node centrality distribution from draws of edge weight
# 5) use MOTNP age distributions to predict network centrality (third part of this script)

#### set up ####
options(future.globals.maxSize = 10000*(1024^2))   # running model with full age distributions = 5.01GB of globals, which exceeds default maximum allowed size. Set to allow up to 10 GB to allow model to run

# library(tidyverse) ; library(cmdstanr) ; library(brms) ; library(Rcpp) ; library(ggdist) ; library(posterior) ; library(bayesplot) ; library(igraph) ; library(LaplacesDemon) ; library(bisonR) ; library(janitor)
library(dplyr, lib.loc = '../packages/')       # library(tidyverse)
library(rstan, lib.loc = '../packages/')           # library(cmdstanr)
#library(janitor, lib.loc = '../packages/')         # library(janitor)

## load work space from calculating edge weights
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')

#### dyadic regression ####
## define PDF output
pdf('../outputs/anpshort1_dyadicregression_conditionalprior.pdf')

## prior predictive check ####
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

## fit multivariate Gaussian distribution to output of edge weight model ####
# get the weights on the logit scale (they are not currently because we used a beta prior and identity link here rather than logistic link)
logit_weights <- apply(edge_samples, 2, qlogis)

# fit a multivariate normal dist to the edges
logit_edge_draws_mu <- apply(logit_weights, 2, mean)
logit_edge_draws_cov <- cov(logit_weights)

## plot to see how well the approximation is working ####
# Randomly selecting samples to examine
num_check <- 20
selected_samples <- sample(1:n_dyads, num_check, replace = FALSE)

# Setting grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

# plot
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
par(mfrow=c(1,1))

# save image so far
save.image('anpshort1_dyadicregression.RData')

# add time marker
print(paste0('multivariate Gaussian approximation fitted at ', Sys.time()))

## fit dyadic regression ####
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
  #age_max = cdf_1$age_max,                  # age of  older  dyad member
  age_diff = cdf_1$age_diff,                # age difference between dyad members
  node_1 = cdf_1$node_1_period,             # node IDs for multimembership effects
  node_2 = cdf_1$node_2_period              # node IDs for multimembership effects
)

## load dyadic regression model
dyadic_regression <- stan_model('models/dyadic_regression.stan')

# add time marker
print(paste0('start model run at ', Sys.time()))

## fit dyadic regression
fit_dyadreg_anp1 <- sampling(
  object = dyadic_regression,
  data = dyad_data,
  chains = n_chains,
  cores = n_chains)

# save image so far
save.image('anpshort1_dyadicregression.RData')

# add time marker
print(paste0('finish model run at ', Sys.time()))

## check outputs ####
#load('anpshort1_dyadicregression.RData')
rm(cdf_1, counts_df, dyad_data, edge_binary, edgelist, edges, i, x, make_edgelist, plot_network_threshold_anp) ; gc()

# obtain summary
fit_dyadreg_anp1$summary()

## extract draws
draws <- fit_dyadreg_anp1$draws(format = 'df')

## extract dyadic regression slopes
b_diff <- draws[,,'beta_age_diff']
b_min <- draws[,,'beta_age_min']
sigma <- draws[,,'sigma']

#b_max <- draws[,,'b_max']
#b_int <- draws[,,'b_int']

## extract overall age effects
dimnames(draws)$variable
draws <- draws[,,4:length(draws[1,1,])]   # WILL CHANGE IF YOU ADD ANY MORE REGRESSION SLOPES

## extract multi-membership samples -- CHECK THIS IS IN THE RIGHT ORDER
mm_matrix <- draws[,,1:n_nodes]
sigma_mm <- draws[,,n_nodes+2]
predictor <- draws[,,(n_nodes+3):length(draws[1,1,])]
length(draws[1,1,]) == length(mm_matrix[1,1,]) + length(sigma_mm[1,1,]) + length(sigma[1,1,]) + length(predictor[1,1,])
rm(draws) ; gc()

# add time marker
print(paste0('parameters extracted at ', Sys.time()))

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
b_diff <- extract_slopes(b_diff)
b_min <- extract_slopes(b_min)
#b_max <- extract_slopes(b_max)
sigma <- extract_slopes(sigma)
#b_int <- extract_slopes(b_int)
parameters <- rbind(b_min, b_diff, #b_max, #b_int,
                    sigma)

## save data 
saveRDS(parameters, '../data_processed/anp1_dyadicregression_slopeparameters.RDS')

# add time marker
print(paste0('parameters saved to file at ', Sys.time()))

pdf('../outputs/anp1_dyadicregression_plots2.pdf')

## traceplots
ggplot(data = parameters)+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme_classic()+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap( . ~ parameter , scales = 'free_y')

# add time marker
print(paste0('traceplots run at ', Sys.time()))

## plot edges against age values ####
#edges <- readRDS('../data_processed/anpshort1_edgedistributions_conditionalprior.RDS')

### create mean data frame to plot average values over full distribution
mean_edges <- data.frame(dyad = cdf_1$dyad_id,
                         node_1 = cdf_1$node_1_period, node_2 = cdf_1$node_2_period,
                         together = cdf_1$event_count, count_dyad = cdf_1$period_count_dyad,
                         age_1 = cdf_1$age_start_1, age_2 = cdf_1$age_start_2,
                         age_min = cdf_1$age_min, age_max = cdf_1$age_max, age_diff = cdf_1$age_diff,
                         edge_mean = NA)
for(i in 1:nrow(mean_edges)) {
  x <- edges[edges$dyad_id == mean_edges$dyad[i],]
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

# add time marker
print(paste0('age vs edge predictions completed at ', Sys.time()))

## plot predictions ####
## posterior predictive check
plot(density(edge_samples[1, ]), main = "Posterior predictive density of responses (edge weights)",
     ylim = c(0, 12), col = rgb(0, 0, 0, 0.25))
for (i in 1:100) {
  j <- sample(1:1000, 1)
  mu_plot <- b_min$slope_draw[j]*dyad_data$age_min + b_diff$slope_draw[j]*dyad_data$age_diff
  sigma_plot <- dyad_data$logit_edge_cov + diag(rep(sigma$slope_draw[j], n_dyads))
  mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)
  
  lines(density(edge_samples[j, ]), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(plogis(mv_norm)), col = rgb(1, 0, 0, 0.25))   # red lines for predictions
}










load('TESTING_DYADIC_REGRESSION_17THAUGUST2023.RData') # delete this once you've loaded it back in tomorrow!
rm(counts_df, dyad_data, edge_binary, edgelist, fit_edges_anp1, parameters, sigma_plot, summary, x, all_node_IDs, i, j, mm_matrix, mu_plot, mv_norm, n_windows, periods, predictor, sigma_mm, extract_slopes, make_edgelist, plot_network_threshold_anp) ; gc










# compute posterior means
ages <- min(nodes$age):max(nodes$age)
age_differences <- 0:30
mu <- array(NA, dim = c(n_samples, length(ages), length(age_differences)),
            dimnames = list(1:n_samples, ages, age_differences))
for(draw in 1:n_samples){
  for(minimum in 1:length(ages)){
    for(difference in 1:length(age_differences)){
      #if( min <= diff ){
      mu[draw,minimum,difference] <- b_min$slope_draw[draw]*ages[minimum] + b_diff$slope_draw[draw]*age_differences[difference]# + b_int$slope_draw[sample]*age_min[min]*age_max[max]
      #}
    }
  }
}
mu <- plogis(mu)     # IS THIS ONLY NECESSARY IF PRODUCING VALUES ON LOGIT NORMAL SCALE RATHER THAN BETA?? DON'T THINK SO BEACUSE WORKING FROM THE APPROXIMATION NOT THE ACTUAL EDGE DISTRIBUTIONS

mu[10,,] # 10th draw for all pairs of min and diff
mu[,16,] # 16 years minimum age for all age differences (max age = 16-46)
mu[,,25] # 25 years age difference (max age = 30-45)

#mu_max_age <- function(edge_prediction_array, age_threshold){
#  minimum_ages <- as.numeric(colnames(edge_prediction_array[1,,]))
#  maximum_ages <- as.numeric(names(edge_prediction_array[1,1,]))
#  predictions <- edge_prediction_array[,which(minimum_ages <= age_threshold),which(maximum_ages == age_threshold)]
#  if(is.null(dim(predictions)) == TRUE) { mean_edge <- mean(predictions) } else { mean_edge <- apply(predictions, 2, mean) }
#  return(mean_edge)
#}

#pi_max_age <- function(edge_prediction_array, age_threshold){
#  minimum_ages <- as.numeric(colnames(edge_prediction_array[1,,]))
#  maximum_ages <- as.numeric(names(edge_prediction_array[1,1,]))
#  predictions <- edge_prediction_array[,which(minimum_ages <= age_threshold),which(maximum_ages == age_threshold)]
#  if(is.null(dim(predictions)) == TRUE) { 
#    mean_edge <- rethinking::HPDI(predictions) 
#  } else { 
#    mean_edge <- apply(predictions, 2, rethinking::HPDI) 
#  }
#  return(mean_edge)
#}

#mu_mean_max20 <- mu_max_age(mu, 20)
#mu_mean_max25 <- mu_max_age(mu, 25)
#mu_mean_max30 <- mu_max_age(mu, 30)
#mu_mean_max35 <- mu_max_age(mu, 35)
#mu_mean_max41 <- mu_max_age(mu, 41)
#mu_mean_max46 <- mu_max_age(mu, 46)

#pi_mean_max20 <- pi_max_age(mu, 20)
#pi_mean_max25 <- pi_max_age(mu, 25)
#pi_mean_max30 <- pi_max_age(mu, 30)
#pi_mean_max35 <- pi_max_age(mu, 35)
#pi_mean_max41 <- pi_max_age(mu, 41)
#pi_mean_max46 <- pi_max_age(mu, 46)

mu_diff_age <- function(edge_prediction_array, min_threshold, diff_threshold){
  minimum_ages <- as.numeric(rownames(edge_prediction_array[1,,]))
  difference_ages <- as.numeric(names(edge_prediction_array[1,1,]))
  predictions <- edge_prediction_array[,which(minimum_ages <= min_threshold),
                                       which(difference_ages == diff_threshold)]
  if(is.null(dim(predictions)) == TRUE) { 
    mean_edge <- mean(predictions) 
  } else { 
    mean_edge <- apply(predictions, 2, mean) 
  }
  return(mean_edge)
}

pi_diff_age <- function(edge_prediction_array, min_threshold, diff_threshold){
  minimum_ages <- as.numeric(rownames(edge_prediction_array[1,,]))
  difference_ages <- as.numeric(names(edge_prediction_array[1,1,]))
  predictions <- edge_prediction_array[,which(minimum_ages <= min_threshold),
                                       which(difference_ages == diff_threshold)]
  if(is.null(dim(predictions)) == TRUE) { 
    hpdi_edge <- rethinking::HPDI(predictions) 
  } else { 
    hpdi_edge <- apply(predictions, 2, rethinking::HPDI) 
  }
  return(hpdi_edge)
}

mu_means_list <- data.frame(min_age = rep(ages, each = length(age_differences)),
                            diff_age = rep(age_differences, length(ages)),
                            mean = NA, hpdi = NA)
for(i in 1:nrow(mu_means_list)){
  mu_means_list$mean[i] <- list(mu_diff_age(edge_prediction_array = mu,
                                            min_threshold = mu_means_list$min_age[i],
                                            diff_threshold = mu_means_list$diff_age[i]))
  mu_means_list$hpdi[i] <- list(pi_diff_age(edge_prediction_array = mu,
                                            min_threshold = mu_means_list$min_age[i],
                                            diff_threshold = mu_means_list$diff_age[i]))
}
hpdi <- as.data.frame(mu_means_list$hpdi[i])
hpdi_lower <- hpdi[1,] %>% 
  pivot_longer(everything(), names_to = 'x_min', values_to = 'lower_bound') %>% 
  separate(x_min, into = c('x','min_age'), remove = T, sep = 1) %>% 
  mutate(min_age = as.numeric(min_age)) %>% 
  select(-x)
hpdi_upper <- hpdi[2,] %>% 
  pivot_longer(everything(), names_to = 'x_min', values_to = 'upper_bound') %>% 
  separate(x_min, into = c('x','min_age'), remove = T, sep = 1) %>% 
  mutate(min_age = as.numeric(min_age)) %>% 
  select(-x)
hpdi <- left_join(hpdi_lower, hpdi_upper, by = 'min_age') %>% 
  mutate(diff_age = mu_means_list$diff_age[i]) %>% 
  relocate(diff_age, .after = min_age)

for(i in 1:nrow(mu_means_list)){
  hpdi_loop <- as.data.frame(mu_means_list$hpdi[i])
  if(ncol(hpdi_loop) > 1){
    hpdi_loop_lower <- hpdi_loop[1,] %>% 
      pivot_longer(everything(), names_to = 'x_min', values_to = 'lower_bound') %>% 
      separate(x_min, into = c('x','min_age'), remove = T, sep = 1) %>% 
      mutate(min_age = as.numeric(min_age)) %>% 
      select(-x)
    hpdi_loop_upper <- hpdi_loop[2,] %>% 
      pivot_longer(everything(), names_to = 'x_min', values_to = 'upper_bound') %>% 
      separate(x_min, into = c('x','min_age'), remove = T, sep = 1) %>% 
      mutate(min_age = as.numeric(min_age)) %>% 
      select(-x)
    hpdi_loop <- left_join(hpdi_loop_lower, hpdi_loop_upper, by = 'min_age') %>% 
      mutate(diff_age = mu_means_list$diff_age[i])
  } else {
    hpdi_loop <- data.frame(min_age = mu_means_list$min_age[i],
                            diff_age = mu_means_list$diff_age[i],
                            lower_bound = hpdi_loop[1,1],
                            upper_bound = hpdi_loop[2,1])
  }
  hpdi <- rbind(hpdi, hpdi_loop)
}

mu_means <- hpdi %>% 
  left_join(mu_means_list[,c('min_age','diff_age','mean')],
            by = c('min_age','diff_age')) %>% 
  unnest(cols = mean) %>% 
  distinct()

rm(hpdi, hpdi_loop, hpdi_loop_lower, hpdi_loop_upper, hpdi_lower, hpdi_upper) ; gc()

for(difference in 1:length(age_differences)){
  
  ## simulate from posterior
  sim_edges <- matrix(nrow = n_samples, ncol = length(ages))
  for(draw in 1:nrow(sim_edges)){
    for(minimum in 1:ncol(sim_edges)){
      mean_sim <- b_min$slope_draw[draw]*ages[minimum] + b_diff$slope_draw[draw]*age_differences[difference]
      sigma_sim <- sigma$slope_draw[draw]
      sim_edges[draw,minimum] <- MASS::mvrnorm(n = 1, mu = mean_sim, Sigma = sigma_sim)
      #sim_edges[i,j] <- LaplacesDemon::rmvnc(n = 1, mu = mean_sim, U = sigma_sim)
    }
  }
  
  hist(sim_edges)  # almost entirely equal
  sim_mu <- apply(sim_edges, 2, mean)
  sim_pi <- apply(sim_edges, 2, rethinking::HPDI, prob = 0.95)
  
  # plot
  plot(data = mean_edges[mean_edges$age_diff == age_differences[difference],],
       edge_mean ~ age_min, col = rgb(0,0,1,0.5), # raw sightings
       pch = 19, las = 1, xlim = c(5,45), ylim = c(0,1),
       xlab = 'minimum age', ylab = 'edge weight',
       main = paste0('difference = ',age_differences[difference],' years'))
  #rethinking::shade(pi_list[[i]],
  #                  age_min[which(age_min <= ages[i])],
  #                  lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean shading
  #polygon(y = plogis(sim_pi),
  #        x = c(ages,ages),
  #        lwd = 1, col = rgb(0.5,0,1,0.2))      # add mean shading
  polygon(y = c(plogis(sim_pi[1,]),plogis(sim_pi[2,])),
          x = c(rev(ages),ages),
          lwd = 1, col = rgb(0.5,0,1,0.2))       # add mean shading
  lines(x = ages,
        y = plogis(sim_mu), lwd = 2, col = 'purple')   # add mean line
  #rethinking::shade(sim_pi, age_min)             # add predictions
}




## simulate from posterior -- IS THIS RIGHT??
sim_edges <- matrix(nrow = n_samples, ncol = length(age_min))
for(i in 1:nrow(sim_edges)){
  for(j in 1:ncol(sim_edges)){
    mean_sim <- b_min$slope_draw[i]*age_min[j] + b_diff$slope_draw[i]*age_diff[j]
    sigma_sim <- sigma$slope_draw[i]
    sim_edges[i,j] <- MASS::mvrnorm(n = 1, mu = mean_sim, Sigma = sigma_sim)
    #sim_edges[i,j] <- LaplacesDemon::rmvnc(n = 1, mu = mean_sim, U = sigma_sim)
  }
}
hist(sim_edges)  # almost entirely equal
sim_pi <- apply(sim_edges, 2, rethinking::HPDI, prob = 0.95)

# plot
par(mfrow = c(length(unique(mu_means$min_age)),length(unique(mu_means$diff_age))))
par(mai = c(1,1,1,0.6))
for(i in 1:length(ages)){
  plot(data = mean_edges[mean_edges$age_min == ages[i],],
       edge_mean ~ age_diff, col = rgb(0,0,1,0.5), # raw sightings
       pch = 19, las = 1, xlim = c(5,45), ylim = c(0,0.5),
       xlab = 'age difference', ylab = 'edge weight',
       main = paste0('younger = ',ages[i],' years old'))
  #rethinking::shade(pi_list[[i]],
  #                  age_min[which(age_min <= ages[i])],
  #                  lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean shading
  polygon(y = sim_pi[,i],
          x = 0:30,
          lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean shading
  #polygon(x = age_diff, y = 
  #                  age_min[which(age_min <= ages[i])],
  #                  lwd = 2, col = rgb(0.5,0,1,0.2))      # add mean shading
  lines(age_min[which(age_min <= age_min[i])],
        mu_list[[i]], lwd = 2, col = 'purple')            # add mean line
  rethinking::shade(sim_pi, age_min)                      # add predictions
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

# save image
save.image('anpshort1_dyadicregression.RData')

### clear workspace
rm(cdf_1, n_dyads, counts_ls, n_eles, fit_edges_anp1, edges, plot_dyads, plot_edges, nodes, edge_weights_matrix, edge_samples)

# add time marker
print(paste0('completed at ', Sys.time()))






























# #### nodal regression -- rstan ####
# ## define PDF output
# pdf('../outputs/anpshort1_nodalregression_conditionalprior.pdf')
# 
# ## prior predictive check ####
# ## simulate
# age <- 1:60
# beta_mu <- 0
# beta_sigma <- 0.005
# mean_age <- mean(nodes$age)
# plot(NULL, xlim = c(10,60), ylim = c(0,1), las = 1,
#      xlab = 'age', ylab = 'eigenvector centrality')
# for(i in 1:100){
#   intercept <- rbeta(1,2,2)   # is this right?? I've gone for a symmetrical one here that in itself explores most of the parameter space and allows it to see whether some of the lines are steep enough to go from top to bottom, but on the assumption that when combined, they will explore only the space relevant to their starting position
#   beta <- rnorm(1, beta_mu, beta_sigma)
#   lines(x = age, y = intercept + (age - mean_age)*beta, col = rgb(0,0,1,0.5)) # vast majority come out somewhere sensible, and those that don't would if they started at a different value for age 10 so fine in combo with posterior intercept
# }
# 
# ## run model -- rstan ####
# load('anp_nodalregression/anpshort1_nodalregression_conditionaledge.RData')
# 
# ## create data list
# eigen_all <- eigen[,,1]
# for(i in 2:(n_samples*n_chains)){
#   eigen_all <- rbind(eigen_all, eigen[,,i])
# }
# eigen_all <- as.data.frame(eigen_all) %>% 
#   left_join(nodes, by = c('node','age','sightings'))
# eigen_wide <- matrix(NA, nrow = n_samples*n_chains, ncol = n_eles,
#                      dimnames = list(1:(n_samples*n_chains),
#                                      ele_ids))
# for(i in 1:n_eles){
#   node_eigen <- eigen_all %>%
#     filter(id == colnames(eigen_wide)[i])
#   eigen_wide[,i] <- node_eigen$eigenvector
# }
# rm(node_eigen) ; gc()
# mean_eigen <- apply(eigen_wide, 2, mean) %>% 
#   as.data.frame()
# mean_eigen$id <- rownames(mean_eigen)
# colnames(mean_eigen)[1] <- 'mean_eigen'
# nodes <- nodes %>% left_join(mean_eigen,  by = 'id')
# 
# stdv_eigen <- apply(eigen_wide, 2, sd) %>% 
#   as.data.frame()
# stdv_eigen$id <- rownames(stdv_eigen)
# colnames(stdv_eigen)[1] <- 'stdv_eigen'
# nodes <- nodes %>% left_join(stdv_eigen,  by = 'id')
# 
# eigen_list <- list(num_nodes = nrow(nodes),
#                    nodes = nodes$node,
#                    centrality_mu = nodes$mean_eigen,
#                    centrality_sd = nodes$stdv_eigen,
#                    node_age = nodes$age,
#                    node_age2 = (nodes$age)^2)
# 
# # set prior
# nodal_regression <- cmdstan_model('models/eigen_regression.stan')
# 
# ## run model
# fit_anp1_eigen <- nodal_regression$sample(
#   data = eigen_list,
#   chains = n_chains,
#   parallel_chains = n_chains
# )
# 
# ## save output
# save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_cmdstan.RData')
# 
# ## posterior check ####
# # load('anp_nodalregression/anpshort1_nodalregression_conditionaledge_cmdstan.RData')
# # extract draws
# post_eigen <- as.data.frame(as_draws_df(fit_anp1_eigen)) %>% 
#   janitor::clean_names() %>% 
#   select(-lp) %>% 
#   pivot_longer(cols = 1:last_col(3))
# 
# # traceplot linear effect size
# b_age <- post_eigen %>% filter(name == 'beta_age')
# plot(data = b_age[b_age$chain == 1,], value ~ iteration, col = 'purple',
#      type = 'l', xlim = c(0,1000))
# lines(data = b_age[b_age$chain == 2,], value ~ iteration, col = 'red')
# lines(data = b_age[b_age$chain == 3,], value ~ iteration, col = 'blue')
# lines(data = b_age[b_age$chain == 4,], value ~ iteration, col = 'green')
# 
# # traceplot quadratic effect size
# b_age2 <- post_eigen %>% filter(name == 'beta_age2')
# plot(data = b_age2[b_age2$chain == 1,], value ~ iteration, col = 'purple',
#      type = 'l', xlim = c(0,1000))
# lines(data = b_age2[b_age2$chain == 2,], value ~ iteration, col = 'red')
# lines(data = b_age2[b_age2$chain == 3,], value ~ iteration, col = 'blue')
# lines(data = b_age2[b_age2$chain == 4,], value ~ iteration, col = 'green')
# 
# hist(b_age$value)   # natural scale? should it be: hist(plogis(b_age$value))
# hist(b_age2$value)  # hist(plogis(b_age2$value))
# 
# # plot raw data
# ggplot()+
#   geom_point(data = eigen_all, mapping = aes(x = age, y = eigenvector),
#              colour = rgb(253/255, 231/255, 37/255, 0.01))+
#   geom_point(data = nodes,
#              mapping = aes(x = age, y = mean_eigen, size = sightings),
#              colour = rgb(68/255, 1/255, 84/255))+
#   scale_x_continuous('age (years)')+
#   scale_y_continuous('eigenvector centrality')+
#   theme_bw()+
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 16))
# 
# # compare empirical distribution to posterior predictive distribution
# nodes$node_rank <- as.integer(as.factor(nodes$node))
# mean_predict <- post_eigen %>% 
#   filter(name != 'beta_age') %>% filter(name != 'beta_age2') %>% 
#   filter(name != 'sigma') %>% 
#   group_by(name) %>% 
#   mutate(mean_predict = mean(value)) %>% 
#   ungroup() %>% 
#   separate(name, into = c('predictor','node_rank'), remove = F, sep = '_') %>% 
#   select(-predictor) %>% 
#   mutate(node_rank = as.integer(node_rank)) %>% 
#   rename(predicted_eigen = value) %>% 
#   left_join(nodes, by = 'node_rank')
# mean_predict %>% 
#   select(mean_eigen, mean_predict, age, sightings) %>% 
#   distinct() %>% 
#   ggplot()+
#   geom_point(aes(x = mean_eigen, y = mean_predict, colour = age, size = sightings))+
#   scale_x_continuous('mean calculated eigenvector')+
#   scale_y_continuous('mean predicted eigenvector')+
#   theme_classic()
# ggplot(mean_predict)+
#   geom_point(aes(x = mean_eigen, y = predicted_eigen), colour = rgb(0,0,1,0.01))+
#   geom_point(aes(x = mean_eigen, y = mean_predict), colour = 'white')+
#   scale_x_continuous('mean calculated eigenvector')+
#   scale_y_continuous('mean predicted eigenvector')+
#   theme_classic()
# 
# ## save output
# save.image('anp_nodalregression/anpshort1_nodalregression_conditionaledge_cmdstan.RData')
# dev.off()
