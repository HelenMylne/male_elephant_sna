#### set up ####
# script to simulate dyadic regression and check actually working
# library(StanHeaders) ; library(rstan) ; library(tidyverse) ; library(car) ; library(LaplacesDemon) ; library(patchwork)
library(StanHeaders, lib.loc = '../packages/')    # library(StanHeaders)
library(rstan, lib.loc = '../packages/')          # library(rstan)
library(tidyverse, lib.loc = '../packages/')      # library(tidyverse)
library(car, lib.loc = '../packages/')            # library(car)
library(LaplacesDemon, lib.loc = '../packages/')  # library(LaplacesDemon)
library(patchwork, lib.loc = '../packages/')  # library(LaplacesDemon)

theme_set(theme_bw())

# pdf('step5_dyadicregression/simulate_dyadic_anp_diffmean_setup.pdf')

set.seed(12345)

#### load raw data as basis for simulation ####
load('step5_dyadicregression/anpshort_dyadicregression_dataimported.RData')
edge_summary <- edge_summary %>%
  #filter(period < 10) %>%
  rename(window = period) %>%
  mutate(dyad_random = as.integer(as.factor(dyad_random))) # retain randomised order, but remove all of the missing dyads caused by filtering out windows 4:36
#node_data <- node_data %>% filter(window < 10)

# ## cut down number of elephants
# sample_nodes <- c(paste0(sample(x = node_data$id[node_data$window == 1],
#                                 size = ceiling(length(which(node_data$window == 1))/3),
#                                 replace = F),
#                          '_1'),
#                   paste0(sample(x = node_data$id[node_data$window == 2],
#                                 size = ceiling(length(which(node_data$window == 2))/3),
#                                 replace = F),
#                          '_2'),
#                   paste0(sample(x = node_data$id[node_data$window == 3],
#                                 size = ceiling(length(which(node_data$window == 3))/3),
#                                 replace = F),
#                          '_3'))
# node_data <- node_data %>%
#   mutate(id_window = paste0(id, '_', window)) %>%
#   filter(id_window %in% sample_nodes)
# edge_summary <- edge_summary %>%
#   select(-sri, -median, -`2.5%`, -`97.5%`) %>%
#   mutate(id1_window = paste0(id_1, '_', window),
#          id2_window = paste0(id_2, '_', window)) %>%
#   filter(id1_window %in% sample_nodes,
#          id2_window %in% sample_nodes)

## reset node_random values so doesn't screw with indexing
nodes_random <- node_data %>%
  select(node, id) %>%
  distinct() %>%
  mutate(id_1 = id, id_2 = id,
         node_1 = node, node_2 = node)
nodes_random$node_1_random <- nodes_random$node_2_random <- sample(1:nrow(nodes_random),
                                                                   replace = F, nrow(nodes_random))
if('node_1_random' %in% colnames(edge_summary)){
  edge_summary <- edge_summary %>%
    select(-node_1_random, -node_2_random) %>%
    left_join(nodes_random[,c('id_1','node_1','node_1_random')], by = c('node_1','id_1')) %>%
    left_join(nodes_random[,c('id_2','node_2','node_2_random')], by = c('node_2','id_2'))
} else {
  edge_summary <- edge_summary %>%
    left_join(nodes_random[,c('id_1','node_1','node_1_random')], by = c('node_1','id_1')) %>%
    left_join(nodes_random[,c('id_2','node_2','node_2_random')], by = c('node_2','id_2'))
}

#### simulate age effect ####
rm(list = ls()[! ls() %in% c('node_data','edge_summary')]) ; gc()

## identify global parameters
(n_data <- nrow(edge_summary))
(n_dyads <- length(unique(edge_summary$dyad_random)))
(n_nodes <- length(unique(node_data$id)))
(n_windows <- length(unique(edge_summary$window)))

## set "true" values
sim_diff <- -0.8
sim_mean <- 0.2
sim_int <- -5
sim_node_effect <- rnorm(n = n_nodes, mean = 0, sd = 0.4)
sim_window_effect <- rnorm(n = n_windows, mean = 0, sd = 0.6)
sim_sigma <- 0.5
print('true parameters set')

# ## standardise age values -- standardise together so that age of M1 when younger == age of M1 when older in the same time window
# mean_age <- mean(node_data$age)
# stdv_age <- sd(node_data$age)
#
# edge_summary <- edge_summary %>%
#   rename(window = period) %>%
#   mutate(age_min_std = (age_min - mean_age) / stdv_age,
#          age_max_std = (age_max - mean_age) / stdv_age)
#

## calculate age average and difference
edge_summary$age_diff <- abs(edge_summary$age_start_1 - edge_summary$age_start_2)
edge_summary$age_mean <- NA
for(i in 1:nrow(edge_summary)){
  edge_summary$age_mean[i] <- mean(c(edge_summary$age_start_1[i], edge_summary$age_start_2[i]))
}
edge_summary$age_diff_std <- (edge_summary$age_diff - mean(edge_summary$age_diff)) / sd(edge_summary$age_diff)
edge_summary$age_mean_std <- (edge_summary$age_mean - mean(edge_summary$age_mean)) / sd(edge_summary$age_mean)
print('predictor values calculated')

## set mean edge weight
str(edge_summary)
edge_summary$mu <- sim_int + edge_summary$age_diff_std * sim_diff + edge_summary$age_mean_std * sim_mean + sim_node_effect[edge_summary$node_1_random] + sim_node_effect[edge_summary$node_2_random] + sim_window_effect[edge_summary$window] # simulate mean centrality on unstandardised scale
ggplot(edge_summary)+
  geom_point(aes(x = age_diff, y = mu, colour = age_mean))
print('mean edge simulated')

## standardise edges
edge_summary$mu_std <- ( edge_summary$mu - mean(edge_summary$mu) ) / sd(edge_summary$mu)
ggplot(edge_summary)+
  geom_point(aes(x = age_diff, y = mu_std, colour = age_mean))
print('mean edge standardised')

## simulate full distribution of samples per node
sim_dat <- matrix(data = NA, nrow = 1000, ncol = n_data,
                  dimnames = list(1:1000, edge_summary$dyad_window))    # create matrix
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = edge_summary$mu_std[j], sd = sim_sigma)  # simulate distribution
}
plot(sim_dat[1,] ~ edge_summary$age_mean)          # plot simulated values against minimum age
plot(sim_dat[1,] ~ edge_summary$age_diff)          # plot simulated values against maximum age

## print progress marker
print('edges simulated')
save.image('step5_dyadicregression/anp_dyadic_simulation_diffmean_noncentredsigma.RData')

#### plot raw data ####
ggplot()+
  geom_point(data = edge_summary, aes(x = age_mean_std,
                             y = mu_std,
                             colour = age_diff),
             shape = 19)+
  scale_x_continuous('mean age of dyad members (standardised)')+
  scale_y_continuous('mean estimated edge weight (standardised)')+
  scale_colour_viridis_c()+
  labs(colour = 'age difference between dyad')

# sim_dat_long <- sim_dat %>%
#   as.data.frame() %>%
#   pivot_longer(cols = everything(), names_to = 'dyad_window', values_to = 'edge_draw') %>%
#   left_join(edge_summary, by = 'dyad_window')
# ggplot()+
#   geom_violin(data = sim_dat_long,
#               aes(x = as.factor(age_cat_mean),
#                   y = edge_draw,
#                   fill = as.factor(age_diff)),
#               #shape = 19,
#               alpha = 0.5)+
#   geom_point(data = edge_summary,
#              aes(x = as.factor(age_cat_mean), y = mu_std,
#                  group = as.factor(age_diff)),
#              shape = 19,
#              #colour = 'white',
#              size = 1)+
#   scale_x_discrete('age category of younger dyad member')+
#   scale_y_continuous('mean estimated edge weight')+
#   scale_fill_viridis_d()+
#   theme(legend.position = 'bottom')+
#   labs(fill = 'age category of older elephant')

## print progress marker
print('edges plotted')

#### fit Gaussian distribution to output of edge weight model ####
logit_edge_draws_mu <- apply(sim_dat, 2, mean)
logit_edge_draws_sd <- apply(sim_dat, 2, sd)

## set up plotting
num_check <- 64
selected_samples <- sample(1:n_dyads, num_check, replace = FALSE)

### Setting grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

### plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- logit_edge_draws_sd[i] #sqrt(logit_edge_draws_cov[i,i])

  fitted_values <- rnorm(1e5, mean = mu, sd = sd)

  hist(unlist(sim_dat[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="red", lw=1.5)
}

for (i in selected_samples) {
  plot(unlist(sim_dat[,i]), unlist(sim_dat[,i+1]),
       col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

## reset plot window and clean up
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))
rm(cols, fitted_values, i, j, mu, num_check, rows, sd, selected_samples) ; gc()

## print progress marker
print('normal approximation complete')

save.image('step5_dyadicregression/anp_dyadic_simulation_dataprep_allwindows.RData')

#################### Age difference ####################
# load('step5_dyadicregression/anp_dyadic_simulation_dataprep_allwindows.RData')
dev.off()
pdf('step5_dyadicregression/simulate_dyadic_anp_diff_noncentredsigma.pdf')

#### prior predictive check ####
n <- 100
intercept <- rnorm(n, 0, 0.5)
beta_age_diff <- rnorm(n, 0, 1)
tau_sigma_raw <- rnorm(n, 0, 2)

# mu_dyad <- rnorm(n, 0 , 0.4)
# rand_dyad <- rnorm(n, 0 ,1)
# sigma_dyad <- rexp(n, 2)

mu_window <- rnorm(n, 0, 0.4)
rand_window <- rnorm(n, 0, 1)
sigma_window <- rexp(n,2)

mu_mm <- rnorm(n, 0, 0.5)
rand_mm <- rnorm(n, 0, 1)
tau_mm <- rcauchy(n, 0, 0.5)
raw_sigma <- rnorm(n, 0, 1)

# dyad_random_effect <- mu_dyad + rand_dyad * sigma_dyad
window_random_effect <- mu_window + rand_window * sigma_window

node_mean <- mu_mm + rand_mm * tau_mm
# node_sigma <- sigma_mm * theta_mm

tau_sigma <- exp(tau_sigma_raw)
node_sigma <- exp(raw_sigma)

mm_nodes <- matrix(NA, nrow = n, ncol = n)
for(i in 1:ncol(mm_nodes)){
  mm_nodes[,i] <- rnorm(n, mean = node_mean[i], sd = node_sigma[i])
}

plot(NULL, las = 1, xlab = 'age difference', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-15, max(logit_edge_draws_mu)+15),
     xlim = c(min(edge_summary$age_diff_std), max(edge_summary$age_diff_std)))
abline(h = min(logit_edge_draws_mu), lty = 2)
abline(h = max(logit_edge_draws_mu), lty = 2)

x <- c(min(edge_summary$age_diff_std), max(edge_summary$age_diff_std))
for(i in 1:n){
  y <- rep(NA, length(x))
  nodes_sample <- sample(1:n, size = 2, replace = F)
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age_diff[i]*x[j] +
      window_random_effect[i] + #dyad_random_effect[i] +
      mm_nodes[i,nodes_sample[1]] + mm_nodes[i,nodes_sample[2]]
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}

rm(intercept, beta_age_diff, tau_sigma_raw, mu_window, rand_window, sigma_window, mu_mm, rand_mm, tau_mm, raw_sigma, window_random_effect, node_mean, tau_sigma, node_sigma, mm_nodes, i, j, n, x, y) ; gc()

## print progress marker
print('prior predictive complete')
save.image('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')

#### fit dyadic regression ####
#load('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')

## create data list
dyad_data <- list(
  # global data size
  num_data = n_data,                         # number of edges measured
  num_dyads = n_dyads,                       # number of dyads
  num_nodes = n_nodes,                       # number of nodes
  num_windows = n_windows,                   # number of time windows measured

  # Gaussian approximation of edge weights
  logit_edge_mu = logit_edge_draws_mu,       # Means of Gaussian approximation of logit edge weights
  logit_edge_sd = logit_edge_draws_sd,       # SD of Gaussian approximation of logit edge weights

  # explanatory variables
  age_diff = edge_summary$age_diff_std,      # age of younger dyad member

  # multimembership terms
  node_1 = edge_summary$node_1_random,       # Node 1 IDs for multimembership terms
  node_2 = edge_summary$node_2_random,       # Node 2 IDs for multimembership terms

  # time window
  window = edge_summary$window,              # Window ID random effect

  # dyad id
  dyad = edge_summary$dyad_random
)

## print progress marker
print('data list created')

## load dyadic regression model
dyadic_regression <- stan_model('models/dyadic_regression_anp_agediff_noncentredsigma.stan')
#dyadic_regression <- cmdstan_model('models/dyadic_regression_anp_agediff_noncentredsigma.stan')

## print progress marker
print('model loaded in')

## define parameters
n_chains <- 4
n_samples <- 2500

## fit dyadic regression
print('start model')
fit_dyadreg_simdiff <- sampling(dyadic_regression,
                                data = dyad_data, 
                                iter = n_samples*3, warmup = n_samples*2,
                                chains = n_chains, cores = n_chains,
                                control=list(adapt_delta=0.9,
                                             max_treedepth = 15)
)
# fit_dyadreg_simdiff <- dyadic_regression$sample(
#   data = dyad_data,
#   iter_warmup = n_samples,
#   iter_sampling = n_samples,
#   chains = n_chains,
#   parallel_chains = n_chains)

save.image('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')

## print progress marker
print('model run')

#### check outputs ####
load('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')

## extract model fit
# summary <- fit_dyadreg_simdiff$summary()
summary <- rstan::summary(fit_dyadreg_simdiff)
(summary <- as.data.frame(summary$summary))
par(mfrow = c(3,1))
# hist(summary$rhat, breaks = 50)
# hist(summary$ess_bulk, breaks = 50)
# hist(summary$ess_tail, breaks = 50)
hist(summary$Rhat, breaks = 50)
hist(summary$n_eff, breaks = 50)
par(mfrow = c(1,1))

## extract draws
# draws <- fit_dyadreg_simdiff$draws(format = 'df')
draws <- rstan::extract(fit_dyadreg_simdiff)

## extract dyadic regression parameters
b_diff <- draws$beta_age_diff
intercept <- draws$intercept
sigma <- draws$sigma
parameters <- data.frame(beta_age_diff = b_diff,
                         intercept = intercept,
                         sigma = sigma) %>%
  mutate(chain = rep(1:n_chains, each = n_samples),
         position = rep(1:n_samples, n_chains)) %>%
  pivot_longer(cols = c('beta_age_diff','sigma','intercept'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplot function
traceplot <- function(parameter_data){
  plot <- ggplot(data = parameter_data)+
    geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
    theme(legend.position = 'none')+
    scale_colour_viridis_d()+
    facet_wrap(. ~ parameter , scales = 'free_y')
  return(plot)
}

## global parameter traceplots
# ggplot(data = parameters)+
#   geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
#   theme(legend.position = 'none')+
#   scale_colour_viridis_d()+
#   facet_wrap(. ~ parameter , scales = 'free_y')
traceplot(parameter_data = parameters)

print('global parameters plotted')

## random effect extraction function
long_random <- function(random_effect, n){
  random_effect_long <- random_effect %>%
    mutate(chain = rep(1:n_chains, each = n_samples),
           position = rep(1:n_samples, n_chains)) %>%
    pivot_longer(cols = all_of(1:n),
                 names_to = 'parameter',
                 values_to = 'slope_draw')
  return(random_effect_long)
}

## extract multimembership draws
mm_nodes <- as.data.frame(draws$mm_nodes)
colnames(mm_nodes) <- unique(c(dyad_data$node_1, dyad_data$node_2)) # I don't know that this is right, but I can't see any other order it would be in?
mm_nodes_long <- long_random(mm_nodes, n_nodes)

## multimembership traceplots
mm_nodes_long %>%
  filter(parameter %in% sample(unique(mm_nodes_long$parameter), 25, replace = F)) %>%
  traceplot()

traceplot2 <- function(parameter_data){
  plot <- ggplot(data = parameter_data)+
    geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
    theme(legend.position = 'none')+
    scale_colour_viridis_d()+
    facet_grid(parameter ~ as.factor(chain), scales = 'free_y')
  return(plot)
}
mm_nodes_long %>%
  filter(parameter %in% sample(unique(mm_nodes_long$parameter), 5, replace = F)) %>%
  traceplot2()

print('multimembership parameters plotted')

## extract window random effect draws
rand_window <- as.data.frame(draws$window_random_effect)
colnames(rand_window) <- unique(dyad_data$window)
rand_window_long <- long_random(rand_window, n_windows)

## window random effect traceplots
traceplot(rand_window_long)
traceplot2(rand_window_long)

print('window parameters plotted')

# # ## extract dyad random effect draws
# # rand_dyad <- as.data.frame(draws$dyad_random_effect)
# # colnames(rand_dyad) <- unique(dyad_data$dyad)
# # rand_dyad_long <- long_random(rand_dyad, n_dyads)
# #
# # ## dyad random effect traceplots
# # rand_dyad_long %>%
# #   filter(parameter %in% sample(unique(rand_dyad_long$parameter), 25, replace = F)) %>%
# #   traceplot()

## print progress marker
print('outputs checked')

#### posterior predictive check ####
plot(density(as.numeric(sim_dat[1, ])),
     main = "Posterior predictive density of edge weights:\nblack = measured edge, colours = predicted edge (1 colour per chain)",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:nrow(sim_dat), 1)

  # mu_plot <- rep(NA, n_dyads)
  # for(k in 1:n_dyads){
  #   mu_plot[k] <- intercept[j] + b_diff[j]*dyad_data$age_diff[k] + mm_nodes[j,dyad_data$node_1[k]] + mm_nodes[j,dyad_data$node_2[k]]+ rand_window[j,dyad_data$window[k]]# + rand_dyad[j,dyad_data$dyad_id[k]]
  # }
  # 
  # sigma_plot <- dyad_data$logit_edge_sd + sigma[j]
  # norm <- rnorm(1000, mu_plot, sigma_plot)
  # 
  # lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  # lines(density(norm), col = rgb(1, 0, 0, 0.25))                     # red lines for predictions
  
  mu_plot <- matrix(NA, ncol = n_dyads, nrow = n_chains)
  
  for(l in c(1:n_chains)){
    n <- (l*n_samples) - n_samples
    for(k in 1:n_dyads){
      mu_plot[l,k] <- intercept[j+n] + 
        b_diff[j+n]*dyad_data$age_diff[k] + 
        mm_nodes[j+n,dyad_data$node_1[k]] + 
        mm_nodes[j+n,dyad_data$node_2[k]] + 
        #rand_dyad[j+n,dyad_data$dyad_id[k]] +
        rand_window[j+n,dyad_data$window[k]]
    }
    sigma_plot <- dyad_data$logit_edge_sd + diag(rep(sigma[j+n], n_data))
    norm <- rnorm(1000, mu_plot[l,], sigma_plot)
    lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25))
    lines(density(norm),
          col = ifelse(l == 1, rgb(1, 0, 0, 0.25),
                       ifelse(l == 2, rgb(0, 1, 0, 0.25),
                              ifelse(l == 3, rgb(0, 0, 1, 0.25),
                                     rgb(1, 1, 0, 0.25)))) )
  }
}

save.image('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')
print('posterior predictive check complete')

# #### predict from raw data ####
# load('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')
# 
# # ## create predictive data frame
# # age_diff_length <- length(unique(dyad_data$age_diff))
# # pred <- data.frame(age_diff_std = rep(unique(dyad_data$age_diff),
# #                                       n_chains*n_samples),
# #                    age_diff = rep(unique(edge_summary$age_diff),
# #                                   n_chains*n_samples),
# #                    intcp = rep(intercept,
# #                                each = age_diff_length),
# #                    beta_diff = rep(b_diff,
# #                                    each = age_diff_length))
# # 
# # ## calculate predictions
# # pred$pred_mu <- pred$intcp + pred$beta_diff * pred$age_diff_std
# 
# ## create function for predicting mean
# get_mean_predictions <- function(prediction_df, draws_list){
#   mu <- matrix(NA, nrow = length(draws_list$intercept), ncol = nrow(prediction_df),
#                dimnames = list(1:length(draws_list$intercept),
#                                paste0('dyad',prediction_df$dyad_random)))
# 
#   multimembership <- draws_list$mm_nodes
#   # window_random <- draws_list$rand_window
#   # #dyad_random <- draws_list$rand_dyad
#   intcp <- draws_list$intercept
#   b_diff <- draws_list$beta_age_diff
# 
#   for(dyad in 1:ncol(mu)){
#     node_1 <- multimembership[,prediction_df$node_1_random[dyad]]
#     node_2 <- multimembership[,prediction_df$node_2_random[dyad]]
#     # window_effect <- window_random[,prediction_df$window[dyad]]
#     # #dyad_effect <- dyad_random[,prediction_df$dyad_random[dyad]]
# 
#     for(draw in 1:nrow(mu)){
#       mu[draw, dyad] <- intcp[draw] + b_diff[draw] * prediction_df$age_diff_std[dyad] +
#         node_1[draw] + node_2[draw]# + dyad_effect[draw]# + window_effect[draw]
#     }
#   }
#   return(mu)
# }
# 
# ## predict mean distribution
# pred_mean <- get_mean_predictions(prediction_df = edge_summary,
#                                   draws_list = draws)
# 
# ## summarise mean distribution
# edge_summary$pred_mu <- apply(pred_mean, 2, mean)
# edge_summary$pred_mean_lwr <- apply(pred_mean, 2, quantile, prob = 0.025)
# edge_summary$pred_mean_upr <- apply(pred_mean, 2, quantile, prob = 0.975)
# 
# ## create function for predicting full distribution
# get_full_predictions <- function(mu_matrix, sigma){
#   predictions <- mu_matrix
#   for(dyad in 1:ncol(predictions)){
#     for(draw in 1:nrow(predictions)){
#       predictions[draw, dyad] <- rnorm(n = 1,
#                                        mean = mu_matrix[draw, dyad],
#                                        sd = sigma[draw])
#     }
#   }
#   return(predictions)
# }
# 
# ## predict full distribution
# pred_full <- get_full_predictions(mu_matrix = pred_mean,
#                                   sigma = draws$sigma)
# 
# ## summarise full distribution
# edge_summary$pred_full_lwr <- apply(pred_full, 2, quantile, prob = 0.025)
# edge_summary$pred_full_upr <- apply(pred_full, 2, quantile, prob = 0.975)
# 
# ## convert to unstandardised scale
# pred_mean_ustd <- pred_mean * sd(edge_summary$mu) + mean(edge_summary$mu)
# pred_full_ustd <- pred_full * sd(edge_summary$mu) + mean(edge_summary$mu)
# 
# ## summarise unstandardised distributions
# edge_summary$pred_mu_ustd <- apply(pred_mean_ustd, 2, mean)
# edge_summary$pred_mean_lwr_ustd <- apply(pred_mean_ustd, 2, quantile, prob = 0.025)
# edge_summary$pred_mean_upr_ustd <- apply(pred_mean_ustd, 2, quantile, prob = 0.975)
# edge_summary$pred_full_lwr_ustd <- apply(pred_full_ustd, 2, quantile, prob = 0.025)
# edge_summary$pred_full_upr_ustd <- apply(pred_full_ustd, 2, quantile, prob = 0.975)
# 
# ## convert to invlogit scale
# pred_mean_invlogit <- invlogit(pred_mean_ustd)
# pred_full_invlogit <- invlogit(pred_full_ustd)
# 
# ## summarise invlogit distributions
# edge_summary$pred_mu_invlogit <- apply(pred_mean_invlogit, 2, mean)
# edge_summary$pred_mean_lwr_invlogit <- apply(pred_mean_invlogit, 2, quantile, prob = 0.025)
# edge_summary$pred_mean_upr_invlogit <- apply(pred_mean_invlogit, 2, quantile, prob = 0.975)
# edge_summary$pred_full_lwr_invlogit <- apply(pred_full_invlogit, 2, quantile, prob = 0.025)
# edge_summary$pred_full_upr_invlogit <- apply(pred_full_invlogit, 2, quantile, prob = 0.975)
# 
# # ## summarise
# # pred_summary <- pred %>%
# #   group_by(age_diff_std, age_mean_std) %>%
# #   mutate(pred_lwr = quantile(pred_mu, probs = 0.025),
# #          pred_mean = mean(pred_mu),
# #          pred_upr = quantile(pred_mu, probs = 0.975),
# #          pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
# #          pred_unstd_mean = mean(pred_unstd),
# #          pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>%
# #   ungroup() %>%
# #   select(age_diff_std, age_mean_std,
# #          pred_lwr, pred_mean, pred_upr,
# #          pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>%
# #   distinct()
# 
# #### plot predictions ####
# # pred_summary <- pred_summary %>%
# #   left_join(distinct(edge_summary[,c('age_diff','age_diff_std')]),
# #             by = 'age_diff_std') %>%
# #   left_join(distinct(edge_summary[,c('age_max','age_max_std')]),
# #             by = 'age_max_std')
# ggplot()+
#   geom_ribbon(data = edge_summary,
#               aes(x = age_diff,
#                   ymin = pred_full_lwr, ymax = pred_full_upr,
#                   fill = as.factor(window)),
#               alpha = 0.3)+
#   geom_ribbon(data = edge_summary,
#               aes(x = age_diff,
#                   ymin = pred_mean_lwr, ymax = pred_mean_upr,
#                   fill = as.factor(window)),
#               alpha = 0.3)+
#   geom_line(data = edge_summary,
#             aes(x = age_diff,
#                 y = pred_mu,
#                 colour = as.factor(window)),
#             linewidth = 1)+
#   scale_colour_viridis_d(direction = -1)+ scale_fill_viridis_d(direction = -1)+
#   geom_point(data = edge_summary,
#              aes(x = age_diff,
#                  y = mu_std,
#                  colour = as.factor(window)),
#              alpha = 0.2
#   )+
#   scale_x_continuous('age difference between dyad members')+
#   scale_y_continuous('mean estimated edge weight')+
#   facet_wrap(. ~ window)+
#   theme(axis.text = element_text(size = 14),
#         legend.position = 'bottom', #c(0.8,0.2),
#         # legend.title = element_text(size = 14),
#         # legend.text = element_text(size = 12)
#         axis.title = element_text(size = 18)) #+
#   # guides(colour = guide_legend(ncol = 6),
#   #        fill = guide_legend(ncol = 6))
# 
# # raw_summary <- edge_summary %>%
# #   select(age_diff, mu, mu_std) %>%
# #   group_by(age_diff) %>%
# #   mutate(mu_mean = mean(mu),
# #          mu_std_mean = mean(mu_std)) %>%
# #   ungroup() %>%
# #   select(-mu, -mu_std) %>%
# #   distinct()
# # compare <- pred_summary %>%
# #   left_join(raw_summary[,c('mu_mean','mu_std_mean','age_pair')], by = 'age_diff') %>%
# #   rename(raw_unstd = mu_mean, raw_std = mu_std_mean) %>%
# #   select(age_diff, raw_unstd, pred_unstd_mean, raw_std, pred_mean, pred_lwr, pred_unstd_lwr, pred_upr, pred_unstd_upr)
# # ggplot(compare)+
# #   geom_ribbon(aes(x = raw_unstd,
# #                   ymin = pred_unstd_lwr, ymax = pred_unstd_upr),
# #               alpha = 0.3)+
# #   geom_line(aes(x = raw_unstd, y = pred_unstd_mean),
# #             linewidth = 1)+
# #   geom_line(data = data.frame(x = c(min(edge_summary$mu),max(edge_summary$mu)),
# #                               y = c(min(edge_summary$mu),max(edge_summary$mu))),
# #             aes(x = x, y = y),
# #             linewidth = 0.5, colour = 'red')+
# #   scale_x_continuous('raw mean')+
# #   scale_y_continuous('predicted mean')+
# #   theme(axis.text = element_text(size = 14),
# #         axis.title = element_text(size = 18),
# #         legend.title = element_text(size = 14),
# #         legend.text = element_text(size = 12))
# 
# ggplot()+
#   geom_point(data = edge_summary,
#              mapping = aes(x = mu_std, y = pred_mu))+
#   geom_abline(slope = 1, intercept = 0)+
#   facet_wrap(. ~ window)
# 
# ggplot()+
#   geom_point(data = edge_summary,
#              mapping = aes(x = mu, y = pred_mu_ustd))+
#   geom_abline(slope = 1, intercept = 0)
# 
# ggplot()+
#   geom_point(data = edge_summary,
#              mapping = aes(x = invlogit(mu), y = pred_mu_invlogit))+
#   geom_abline(slope = 1, intercept = 0)
# 
# dev.off()
# save.image('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')
# 
# ## print progress marker
# print('predictions complete')
# 
# #### extract original parameters ####
# load('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')
# 
# ## create altered dataframe to predict from
# pred_new <- edge_summary %>%
#   rename(age_diff_org = age_diff,
#          age_diff_org_std = age_diff_std) %>%
#   mutate(age_diff_new = age_diff_org)
# # pred_new_diff <- pred_new %>%
# #   mutate(age_diff_new = age_diff_new + 1) %>%
# #   mutate(age_diff_new_std = (age_diff_new - mean(edge_summary$age_diff))/sd(edge_summary$age_diff))
# pred_new_diff <- pred_new %>%
#   mutate(age_diff_new_std = age_diff_org_std + 1) %>%
#   mutate(age_diff_new = age_diff_new_std * sd(edge_summary$age_diff) + mean(edge_summary$age_diff))
# 
# ## adapt function for predicting mean
# get_mean_predictions <- function(prediction_df, draws_list, new){
#   mu <- matrix(NA, nrow = length(draws_list$intercept), ncol = nrow(prediction_df),
#                dimnames = list(1:length(draws_list$intercept),
#                                paste0('dyad',prediction_df$dyad_random)))
# 
#   multimembership <- draws_list$mm_nodes
#   # window_random <- draws_list$rand_window
#   # #dyad_random <- draws_list$rand_dyad
#   intcp <- draws_list$intercept
#   b_diff <- draws_list$beta_age_diff
# 
#   if(is.null(new) == FALSE){
#     if(new == TRUE){
#       prediction_df$age_diff <- prediction_df$age_diff_new
#       prediction_df$age_diff_std <- prediction_df$age_diff_new_std
#     }
#     if(new == FALSE){
#       prediction_df$age_diff <- prediction_df$age_diff_org
#       prediction_df$age_diff_std <- prediction_df$age_diff_org_std
#     }
# 
#   }
# 
#   for(dyad in 1:ncol(mu)){
#     node_1 <- multimembership[,prediction_df$node_1_random[dyad]]
#     node_2 <- multimembership[,prediction_df$node_2_random[dyad]]
#     # window_effect <- window_random[,prediction_df$window[dyad]]
#     # #dyad_effect <- dyad_random[,prediction_df$dyad_random[dyad]]
# 
#     for(draw in 1:nrow(mu)){
#       mu[draw, dyad] <- intcp[draw] + b_diff[draw] * prediction_df$age_diff_std[dyad] +  node_1[draw] + node_2[draw]# + window_effect[draw]# + dyad_effect[draw]
#         #
#     }
#   }
#   return(mu)
# }
# 
# ## predict mean distributions
# pred_org_diff_mean <- get_mean_predictions(prediction_df = pred_new_diff,
#                                            draws_list = draws,
#                                            new = FALSE)
# pred_new_diff_mean <- get_mean_predictions(prediction_df = pred_new_diff,
#                                            draws_list = draws,
#                                            new = TRUE)
# 
# ## predict full distribution
# pred_org_diff_full <- get_full_predictions(mu_matrix = pred_org_diff_mean,
#                                            sigma = draws$sigma)
# pred_new_diff_full <- get_full_predictions(mu_matrix = pred_new_diff_mean,
#                                            sigma = draws$sigma)
# 
# ## convert to unstandardised scale
# pred_org_diff_full_ustd <- pred_org_diff_full * sd(edge_summary$mu) + mean(edge_summary$mu)
# pred_new_diff_full_ustd <- pred_new_diff_full * sd(edge_summary$mu) + mean(edge_summary$mu)
# 
# ## calculate contrasts: age difference increases by 1 stdv, standardised scale
# contrast_diff <- pred_new_diff_full - pred_org_diff_full
# print(paste0('contrast for age_diff + 1 (std): ',
#              round(mean(contrast_diff), 5), ' ± ', round(sd(contrast_diff), 5),
#              ' [', round(quantile(contrast_diff, prob = 0.025), 5),
#              ' - ',round(quantile(contrast_diff, prob = 0.975), 5), ']'))
# 
# ## calculate contrasts: minimum age increases by 1, outcome scale
# contrast_diff_ustd <- pred_new_diff_full_ustd - pred_org_diff_full_ustd
# print(paste0('contrast for age_diff + 1 (ustd): ',
#              round(mean(contrast_diff_ustd), 5), ' ± ', round(sd(contrast_diff_ustd), 5),
#              ' [', round(quantile(contrast_diff_ustd, prob = 0.025), 5),
#              ' - ',round(quantile(contrast_diff_ustd, prob = 0.975), 5), ']'))
# print(paste0('true age diff effect: ', sim_diff))
# 
# ## save workspace
# save.image('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')
# 
# #### final clean plots ####
# load('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')
# 
# ## create counterfactual data frame for plotting from: minimum age on x
# fake_data <- expand.grid(window = 1:n_windows,
#                          dyad_random = unique(edge_summary$dyad_random),
#                          age_diff = min(edge_summary$age_diff):max(edge_summary$age_diff)) %>% 
#   left_join(distinct(edge_summary[,c('dyad_random','node_1','node_2')]),
#             by = 'dyad_random') %>%
#   mutate(age_diff_std = (age_diff - mean(edge_summary$age_diff)) / sd(edge_summary$age_diff))
# 
# ## predict mean values
# fake_mean <- get_mean_predictions(prediction_df = fake_data, draws_list = draws, new = NULL)
# length(which(is.na(fake_mean) == TRUE))
# print('mean values extracted')
# 
# ## predict full distribution
# fake_full <- get_full_predictions(mu_matrix = fake_mean, sigma = draws$sigma)
# warnings()
# print('full values extracted')
# length(which(is.na(fake_full) == TRUE))
# 
# ## summarise values
# fake_data$mean_prediction <- apply(fake_mean, 2, mean, na.rm = T)
# fake_data$mu_lwr_prediction <- apply(fake_mean, 2, quantile, prob = 0.025, na.rm = T)
# fake_data$mu_upr_prediction <- apply(fake_mean, 2, quantile, prob = 0.975, na.rm = T)
# fake_data$full_lwr_prediction <- apply(fake_full, 2, quantile, prob = 0.025, na.rm = T)
# fake_data$full_upr_prediction <- apply(fake_full, 2, quantile, prob = 0.975, na.rm = T)
# print('predictions summarised')
# save.image('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')
# 
# ## summarise summary
# summary_fake <- fake_data %>% 
#   group_by(age_diff) %>% 
#   summarise(mean = mean(mean_prediction),
#             mu_lwr = mean(mu_lwr_prediction),
#             mu_upr = mean(mu_upr_prediction),
#             full_lwr = mean(full_lwr_prediction),
#             full_upr = mean(full_upr_prediction))
# print('summary summarised')
# 
# ## plot altogether
# (diff_overall <- ggplot()+
#     geom_ribbon(data = summary_fake,
#                 mapping = aes(x = age_diff,
#                               ymin = full_lwr, ymax = full_upr),
#                 alpha = 0.2)+
#     geom_ribbon(data = summary_fake,
#                 mapping = aes(x = age_diff,
#                               ymin = mu_lwr, ymax = mu_upr),
#                 alpha = 0.4)+
#     geom_line(data = summary_fake,
#               mapping = aes(x = age_diff,
#                             y = mean))+
#     geom_point(data = edge_summary,
#                mapping = aes(x = age_diff,
#                              y = mu_std),
#                alpha = 0.2)+
#     scale_colour_viridis_c(begin = (1/8), end = 1)+
#     scale_fill_viridis_c(begin = (1/8), end = 1)+
#     theme(legend.position = 'bottom')+
#     labs(fill = 'age of older male',
#          colour = 'age of older male'))
# ggsave(plot = diff_overall, device = 'png',
#        filename = 'anp_dyadic_simulation_diffageeffect_overall.png',
#        path = 'step5_dyadicregression/',
#        width = 1200, height = 800, units = 'px')
# print('overall plot completed with minimum age')
# 
# ## resummarise summary
# summary_fake <- fake_data %>% 
#   group_by(age_diff, window) %>% 
#   summarise(mean = mean(mean_prediction),
#             mu_lwr = mean(mu_lwr_prediction),
#             mu_upr = mean(mu_upr_prediction),
#             full_lwr = mean(full_lwr_prediction),
#             full_upr = mean(full_upr_prediction))
# 
# ## plot separately by window
# (diff_bywindow <- ggplot()+
#     geom_ribbon(data = summary_fake,
#                 mapping = aes(x = age_diff,
#                               ymin = full_lwr, ymax = full_upr),
#                 alpha = 0.2)+
#     geom_ribbon(data = summary_fake,
#                 mapping = aes(x = age_diff,
#                               ymin = mu_lwr, ymax = mu_upr),
#                 alpha = 0.4)+
#     geom_line(data = summary_fake,
#               mapping = aes(x = age_diff,
#                             y = mean))+
#     geom_point(data = edge_summary,
#                mapping = aes(x = age_diff,
#                              y = mu_std),
#                alpha = 0.2)+
#     scale_colour_viridis_c(begin = (1/8), end = 1)+
#     scale_fill_viridis_c(begin = (1/8), end = 1)+
#     theme(legend.position = 'bottom')+
#     labs(fill = 'age of older male',
#          colour = 'age of older male')+
#     facet_wrap(. ~ window))
# ggsave(plot = diff_bywindow, device = 'png',
#        filename = 'anp_dyadic_simulation_diffageeffect_bywindow.png',
#        path = 'step5_dyadicregression/',
#        width = 1600, height = 2400, units = 'px')
# print('faceted plot completed with minimum age')
# 
# ## save output
# save.image('step5_dyadicregression/anp_dyadic_simulation_agediff_noncentredsigma.RData')
# dev.off()
# print('age difference complete')
# 
