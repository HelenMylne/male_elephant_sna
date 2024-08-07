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

pdf('step5_dyadicregression/simulate_dyadic_anp_tighterpriors.pdf')

set.seed(12345)

#### load raw data as basis for simulation -- simulating using first 3 short windows only for sake of speed and computing power ####
# ## load time window 1 and remove additional data
# load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
# rm(list = ls()[! ls() %in% c('nodes','cdf_1','edges','edge_samples')]) ; gc()
#
# ## create full data frame with window variable
# node_data <- nodes %>%
#   mutate(window = 1)
# edge_summary <- cdf_1 %>%
#   dplyr::select(-dyad_rank) %>%
#   mutate(window = 1)
# edges_all <- list()
# edges_all[[1]] <- edge_samples
#
# ## attach all other data frames
# rm(list = ls()[! ls() %in% c('node_data','edge_summary','edges_all')]) ; gc()
# n_windows <- 3
# for(time_window in 2:n_windows){
#   ## import workspace image for time window
#   load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))
#
#   ## ensure all have the same names: edge weights
#   if('edge_weights_matrix' %in% ls()){
#     edge_samples <- edge_weights_matrix
#   } else {
#     if('edge_weight_matrix' %in% ls()){
#       edge_samples <- edge_weight_matrix
#     }
#   }
#
#   ## ensure all have the same names: dyads data frame
#   if('counts_df' %in% ls()){
#     cdf <- counts_df
#   }
#
#   ## add window variable
#   nodes <- nodes %>%
#     mutate(window = time_window) %>%
#     dplyr::select(colnames(node_data))
#   cdf <- cdf %>%
#     mutate(window = time_window) %>%
#     dplyr::select(colnames(edge_summary))
#
#   ## attach to full data frame
#   node_data <- rbind(node_data, nodes)
#   edge_summary <- rbind(edge_summary, cdf)
#   edges_all[[time_window]] <- edges_all
#
#   ## print progress marker
#   print(time_window)
#
#   ## clean environment
#   rm(list = ls()[! ls() %in% c('node_data','edge_summary','edges_all','time_window')]) ; gc()
# }

load('step5_dyadicregression/anpshort_dyadicregression_dataimported.RData')
edge_summary <- edge_summary %>%
  filter(period < 4) %>%
  mutate(dyad_random = as.integer(as.factor(dyad_random))) # retain randomised order, but remove all of the missing dyads caused by filtering out windows 4:36
node_data <- node_data %>% filter(window < 4)
# edge_weights <- list(edge_weights[[1]],edge_weights[[2]],edge_weights[[3]])
# edge_samples <- edge_weights[[1]]
# colnames(edge_samples) <- edge_summary$dyad_window[edge_summary$period == 1]
# for(i in 1:n_windows){
#   edges <- edge_weights[[i]]
#   colnames(edges) <- edge_summary$dyad_window[edge_summary$period == i]
#   edge_samples <- cbind(edge_samples, edges)
# }

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
edge_summary <- edge_summary %>%
  select(-sri, -median, -`2.5%`, -`97.5%`)

## identify global parameters
n_data <- nrow(edge_summary)
n_dyads <- length(unique(edge_summary$dyad_random))
n_nodes <- length(unique(node_data$id))
n_windows <- length(unique(edge_summary$period))

## set "true" values
sim_min <- -0.2
sim_max <- 0.8
sim_int <- -5
sim_node_effect <- rnorm(n = n_nodes, mean = 0, sd = 0.4)
sim_window_effect <- rnorm(n = n_windows, mean = 0, sd = 0.6)
sim_sigma <- 0.5

## standardise age values -- standardise together so that age of M1 when younger == age of M1 when older in the same time window
mean_age <- mean(node_data$age)
stdv_age <- sd(node_data$age)

edge_summary <- edge_summary %>%
  rename(window = period) %>%
  mutate(age_min_std = (age_min - mean_age) / stdv_age,
         age_max_std = (age_max - mean_age) / stdv_age)

# edge_summary$age_min_std <- NA ; edge_summary$age_max_std <- NA
# node_ages <- node_data %>%
#   filter(window == 1) %>%
#   mutate(id_1 = id, id_2 = id,
#          node_1 = node, node_2 = node)
# node_ages$age2_std <- node_ages$age1_std <- (node_ages$age - mean(node_ages$age)) / sd(node_ages$age)
# for(time_window in 2:n_windows){
#   ages <- node_data %>%
#     filter(window == time_window) %>%
#     mutate(id_1 = id, id_2 = id,
#            node_1 = node, node_2 = node)
#   ages$age2_std <- ages$age1_std <- (ages$age - mean(ages$age)) / sd(ages$age)
#
#   node_ages <- rbind(node_ages, ages)
# }
# rm(ages) ; gc()
#
# edge_summary <- edge_summary %>%
#   left_join(node_ages[,c('id_1','node_1','age1_std','window')], by = c('id_1','node_1','window')) %>%
#   left_join(node_ages[,c('id_2','node_2','age2_std','window')], by = c('id_2','node_2','window')) %>%
#   mutate(age_min_std = NA, age_max_std = NA)
#
# for(i in 1:nrow(edge_summary)){
#   edge_summary$age_min_std[i] <- min(c(edge_summary$age1_std[i], edge_summary$age2_std[i]))
#   edge_summary$age_max_std[i] <- max(c(edge_summary$age1_std[i], edge_summary$age2_std[i]))
# }

## set mean edge weight
edge_summary$mu <- sim_int + edge_summary$age_min_std * sim_min + edge_summary$age_max_std * sim_max + sim_node_effect[edge_summary$node_1_random] + sim_node_effect[edge_summary$node_2_random] + sim_window_effect[edge_summary$window] # simulate mean centrality on normal scale
ggplot(edge_summary)+
  geom_point(aes(x = age_max, y = mu, colour = age_min))

## standardise edges
edge_summary$mu_std <- ( edge_summary$mu - mean(edge_summary$mu) ) / sd(edge_summary$mu)
ggplot(edge_summary)+
  geom_point(aes(x = age_max, y = mu_std, colour = age_min))

## simulate full distribution of samples per node
sim_dat <- matrix(data = NA, nrow = 1000, ncol = n_data,
                  dimnames = list(1:1000, edge_summary$dyad_window))    # create matrix
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = edge_summary$mu_std[j], sd = sim_sigma)  # simulate distribution
}
plot(sim_dat[1,] ~ edge_summary$age_min)          # plot simulated values against minimum age
plot(sim_dat[1,] ~ edge_summary$age_max)          # plot simulated values against maximum age

## print progress marker
print('edges simulated')
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

#### plot raw data ####
edge_summary$age_cat_max <- ifelse(edge_summary$age_max < 16, '10-15',
                                   ifelse(edge_summary$age_max < 20, '15-19',
                                          ifelse(edge_summary$age_max < 26, '21-25',
                                                 ifelse(edge_summary$age_max < 41, '25-40',
                                                        '40+'))))
edge_summary$age_cat_min <- ifelse(edge_summary$age_min < 16, '10-15',
                                   ifelse(edge_summary$age_min < 20, '15-19',
                                          ifelse(edge_summary$age_min < 26, '21-25',
                                                 ifelse(edge_summary$age_min < 41, '25-40',
                                                        '40+'))))
ggplot()+
  geom_point(data = edge_summary, aes(x = age_min_std,
                             y = mu_std,
                             colour = age_cat_max),
             shape = 19)+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  scale_colour_viridis_d()+
  labs(colour = 'age category of\nolder male')

ggplot()+
  geom_boxplot(data = edge_summary, aes(x = age_cat_min,
                               y = mu_std,
                               fill = age_cat_max),
               shape = 19)+
  scale_x_discrete('age category of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(legend.position = 'bottom')+
  scale_fill_viridis_d()+
  labs(fill = 'age category of older male')

sim_dat_long <- sim_dat %>%
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = 'dyad_window', values_to = 'edge_draw') %>%
  left_join(edge_summary, by = 'dyad_window')
ggplot()+
  geom_violin(data = sim_dat_long,
              aes(x = as.factor(age_cat_min),
                  y = edge_draw,
                  fill = as.factor(age_cat_max)),
              #shape = 19,
              alpha = 0.5)+
  geom_point(data = edge_summary,
             aes(x = as.factor(age_cat_min), y = mu_std,
                 group = as.factor(age_cat_max)),
             shape = 19,
             #colour = 'white',
             size = 1)+
  scale_x_discrete('age category of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  scale_fill_viridis_d()+
  theme(legend.position = 'bottom')+
  labs(fill = 'age category of older elephant')

## print progress marker
print('edges plotted')

#### fit multivariate Gaussian distribution to output of edge weight model ####
# ## calculate means and covariance matrix per window
# logit_edge_draws_mu <- list()
# logit_edge_draws_cov <- list()
# for(time_window in 1:n_windows){
#   window_edges <- sim_dat[,which(edge_summary$window == time_window)]
#   logit_edge_draws_mu[[time_window]] <- apply(window_edges, 2, mean)
#   # logit_edge_draws_cov[[time_window]] <- cov(window_edges)
#   print(paste0('mean calculated for window ',time_window, ' at ',Sys.time()))
# }
#
# ## parallelise making covariance matrix
# library(foreach)
# library(doParallel)
#
# #setup parallel backend to use many processors
# cores = detectCores()
# cl <- makeCluster(cores)
# registerDoParallel(cl)
#
# finalMatrix <- foreach(time_window = 1:n_windows,
#                        .combine=cbind) %dopar% {
#   tempMatrix = cov(sim_dat[,which(edge_summary$window == time_window)])
#   colnames(tempMatrix) = paste0(edge_summary$dyad_id[edge_summary$window == time_window],
#                                 '_',time_window)
#   tempMatrix
#                        }
#
# #stop cluster
# stopCluster(cl)
# print('multivariate Gaussian complete')

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

#### prior predictive check ####
n <- 100
beta_age_min <- rnorm(n, 0, 0.8)
beta_age_max <- rnorm(n, 0, 0.8)
intercept <- rnorm(n, 0, 2)
plot(NULL, las = 1, xlab = 'minimum age', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5),
     xlim = c(min(edge_summary$age_min_std), max(edge_summary$age_max_std)))
abline(h = min(logit_edge_draws_mu), lty = 2)
abline(h = max(logit_edge_draws_mu), lty = 2)
x <- c(min(edge_summary$age_min_std), max(edge_summary$age_max_std))
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age_min[i]*x[j] + beta_age_max[i]*x[j]
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age_max, beta_age_min, intercept, i, y, j, x) ; gc()

## print progress marker
print('prior predictive complete')
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

#### fit dyadic regression ####
load('step5_dyadicregression/anp_dyadic_simulation.RData')

## create data list
dyad_data <- list(
  # global data size
  num_data = n_data,                                         # number of edges measured
  num_dyads = n_dyads,                                       # number of dyads
  num_nodes = n_nodes,                                       # number of nodes
  num_windows = n_windows,                                   # number of time windows measured

  # # per time window data size
  # num_dyads_window1 = length(which(edge_summary$window == 1)),    # rows in data, window 1
  # num_dyads_window2 = length(which(edge_summary$window == 2)),    # rows in data, window 2
  # num_dyads_window3 = length(which(edge_summary$window == 3)),    # rows in data, window 3

  # # number of dyads in all preceding time windows for age indexing
  # num_dyads_prev_windows = c(0,                              # Number of rows in data in time window -1
  #                            length(which(edge_summary$window == 1)),
  #                            length(which(edge_summary$window == 2))),

  # time window
  window = edge_summary$window,                              # window ID

  # # edge means per time window
  # logit_edge_mu_1 = logit_edge_draws_mu[[1]],                # Mean logit edge weights, window 1
  # logit_edge_mu_2 = logit_edge_draws_mu[[2]],                # Mean logit edge weights, window 2
  # logit_edge_mu_3 = logit_edge_draws_mu[[3]],                # Mean logit edge weights, window 3
  #
  # # covariance matrices per time window
  # logit_edge_cov_1 = logit_edge_draws_cov[[1]],              # SD logit edge weights, window 1
  # logit_edge_cov_2 = logit_edge_draws_cov[[2]],              # SD logit edge weights, window 2
  # logit_edge_cov_3 = logit_edge_draws_cov[[3]],              # SD logit edge weights, window 3

  # Gaussian approximation of edge weights
  logit_edge_mu = logit_edge_draws_mu,                         # Means of Gaussian approximation of logit edge weights
  logit_edge_sd = logit_edge_draws_sd,                         # SD of Gaussian approximation of logit edge weights

  # explanatory variables
  age_min = edge_summary$age_min_std,                          # age of younger dyad member
  age_max = edge_summary$age_max_std,                          # age of older dyad member

  # multimembership terms
  node_1 = edge_summary$node_1_random,                         # Node 1 IDs for multimembership terms
  node_2 = edge_summary$node_2_random,                         # Node 2 IDs for multimembership terms

  # # dyad IDs for each time windows
  # dyads_window1 = edge_summary$dyad_random[edge_summary$window == 1],  # Dyad IDs, window 1
  # dyads_window2 = edge_summary$dyad_random[edge_summary$window == 2],  # Dyad IDs, window 2
  # dyads_window3 = edge_summary$dyad_random[edge_summary$window == 3],  # Dyad IDs, window 3

  # dyad IDs for all time windows
  dyad_id = edge_summary$dyad_random
)

## print progress marker
print('data list created')

## load dyadic regression model
dyadic_regression <- stan_model('models/dyadic_regression_anp_sd.stan')
#dyadic_regression <- cmdstan_model('models/dyadic_regression_anp_sd.stan')

## print progress marker
print('model loaded in')

## define parameters
n_chains <- 4
n_samples <- 1000

## fit dyadic regression
print('start model')
fit_dyadreg_sim <- sampling(dyadic_regression,
                            data = dyad_data,
                            iter = n_samples*2, warmup = n_samples,
                            chains = n_chains, cores = n_chains
                            )
# fit_dyadreg_sim <- dyadic_regression$sample(
#   data = dyad_data,
#   iter_warmup = n_samples,
#   iter_sampling = n_samples,
#   chains = n_chains,
#   parallel_chains = n_chains)

save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

## print progress marker
print('model run')

#### check outputs ####
load('step5_dyadicregression/anp_dyadic_simulation.RData')

## extract model fit
# summary <- fit_dyadreg_sim$summary()
summary <- rstan::summary(fit_dyadreg_sim)
(summary <- as.data.frame(summary$summary))
par(mfrow = c(3,1))
# hist(summary$rhat, breaks = 50)
# hist(summary$ess_bulk, breaks = 50)
# hist(summary$ess_tail, breaks = 50)
hist(summary$Rhat, breaks = 50)
hist(summary$n_eff, breaks = 50)
par(mfrow = c(1,1))

## extract draws
# draws <- fit_dyadreg_sim$draws(format = 'df')
draws <- rstan::extract(fit_dyadreg_sim)

## extract dyadic regression parameters
b_max <- draws$beta_age_max
b_min <- draws$beta_age_min
intercept <- draws$intercept
sigma <- draws$sigma
parameters <- data.frame(beta_age_max = b_max,
                         beta_age_min = b_min,
                         intercept = intercept,
                         sigma = sigma) %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = c('beta_age_max','beta_age_min','sigma','intercept'),
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
    mutate(chain = rep(1:4, each = 1000),
           position = rep(1:1000, 4)) %>%
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

print('multimembership parameters plotted')

## extract window random effect draws
rand_window <- as.data.frame(draws$rand_window)
colnames(rand_window) <- unique(dyad_data$window)
rand_window_long <- long_random(rand_window, n_windows)

## window random effect traceplots
traceplot(rand_window_long)

print('window parameters plotted')

## extract dyad random effect draws
rand_dyad <- as.data.frame(draws$rand_dyad)
colnames(rand_dyad) <- unique(dyad_data$dyad)
rand_dyad_long <- long_random(rand_dyad, n_dyads)

## dyad random effect traceplots
rand_dyad_long %>%
  filter(parameter %in% sample(unique(rand_dyad_long$parameter), 25, replace = F)) %>%
  traceplot()

## print progress marker
print('outputs checked')

#### posterior predictive check ####
plot(density(as.numeric(sim_dat[1, ])),
     main = "Posterior predictive density of edge weights:\nblack = measured edge, red = predicted edge",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:nrow(sim_dat), 1)

  mu_plot <- rep(NA, n_dyads)
  for(k in 1:n_dyads){
    mu_plot[k] <- intercept[j] + b_min[j]*dyad_data$age_min[k] + b_max[j]*dyad_data$age_max[k] + mm_nodes[j,dyad_data$node_1[k]] + mm_nodes[j,dyad_data$node_2[k]]
  }

  sigma_plot <- dyad_data$logit_edge_sd + diag(rep(sigma[j], n_data))
  norm <- rnorm(1000, mu_plot, sigma_plot)

  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions

}

save.image('step5_dyadicregression/anp_dyadic_simulation.RData')
print('posterior predictive check complete')

#### predict from raw data ####
load('step5_dyadicregression/anp_dyadic_simulation.RData')

# ## create predictive data frame
# age_min_length <- length(unique(dyad_data$age_min))
# age_max_length <- length(unique(dyad_data$age_max))
# pred <- data.frame(age_min_std = rep(rep(unique(dyad_data$age_min),
#                                      each = age_max_length),
#                                  n_chains*n_samples),
#                    age_max_std = rep(rep(unique(dyad_data$age_max),
#                                      age_min_length),
#                                  n_chains*n_samples),
#                    age_min = rep(rep(unique(edge_summary$age_min),
#                                          each = age_max_length),
#                                      n_chains*n_samples),
#                    age_max = rep(rep(unique(edge_summary$age_max),
#                                          age_min_length),
#                                      n_chains*n_samples),
#                    intcp = rep(intercept,
#                                each = age_max_length*age_min_length),
#                    beta_min = rep(b_min,
#                                   each = age_max_length*age_min_length),
#                    beta_max = rep(b_max,
#                                   each = age_max_length*age_min_length)) %>%
#   filter(age_min_std <= age_max_std)
#
# ## calculate predictions
# pred$pred_mu <- pred$intcp + pred$beta_min * pred$age_min_std + pred$beta_max * pred$age_max_std

## create function for predicting mean
get_mean_predictions <- function(prediction_df, draws_list){
  mu <- matrix(NA, nrow = length(draws_list$intercept), ncol = nrow(prediction_df),
               dimnames = list(1:length(draws_list$intercept),
                               paste0('dyad',prediction_df$dyad_random)))

  multimembership <- draws_list$mm_nodes
  window_random <- draws_list$rand_window
  dyad_random <- draws_list$rand_dyad
  intcp <- draws_list$intercept
  b_max <- draws_list$beta_age_max
  b_min <- draws_list$beta_age_min

  for(dyad in 1:ncol(mu)){
    node_1 <- multimembership[,prediction_df$node_1_random[dyad]]
    node_2 <- multimembership[,prediction_df$node_2_random[dyad]]
    window_effect <- window_random[,prediction_df$window[dyad]]
    dyad_effect <- dyad_random[,prediction_df$dyad_random[dyad]]

    for(draw in 1:nrow(mu)){
      mu[draw, dyad] <- intcp[draw] +
        b_max[draw] * prediction_df$age_max_std[dyad] +
        b_min[draw] * prediction_df$age_min_std[dyad] +
        node_1[draw] + node_2[draw] +
        window_effect[draw] + dyad_effect[draw]
    }
  }
  return(mu)
}

## predict mean distribution
pred_mean <- get_mean_predictions(prediction_df = edge_summary,
                                  draws_list = draws)

## summarise mean distribution
edge_summary$pred_mu <- apply(pred_mean, 2, mean)
edge_summary$pred_mean_lwr <- apply(pred_mean, 2, quantile, prob = 0.025)
edge_summary$pred_mean_upr <- apply(pred_mean, 2, quantile, prob = 0.975)

## create function for predicting full distribution
get_full_predictions <- function(mu_matrix, sigma){
  predictions <- mu_matrix
  for(dyad in 1:ncol(predictions)){
    for(draw in 1:nrow(predictions)){
      predictions[draw, dyad] <- rnorm(n = 1,
                                       mean = mu_matrix[draw, dyad],
                                       sd = sigma[draw])
    }
  }
  return(predictions)
}

## predict full distribution
pred_full <- get_full_predictions(mu_matrix = pred_mean,
                                  sigma = draws$sigma)

## summarise full distribution
edge_summary$pred_full_lwr <- apply(pred_full, 2, quantile, prob = 0.025)
edge_summary$pred_full_upr <- apply(pred_full, 2, quantile, prob = 0.975)

## convert to unstandardised scale
pred_mean_ustd <- pred_mean * sd(edge_summary$mu) + mean(edge_summary$mu)
pred_full_ustd <- pred_full * sd(edge_summary$mu) + mean(edge_summary$mu)

## summarise unstandardised distributions
edge_summary$pred_mu_ustd <- apply(pred_mean_ustd, 2, mean)
edge_summary$pred_mean_lwr_ustd <- apply(pred_mean_ustd, 2, quantile, prob = 0.025)
edge_summary$pred_mean_upr_ustd <- apply(pred_mean_ustd, 2, quantile, prob = 0.975)
edge_summary$pred_full_lwr_ustd <- apply(pred_full_ustd, 2, quantile, prob = 0.025)
edge_summary$pred_full_upr_ustd <- apply(pred_full_ustd, 2, quantile, prob = 0.975)

## convert to invlogit scale
pred_mean_invlogit <- invlogit(pred_mean_ustd)
pred_full_invlogit <- invlogit(pred_full_ustd)

## summarise invlogit distributions
edge_summary$pred_mu_invlogit <- apply(pred_mean_invlogit, 2, mean)
edge_summary$pred_mean_lwr_invlogit <- apply(pred_mean_invlogit, 2, quantile, prob = 0.025)
edge_summary$pred_mean_upr_invlogit <- apply(pred_mean_invlogit, 2, quantile, prob = 0.975)
edge_summary$pred_full_lwr_invlogit <- apply(pred_full_invlogit, 2, quantile, prob = 0.025)
edge_summary$pred_full_upr_invlogit <- apply(pred_full_invlogit, 2, quantile, prob = 0.975)

# ## summarise
# pred_summary <- pred %>%
#   group_by(age_min_std, age_max_std) %>%
#   mutate(pred_lwr = quantile(pred_mu, probs = 0.025),
#          pred_mean = mean(pred_mu),
#          pred_upr = quantile(pred_mu, probs = 0.975),
#          pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
#          pred_unstd_mean = mean(pred_unstd),
#          pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>%
#   ungroup() %>%
#   select(age_min_std, age_max_std,
#          pred_lwr, pred_mean, pred_upr,
#          pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>%
#   distinct()

#### plot predictions ####
# pred_summary <- pred_summary %>%
#   left_join(distinct(edge_summary[,c('age_min','age_min_std')]),
#             by = 'age_min_std') %>%
#   left_join(distinct(edge_summary[,c('age_max','age_max_std')]),
#             by = 'age_max_std')
ggplot()+
  geom_ribbon(data = edge_summary,
              aes(x = age_min,
                  ymin = pred_full_lwr, ymax = pred_full_upr,
                  group = as.factor(age_max), fill = as.factor(age_max)),
              alpha = 0.3)+
  geom_line(data = edge_summary,
            aes(x = age_min,
                y = pred_mu,
                colour = as.factor(age_max), group = as.factor(age_max)),
            linewidth = 1)+
  scale_colour_viridis_d(direction = -1)+ scale_fill_viridis_d(direction = -1)+
  geom_point(data = edge_summary,
             aes(x = age_min,
                 y = mu_std,
                 colour = as.factor(age_max))
             )+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        legend.position = 'bottom', #c(0.8,0.2),
        # legend.title = element_text(size = 14),
        # legend.text = element_text(size = 12)
        axis.title = element_text(size = 18))+
  labs(colour = 'maximum age', fill = 'maximum age')+
  guides(colour = guide_legend(ncol = 6),
         fill = guide_legend(ncol = 6))

# raw_summary <- edge_summary %>%
#   select(age_min, age_max, mu, mu_std) %>%
#   group_by(age_min, age_max) %>%
#   mutate(mu_mean = mean(mu),
#          mu_std_mean = mean(mu_std)) %>%
#   ungroup() %>%
#   select(-mu, -mu_std) %>%
#   distinct() %>%
#   mutate(age_pair = paste0(age_min, '_', age_max))
# compare <- pred_summary %>%
#   mutate(age_pair = paste0(age_min, '_', age_max)) %>%
#   left_join(raw_summary[,c('mu_mean','mu_std_mean','age_pair')], by = 'age_pair') %>%
#   rename(raw_unstd = mu_mean, raw_std = mu_std_mean) %>%
#   select(age_min, age_max, raw_unstd, pred_unstd_mean, raw_std, pred_mean, pred_lwr, pred_unstd_lwr, pred_upr, pred_unstd_upr)
# ggplot(compare)+
#   geom_ribbon(aes(x = raw_unstd,
#                   ymin = pred_unstd_lwr, ymax = pred_unstd_upr),
#               alpha = 0.3)+
#   geom_line(aes(x = raw_unstd, y = pred_unstd_mean),
#             linewidth = 1)+
#   geom_line(data = data.frame(x = c(min(edge_summary$mu),max(edge_summary$mu)),
#                               y = c(min(edge_summary$mu),max(edge_summary$mu))),
#             aes(x = x, y = y),
#             linewidth = 0.5, colour = 'red')+
#   scale_x_continuous('raw mean')+
#   scale_y_continuous('predicted mean')+
#   theme(axis.text = element_text(size = 14),
#         axis.title = element_text(size = 18),
#         legend.title = element_text(size = 14),
#         legend.text = element_text(size = 12))

ggplot()+
  geom_point(data = edge_summary,
             mapping = aes(x = mu_std, y = pred_mu))+
  geom_abline(slope = 1, intercept = 0)

ggplot()+
  geom_point(data = edge_summary,
             mapping = aes(x = mu, y = pred_mu_ustd))+
  geom_abline(slope = 1, intercept = 0)

ggplot()+
  geom_point(data = edge_summary,
             mapping = aes(x = invlogit(mu), y = pred_mu_invlogit))+
  geom_abline(slope = 1, intercept = 0)

dev.off()
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

## print progress marker
print('predictions complete')

#### extract original parameters ####
load('step5_dyadicregression/anp_dyadic_simulation.RData')

## create altered dataframes to predict from
pred_new <- edge_summary %>%
  rename(age_min_org = age_min,
         age_max_org = age_max,
         age_min_org_std = age_min_std,
         age_max_org_std = age_max_std) %>%
  mutate(age_min_new = age_min_org,
         age_max_new = age_max_org)
# pred_new_min <- pred_new %>%
#   mutate(age_min_new = age_min_new + 1) %>%
#   mutate(age_min_new_std = (age_min_new - mean_age)/stdv_age,
#          age_max_new_std = (age_max_new - mean_age)/stdv_age) %>%
#   filter(age_min_new <= age_max_new)  # filter out any where min is now bigger than max
# pred_new_max <- pred_new %>%
#   mutate(age_max_new = age_max_new + 1) %>%
#   mutate(age_min_new_std = (age_min_new - mean_age)/stdv_age,
#          age_max_new_std = (age_max_new - mean_age)/stdv_age)
# pred_new_both <- pred_new %>%
#   mutate(age_min_new = age_min_new + 1,
#          age_max_new = age_max_new + 1) %>%
#   mutate(age_min_new_std = (age_min_new - mean_age)/stdv_age,
#          age_max_new_std = (age_max_new - mean_age)/stdv_age)
pred_new_min <- pred_new %>%
  mutate(age_min_new_std = ((age_min_new - mean_age)/stdv_age) + 1,
         age_max_new_std = (age_max_new - mean_age)/stdv_age) %>%
  mutate(age_min_new = age_min_new_std * stdv_age + mean_age,
         age_max_new = age_max_new_std * stdv_age + mean_age) %>%
  filter(age_min_new <= age_max_new)  # filter out any where min is now bigger than max

pred_new_max <- pred_new %>%
  mutate(age_min_new_std = ((age_min_new - mean_age)/stdv_age),
         age_max_new_std = ((age_max_new - mean_age)/stdv_age) + 1) %>%
  mutate(age_min_new = age_min_new_std * stdv_age + mean_age,
         age_max_new = age_max_new_std * stdv_age + mean_age)

pred_new_both <- pred_new %>%
  mutate(age_min_new_std = ((age_min_new - mean_age)/stdv_age) + 1,
         age_max_new_std = ((age_max_new - mean_age)/stdv_age) + 1) %>%
  mutate(age_min_new = age_min_new_std * stdv_age + mean_age,
         age_max_new = age_max_new_std * stdv_age + mean_age)

rm(pred_new) ; gc()

## adapt function for predicting mean
get_mean_predictions <- function(prediction_df, draws_list, new){
  mu <- matrix(NA, nrow = length(draws_list$intercept), ncol = nrow(prediction_df),
               dimnames = list(1:length(draws_list$intercept),
                               paste0('dyad',prediction_df$dyad_random)))

  multimembership <- draws_list$mm_nodes
  window_random <- draws_list$rand_window
  dyad_random <- draws_list$rand_dyad
  intcp <- draws_list$intercept
  b_max <- draws_list$beta_age_max
  b_min <- draws_list$beta_age_min

  if(is.null(new) == FALSE){
    if(new == TRUE){
      prediction_df$age_min <- prediction_df$age_min_new
      prediction_df$age_max <- prediction_df$age_max_new
      prediction_df$age_min_std <- prediction_df$age_min_new_std
      prediction_df$age_max_std <- prediction_df$age_max_new_std
    }
    if(new == FALSE){
      prediction_df$age_min <- prediction_df$age_min_org
      prediction_df$age_max <- prediction_df$age_max_org
      prediction_df$age_min_std <- prediction_df$age_min_org_std
      prediction_df$age_max_std <- prediction_df$age_max_org_std
    }

  }

  for(dyad in 1:ncol(mu)){
    node_1 <- multimembership[,prediction_df$node_1_random[dyad]]
    node_2 <- multimembership[,prediction_df$node_2_random[dyad]]
    window_effect <- window_random[,prediction_df$window[dyad]]
    dyad_effect <- dyad_random[,prediction_df$dyad_random[dyad]]

    for(draw in 1:nrow(mu)){
      mu[draw, dyad] <- intcp[draw] +
        b_max[draw] * prediction_df$age_max_std[dyad] +
        b_min[draw] * prediction_df$age_min_std[dyad] +
        node_1[draw] + node_2[draw] +
        window_effect[draw] + dyad_effect[draw]
    }
  }
  return(mu)
}

## predict mean distributions
pred_org_min_mean <- get_mean_predictions(prediction_df = pred_new_min,
                                          draws_list = draws,
                                          new = FALSE)
pred_new_min_mean <- get_mean_predictions(prediction_df = pred_new_min,
                                          draws_list = draws,
                                          new = TRUE)
pred_new_max_mean <- get_mean_predictions(prediction_df = pred_new_max,
                                          draws_list = draws,
                                          new = TRUE)
pred_new_both_mean <- get_mean_predictions(prediction_df = pred_new_both,
                                           draws_list = draws,
                                           new = TRUE)

## predict full distribution
pred_org_min_full <- get_full_predictions(mu_matrix = pred_org_min_mean,
                                          sigma = draws$sigma)
pred_new_min_full <- get_full_predictions(mu_matrix = pred_new_min_mean,
                                          sigma = draws$sigma)
pred_new_max_full <- get_full_predictions(mu_matrix = pred_new_max_mean,
                                          sigma = draws$sigma)
pred_new_both_full <- get_full_predictions(mu_matrix = pred_new_both_mean,
                                           sigma = draws$sigma)

## convert to unstandardised scale
pred_org_min_full_ustd <- pred_org_min_full * sd(edge_summary$mu) + mean(edge_summary$mu)
pred_new_min_full_ustd <- pred_new_min_full * sd(edge_summary$mu) + mean(edge_summary$mu)
pred_new_max_full_ustd <- pred_new_max_full * sd(edge_summary$mu) + mean(edge_summary$mu)
pred_new_both_full_ustd <- pred_new_both_full * sd(edge_summary$mu) + mean(edge_summary$mu)

## calculate contrasts: minimum age increases by 1, standardised scale
contrast_min <- pred_new_min_full - pred_org_min_full
print(paste0('contrast for age_min + 1 (std): ',
             round(mean(contrast_min), 5), ' ± ', round(sd(contrast_min), 5),
             ' [', round(quantile(contrast_min, prob = 0.025), 5),
             ' - ',round(quantile(contrast_min, prob = 0.975), 5), ']'))

## calculate contrasts: minimum age increases by 1, outcome scale
contrast_min_ustd <- pred_new_min_full_ustd - pred_org_min_full_ustd
print(paste0('contrast for age_min + 1 (ustd): ',
             round(mean(contrast_min_ustd), 5), ' ± ', round(sd(contrast_min_ustd), 5),
             ' [', round(quantile(contrast_min_ustd, prob = 0.025), 5),
             ' - ',round(quantile(contrast_min_ustd, prob = 0.975), 5), ']'))
print(paste0('true minimum effect: ', sim_min))

## calculate contrasts: maximum age increases by 1, standardised scale
contrast_max <- pred_new_max_full - pred_full
print(paste0('contrast for age_max + 1 (std): ',
             round(mean(contrast_max), 5), ' ± ', round(sd(contrast_max), 5),
             ' [', round(quantile(contrast_max, prob = 0.025), 5),
             ' - ',round(quantile(contrast_max, prob = 0.975), 5), ']'))

## calculate contrasts: maximum age increases by 1, outcome scale
contrast_max_ustd <- pred_new_max_full_ustd - pred_full_ustd
print(paste0('contrast for age_max + 1 (ustd): ',
             round(mean(contrast_max_ustd), 5), ' ± ', round(sd(contrast_max_ustd), 5),
             ' [', round(quantile(contrast_max_ustd, prob = 0.025), 5),
             ' - ',round(quantile(contrast_max_ustd, prob = 0.975), 5), ']'))
print(paste0('true minimum effect: ', sim_max))

## calculate contrasts: both ages increase by 1, standardised scale
contrast_both <- pred_new_both_full - pred_full
print(paste0('contrast for both ages + 1 (std): ',
             round(mean(contrast_both), 5), ' ± ', round(sd(contrast_both), 5),
             ' [', round(quantile(contrast_both, prob = 0.025), 5),
             ' - ',round(quantile(contrast_both, prob = 0.975), 5), ']'))

## calculate contrasts: both ages increase by 1, outcome scale
contrast_both_ustd <- pred_new_both_full_ustd - pred_full_ustd
print(paste0('contrast for both ages + 1 (ustd): ',
             round(mean(contrast_both_ustd), 5), ' ± ', round(sd(contrast_both_ustd), 5),
             ' [', round(quantile(contrast_both_ustd, prob = 0.025), 5),
             ' - ',round(quantile(contrast_both_ustd, prob = 0.975), 5), ']'))

## save workspace
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

#### final clean plots ####
load('step5_dyadicregression/anp_dyadic_simulation.RData')

## create counterfactual data frame for plotting from: minimum age on x
fake_data <- expand.grid(window = 1:n_windows,
                         dyad_random = unique(edge_summary$dyad_random),
                         age_min = min(node_data$age):max(node_data$age),
                         age_max = c(15, 20, 25, 30, 35, 40, 45, 50)) %>% 
  left_join(distinct(edge_summary[,c('dyad_random','node_1','node_2')]),
            by = 'dyad_random') %>% 
  filter(age_max >= age_min) %>% 
  mutate(age_min_std = (age_min - mean_age) / stdv_age,
         age_max_std = (age_max - mean_age) / stdv_age)

## predict mean values
fake_mean <- get_mean_predictions(prediction_df = fake_data, draws_list = draws, new = NULL)
warnings()
length(which(is.na(fake_mean) == TRUE))
print('mean values extracted')

## predict full distribution
fake_full <- get_full_predictions(mu_matrix = fake_mean, sigma = draws$sigma)
print('full values extracted')
length(which(is.na(fake_full) == TRUE))

## summarise values
fake_data$mean_prediction <- apply(fake_mean, 2, mean, na.rm = T)
fake_data$mu_lwr_prediction <- apply(fake_mean, 2, quantile, prob = 0.025, na.rm = T)
fake_data$mu_upr_prediction <- apply(fake_mean, 2, quantile, prob = 0.975, na.rm = T)
fake_data$full_lwr_prediction <- apply(fake_full, 2, quantile, prob = 0.025, na.rm = T)
fake_data$full_upr_prediction <- apply(fake_full, 2, quantile, prob = 0.975, na.rm = T)
print('predictions summarised')
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

## summarise summary
summary_fake <- fake_data %>% 
  group_by(age_min, age_max) %>% 
  summarise(mean = mean(mean_prediction),
            mu_lwr = mean(mu_lwr_prediction),
            mu_upr = mean(mu_upr_prediction),
            full_lwr = mean(full_lwr_prediction),
            full_upr = mean(full_upr_prediction))
print('summary summarised')

## plot altogether
(min_overall <- ggplot()+
    geom_ribbon(data = summary_fake,
                mapping = aes(x = age_min,
                              ymin = full_lwr, ymax = full_upr,
                              fill = age_max),
                alpha = 0.4)+
    geom_ribbon(data = summary_fake,
                mapping = aes(x = age_min,
                              ymin = mu_lwr, ymax = mu_upr,
                              fill = age_max),
                alpha = 0.4)+
    geom_line(data = summary_fake,
              mapping = aes(x = age_min,
                            y = mean,
                            colour = age_max))+
    geom_point(data = edge_summary,
               mapping = aes(x = age_min,
                             y = mu_std,
                             colour = age_max))+
    scale_colour_viridis_c(begin = (1/8), end = 1)+
    scale_fill_viridis_c(begin = (1/8), end = 1)+
    theme(legend.position = 'bottom')+
    labs(fill = 'age of older male',
         colour = 'age of older male'))
ggsave(plot = min_overall, device = 'png',
       filename = 'anp_dyadic_simulation_minageeffect_overall.png',
       path = 'step5_dyadicregression/',
       width = 1200, height = 800, units = 'px')
print('overall plot completed with minimum age')

## resummarise summary
summary_fake <- fake_data %>% 
  group_by(age_min, age_max, window) %>% 
  summarise(mean = mean(mean_prediction),
            mu_lwr = mean(mu_lwr_prediction),
            mu_upr = mean(mu_upr_prediction),
            full_lwr = mean(full_lwr_prediction),
            full_upr = mean(full_upr_prediction))

## plot separately by window
(min_bywindow <- ggplot()+
    geom_ribbon(data = summary_fake,
                mapping = aes(x = age_min,
                              ymin = full_lwr, ymax = full_upr,
                              fill = age_max),
                alpha = 0.4)+
    geom_ribbon(data = summary_fake,
                mapping = aes(x = age_min,
                              ymin = mu_lwr, ymax = mu_upr,
                              fill = age_max),
                alpha = 0.4)+
    geom_line(data = summary_fake,
              mapping = aes(x = age_min,
                            y = mean,
                            colour = age_max))+
    geom_point(data = edge_summary,
               mapping = aes(x = age_min,
                             y = mu_std,
                             colour = age_max))+
    scale_colour_viridis_c(begin = (1/8), end = 1)+
    scale_fill_viridis_c(begin = (1/8), end = 1)+
    theme(legend.position = 'bottom')+
    labs(fill = 'age of older male',
         colour = 'age of older male')+
    facet_wrap(. ~ window))
ggsave(plot = min_bywindow, device = 'png',
       filename = 'anp_dyadic_simulation_minageeffect_bywindow.png',
       path = 'step5_dyadicregression/',
       width = 1600, height = 2400, units = 'px')
print('faceted plot completed with minimum age')
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

## create counterfactual data frame for plotting from: maximum age on x
fake_data <- expand.grid(window = 1:n_windows,
                         dyad_random = unique(edge_summary$dyad_random),
                         age_min = c(10, 15, 20, 25, 30, 35, 40, 45),
                         age_max = min(node_data$age):max(node_data$age)) %>% 
  left_join(distinct(edge_summary[,c('dyad_random','node_1','node_2')]),
            by = 'dyad_random') %>% 
  filter(age_max >= age_min) %>% 
  mutate(age_min_std = (age_min - mean_age) / stdv_age,
         age_max_std = (age_max - mean_age) / stdv_age)

## predict mean values
fake_mean <- get_mean_predictions(prediction_df = fake_data, draws_list = draws)
print('mean values extracted')

## predict full distribution
fake_full <- get_full_predictions(mu_matrix = fake_mean, sigma = draws$sigma)
print('full values extracted')

## summarise values
fake_data$mean_prediction <- apply(fake_mean, 2, mean)
fake_data$mu_lwr_prediction <- apply(fake_mean, 2, quantile, prob = 0.025)
fake_data$mu_upr_prediction <- apply(fake_mean, 2, quantile, prob = 0.975)
fake_data$full_lwr_prediction <- apply(fake_full, 2, quantile, prob = 0.025)
fake_data$full_upr_prediction <- apply(fake_full, 2, quantile, prob = 0.975)
print('predictions summarised')
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

## summarise summary
summary_fake <- fake_data %>% 
  group_by(age_min, age_max) %>% 
  summarise(mean = mean(mean_prediction),
            mu_lwr = mean(mu_lwr_prediction),
            mu_upr = mean(mu_upr_prediction),
            full_lwr = mean(full_lwr_prediction),
            full_upr = mean(full_upr_prediction))
print('summary summarised')

## plot altogether
(max_overall <- ggplot()+
    geom_ribbon(data = summary_fake,
                mapping = aes(x = age_max,
                              ymin = full_lwr, ymax = full_upr,
                              fill = age_min),
                alpha = 0.4)+
    geom_ribbon(data = summary_fake,
                mapping = aes(x = age_max,
                              ymin = mu_lwr, ymax = mu_upr,
                              fill = age_min),
                alpha = 0.4)+
    geom_line(data = summary_fake,
              mapping = aes(x = age_max,
                            y = mean,
                            colour = age_min))+
    geom_point(data = edge_summary,
               mapping = aes(x = age_max,
                             y = mu_std,
                             colour = age_min))+
    scale_colour_viridis_c(begin = 0, end = (7/8))+
    scale_fill_viridis_c(begin = 0, end = (7/8))+
    theme(legend.position = 'bottom')+
    labs(fill = 'age of younger male',
         colour = 'age of younger male'))
ggsave(plot = max_overall, device = 'png',
       filename = 'anp_dyadic_simulation_maxageeffect_overall.png',
       path = 'step5_dyadicregression/',
       width = 1200, height = 800, units = 'px')
print('overall plot completed with minimum age')

## resummarise summary
summary_fake <- fake_data %>% 
  group_by(age_min, age_max, window) %>% 
  summarise(mean = mean(mean_prediction),
            mu_lwr = mean(mu_lwr_prediction),
            mu_upr = mean(mu_upr_prediction),
            full_lwr = mean(full_lwr_prediction),
            full_upr = mean(full_upr_prediction))

## plot separately by window
(max_bywindow <- ggplot()+
    geom_ribbon(data = summary_fake,
                mapping = aes(x = age_max,
                              ymin = full_lwr, ymax = full_upr,
                              fill = age_min),
                alpha = 0.4)+
    geom_ribbon(data = summary_fake,
                mapping = aes(x = age_max,
                              ymin = mu_lwr, ymax = mu_upr,
                              fill = age_min),
                alpha = 0.4)+
    geom_line(data = summary_fake,
              mapping = aes(x = age_max,
                            y = mean,
                            colour = age_min))+
    geom_point(data = edge_summary,
               mapping = aes(x = age_max,
                             y = mu_std,
                             colour = age_min))+
    scale_colour_viridis_c(begin = 0, end = (7/8))+
    scale_fill_viridis_c(begin = 0, end = (7/8))+
    theme(legend.position = 'bottom')+
    labs(fill = 'age of younger male',
         colour = 'age of younger male')+
    facet_wrap(. ~ window))
ggsave(plot = max_bywindow, device = 'png',
       filename = 'anp_dyadic_simulation_maxageeffect_bywindow.png',
       path = 'step5_dyadicregression/',
       width = 1600, height = 2400, units = 'px')
print('faceted plot completed with minimum age')

## save output
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')
dev.off()
