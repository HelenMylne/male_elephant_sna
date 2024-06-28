#### set up ####
# script to simulate dyadic regression and check actually working
# library(StanHeaders) ; library(rstan) ; library(tidyverse) ; library(car) ; library(LaplacesDemon)
library(StanHeaders, lib.loc = '../packages/')    # library(StanHeaders)
library(rstan, lib.loc = '../packages/')          # library(rstan)
library(tidyverse, lib.loc = '../packages/')      # library(tidyverse)
library(car, lib.loc = '../packages/')            # library(car)
library(LaplacesDemon, lib.loc = '../packages/')  # library(LaplacesDemon)

theme_set(theme_classic())

pdf('step5_dyadicregression/simulate_dyadic_anp.pdf')

#### load raw data as basis for simulation -- simulating using first 3 short windows only for sake of speed and computing power ####
## load time window 1 and remove additional data
load('anp_edgecalculations/anpshort1_edgeweights_conditionalprior.RData')
rm(list = ls()[! ls() %in% c('nodes','cdf_1','edges','edge_samples')]) ; gc()

## create full data frame with window variable
nodes_all <- nodes %>%
  mutate(window = 1)
cdf_all <- cdf_1 %>%
  dplyr::select(-dyad_rank) %>%
  mutate(window = 1)
edges_all <- list()
edges_all[[1]] <- edge_samples

## attach all other data frames
rm(list = ls()[! ls() %in% c('nodes_all','cdf_all','edges_all')]) ; gc()
n_windows <- 3
for(time_window in 2:n_windows){
  ## import workspace image for time window
  load(paste0('anp_edgecalculations/anpshort',time_window,'_edgeweights_conditionalprior.RData'))

  ## ensure all have the same names: edge weights
  if('edge_weights_matrix' %in% ls()){
    edge_samples <- edge_weights_matrix
  } else {
    if('edge_weight_matrix' %in% ls()){
      edge_samples <- edge_weight_matrix
    }
  }

  ## ensure all have the same names: dyads data frame
  if('counts_df' %in% ls()){
    cdf <- counts_df
  }

  ## add window variable
  nodes <- nodes %>%
    mutate(window = time_window) %>%
    dplyr::select(colnames(nodes_all))
  cdf <- cdf %>%
    mutate(window = time_window) %>%
    dplyr::select(colnames(cdf_all))

  ## attach to full data frame
  nodes_all <- rbind(nodes_all, nodes)
  cdf_all <- rbind(cdf_all, cdf)
  edges_all[[time_window]] <- edges_all

  ## print progress marker
  print(time_window)

  ## clean environment
  rm(list = ls()[! ls() %in% c('nodes_all','cdf_all','edges_all','time_window')]) ; gc()
}

## create data frame in which dyad_id is randomised
dyads_random <- cdf_all %>%
  dplyr::select(dyad_id, node_1, node_2) %>%
  distinct()
dyads_random$dyad_random <- sample(1:nrow(dyads_random), replace = FALSE)

## join to model data frame
cdf_all <- cdf_all %>%
  left_join(dyads_random, by = c('node_1','node_2','dyad_id')) %>%
  relocate(dyad_random, .after = dyad_id)

## combine with window ID
cdf_all$dyad_window <- paste0(cdf_all$dyad_random, '_', cdf_all$window)

## plot dyad ID vs window -- no pattern of dyad ID with window
ggplot(cdf_all)+
  geom_point(aes(x = dyad_random, y = window))

## clean up
rm(dyads_random) ; gc()
save.image('step5_dyadicregression/anp_simulation.RData')
print('dyadic data created')

#### simulate age effect ####
## cut down to only 3 windows for sake of time and effort
cdf_all <- cdf_all %>% filter(window < 4)
nodes_all <- nodes_all %>% filter(window < 4)
edges_all <- edges_all[1:3]

## set "true" values
sim_min <- -0.2
sim_max <- 0.8
sim_int <- 5

## identify global parameters
n_data <- nrow(cdf_all)
n_dyads <- length(unique(cdf_all$dyad_random))
n_nodes <- length(unique(nodes_all$id))

## identify older and younger dyad members
cdf_all$age_min <- ifelse(cdf_all$age_start_1 < cdf_all$age_start_2,
                          cdf_all$age_start_1, cdf_all$age_start_2)
cdf_all$age_max <- ifelse(cdf_all$age_start_1 > cdf_all$age_start_2,
                          cdf_all$age_start_1, cdf_all$age_start_2)

## set mean edge weight
cdf_all$mu <- sim_int + cdf_all$age_min * sim_min + cdf_all$age_max * sim_max   # simulate mean centrality on normal scale
# ggplot(cdf_all)+
#   geom_point(aes(x = age_max, y = mu, colour = age_min))

## standardise edges
cdf_all$mu_std <- ( cdf_all$mu - mean(cdf_all$mu) ) / sd(cdf_all$mu)
# ggplot(cdf_all)+
#   geom_point(aes(x = age_max, y = mu_std, colour = age_min))

## simulate full distribution of samples per node
cdf_all$sd <- abs(sim_max/3)           # make small to start with to be sure model should be able to detect difference
sim_dat <- matrix(data = NA, nrow = 1000, ncol = n_data,
                  dimnames = list(1:1000, cdf_all$dyad_id))    # create matrix
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = cdf_all$mu_std[j], sd = cdf_all$sd[j])  # simulate distribution
}
plot(sim_dat[1,] ~ cdf_all$age_min)          # plot simulated values against minimum age
plot(sim_dat[1,] ~ cdf_all$age_max)          # plot simulated values against maximum age

## print progress marker
print('edges simulated')
save.image('step5_dyadicregression/anp_simulation.RData')

#### plot raw data ####

cdf_all$age_cat_max <- NA ; cdf_all$age_cat_min <- NA
for(i in 1:nrow(cdf_all)){
  cdf_all$age_cat_max[i] <- max(cdf_all$age_num_1[i], cdf_all$age_num_2[i])
  cdf_all$age_cat_min[i] <- min(cdf_all$age_num_1[i], cdf_all$age_num_2[i])
}

# ggplot()+
#   geom_point(data = cdf_all, aes(x = age_min,
#                              y = mu_std,
#                              colour = as.factor(age_cat_max)),
#              shape = 19)+
#   scale_x_continuous('age of younger dyad member')+
#   scale_y_continuous('mean estimated edge weight')
# 
# ggplot()+
#   geom_boxplot(data = cdf_all, aes(x = as.factor(age_cat_min),
#                                y = mu_std,
#                                fill = as.factor(age_cat_max)),
#                shape = 19)+
#   scale_x_discrete('age category of younger dyad member')+
#   scale_y_continuous('mean estimated edge weight')+
#   theme(legend.position = 'bottom')+
#   labs(fill = 'age category of older male')
# 
# sim_dat_long <- sim_dat %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = everything(), names_to = 'dyad_id', values_to = 'edge_draw') %>% 
#   mutate(dyad_id = as.numeric(dyad_id)) %>% 
#   left_join(cdf_all, by = 'dyad_id')
# ggplot()+
#   geom_violin(data = cdf_all_dat_long,
#               aes(x = as.factor(age_cat_min),
#                   y = edge_draw,
#                   fill = as.factor(age_cat_max)),
#               #shape = 19,
#               alpha = 0.5)+
#   geom_point(data = cdf_all,
#              aes(x = as.factor(age_cat_min), y = mu_std,
#                  group = as.factor(age_cat_max)),
#              shape = 19,
#              #colour = 'white',
#              size = 1)+
#   scale_x_discrete('age category of younger dyad member')+
#   scale_y_continuous('mean estimated edge weight')+
#   theme(legend.position = 'bottom')+
#   labs(fill = 'age category of older elephant')
# 
## print progress marker
print('edges plotted')

#### fit multivariate Gaussian distribution to output of edge weight model ####
n_windows <- length(unique(cdf_all$window))

## calculate means and covariance matrix per window
logit_edge_draws_mu <- list()
logit_edge_draws_cov <- list()
for(time_window in 1:n_windows){
  window_edges <- sim_dat[,which(cdf_all$window == time_window)]
  logit_edge_draws_mu[[time_window]] <- apply(window_edges, 2, mean)
  # logit_edge_draws_cov[[time_window]] <- cov(window_edges)
  print(paste0('mean calculated for window ',time_window, ' at ',Sys.time()))
}

## parallelise making covariance matrix
library(foreach)
library(doParallel)

#setup parallel backend to use many processors
cores = detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

finalMatrix <- foreach(time_window = 1:n_windows,
                       .combine=cbind) %dopar% {
  tempMatrix = cov(sim_dat[,which(cdf_all$window == time_window)])
  colnames(tempMatrix) = paste0(cdf_all$dyad_id[cdf_all$window == time_window],
                                '_',time_window)
  tempMatrix
                       }

#stop cluster
stopCluster(cl)
print('multivariate Gaussian complete')

## randomly select samples to examine
num_check <- 20
selected_samples <- sample(1:(n_dyads-1), num_check, replace = FALSE)

## set grid layout
rows <- floor(sqrt(num_check))
cols <- ceiling(num_check / rows)
par(mfrow=c(rows, cols), mar=c(2,2,2,1))

## plot
for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values <- rnorm(1e5, mean = mu, sd = sd)
  
  hist(unlist(sim_dat[,i]), probability = TRUE, las = 1,
       main = paste("Dyad", i), xlab = "Value", breaks = 50)
  lines(density(fitted_values), col="blue", lw=1.5)
}

for (i in selected_samples) {
  mu <- logit_edge_draws_mu[i]
  sd <- sqrt(logit_edge_draws_cov[i,i])
  
  fitted_values <- rnorm(1e5, mean = mu, sd = sd)
  
  plot(unlist(sim_dat[,i]), unlist(sim_dat[,i+1]), col = rgb(0,0,1,0.05), las = 1,
       main = paste("cov ", i ,"&",i+1))
}

## reset plot window and clean up
par(mfrow=c(1,1), mai = c(1,1,0.5,0.5))
rm(cols, fitted_values, i, j, mu, num_check, rows, sd, selected_samples) ; gc()

## print progress marker
print('normal approximation complete')

#### prior predictive check ####
cdf_all$age_min_std <- (cdf_all$age_min - mean(cdf_all$age_min)) / sd(cdf_all$age_min)
cdf_all$age_max_std <- (cdf_all$age_max - mean(cdf_all$age_max)) / sd(cdf_all$age_max)
n <- 100
beta_age_min <- rnorm(n, 0, 0.8)
beta_age_max <- rnorm(n, 0, 0.8)
intercept <- rnorm(n, 0, 2)
plot(NULL, las = 1, xlab = 'minimum age', ylab = 'logit edge weight (standardised)',
     ylim = c(min(logit_edge_draws_mu)-5, max(logit_edge_draws_mu)+5),
     xlim = c(min(cdf_all$age_min_std), max(cdf_all$age_max_std)))
abline(h = min(logit_edge_draws_mu), lty = 2)
abline(h = max(logit_edge_draws_mu), lty = 2)
x <- c(min(cdf_all$age_min_std), max(cdf_all$age_max_std))
for(i in 1:n){
  y <- rep(NA, length(x))
  for(j in 1:length(x)){
    y[j] <- intercept[i] + beta_age_min[i]*x[j] + beta_age_max[i]*x[j]
  }
  lines(x = x, y = y, col = rgb(0,0,1,0.4))
}
rm(n, beta_age_max, beta_age_min, intercept, min_dirichlet, max_dirichlet, i, y, j, x) ; gc()

## print progress marker
print('prior predictive complete')

#### fit dyadic regression ####
## create data list
dyad_data <- list(
  # global data size
  num_data = n_data,                                         # number of edges measured
  num_dyads = n_dyads,                                       # number of dyads
  num_nodes = n_nodes,                                       # number of nodes
  num_windows = n_windows,                                   # number of time windows measured
  
  # per time window data size
  num_dyads_window1 = length(which(cdf_all$window == 1)),    # rows in data, window 1
  num_dyads_window2 = length(which(cdf_all$window == 2)),    # rows in data, window 2
  num_dyads_window3 = length(which(cdf_all$window == 3)),    # rows in data, window 3
  
  # number of dyads in all preceding time windows for age indexing
  num_dyads_prev_windows = c(0,                              # Number of rows in data in time window -1
                             length(which(cdf_all$window == 1)),
                             length(which(cdf_all$window == 2))),
  
  # edge means per time window
  logit_edge_mu_1 = logit_edge_draws_mu[[1]],                # Mean logit edge weights, window 1
  logit_edge_mu_2 = logit_edge_draws_mu[[2]],                # Mean logit edge weights, window 2
  logit_edge_mu_3 = logit_edge_draws_mu[[3]],                # Mean logit edge weights, window 3
  
  # covariance matrices per time window
  logit_edge_cov_1 = logit_edge_draws_cov[[1]],              # SD logit edge weights, window 1
  logit_edge_cov_2 = logit_edge_draws_cov[[2]],              # SD logit edge weights, window 2
  logit_edge_cov_3 = logit_edge_draws_cov[[3]],              # SD logit edge weights, window 3
  
  # explanatory variables
  age_min = cdf_all$age_min,                                 # age of younger dyad member
  age_max = cdf_all$age_max,                                 # age of older dyad member
  
  # multimembership terms
  node_1 = cdf_all$node_1,                                   # Node 1 IDs for multimembership terms
  node_2 = cdf_all$node_2,                                   # Node 2 IDs for multimembership terms
  
  # dyad IDs for all time windows
  dyads_window1 = cdf_all$dyad_random[cdf_all$window == 1],  # Dyad IDs, window 1
  dyads_window2 = cdf_all$dyad_random[cdf_all$window == 2],  # Dyad IDs, window 2
  dyads_window3 = cdf_all$dyad_random[cdf_all$window == 3],  # Dyad IDs, window 3
  
)

## print progress marker
print('data list created')

## load dyadic regression model
dyadic_regression <- stan_model('models/dyadic_regression_anp.stan')
#dyadic_regression <- cmdstan_model('models/dyadic_regression_anp.stan')

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
                            chains = n_chains, cores = n_chains)
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
# obtain summary
fit_dyadreg_cdf_all$summary() %>%
  filter(variable %in% c('beta_age_max','beta_age_min','intercept',
                         'delta_min[1]','delta_min[2]','delta_min[3]','delta_min[4]','delta_min[5]',
                         'delta_max[1]','delta_max[2]','delta_max[3]','delta_max[4]','delta_max[5]'))

## extract model fit
summary <- fit_dyadreg_cdf_all$summary()
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract draws
draws <- fit_dyadreg_cdf_all$draws(format = 'df')

## extract dyadic regression slopes
b_max <- draws$beta_age_max
b_min <- draws$beta_age_min
intercept <- draws$intercept
sigma <- draws$sigma
delta_min <- draws[,c('delta_min[1]','delta_min[2]','delta_min[3]','delta_min[4]','delta_min[5]')] ; colnames(delta_min) <- 1:n_age_cat
delta_max <- draws[,c('delta_max[1]','delta_max[2]','delta_max[3]','delta_max[4]','delta_max[5]')] ; colnames(delta_min) <- 1:n_age_cat
delta_j_min <- draws[,c('delta_j_min[1]','delta_j_min[2]','delta_j_min[3]','delta_j_min[4]','delta_j_min[5]','delta_j_min[6]')] ; colnames(delta_j_min) <- 1:(n_age_cat+1)
delta_j_max <- draws[,c('delta_j_max[1]','delta_j_max[2]','delta_j_max[3]','delta_j_max[4]','delta_j_max[5]','delta_j_max[6]')] ; colnames(delta_j_max) <- 1:(n_age_cat+1)
parameters <- data.frame(beta_age_max = b_max,
                         beta_age_min = b_min,
                         intercept = intercept,
                         sigma = sigma) %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = c('beta_age_max','beta_age_min','sigma','intercept'),
               names_to = 'parameter', values_to = 'slope_draw')

## traceplots -- quite wandery, but well mixed with one another -- could just be because there isn't much data?
ggplot(data = parameters)+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_min %>%  as.data.frame() %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = all_of(1:n_age_cat),
               names_to = 'parameter', values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_max %>% as.data.frame() %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = all_of(1:n_age_cat),
               names_to = 'parameter', values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_j_min %>% as.data.frame() %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = all_of(2:(n_age_cat+1)),
               names_to = 'parameter', values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')
delta_j_max %>% as.data.frame(delta_j_max) %>%
  mutate(chain = rep(1:4, each = 1000),
         position = rep(1:1000, 4)) %>%
  pivot_longer(cols = all_of(2:(n_age_cat+1)),
               names_to = 'parameter', values_to = 'slope_draw') %>%
  ggplot()+
  geom_line(aes(x = position, y = slope_draw, colour = as.factor(chain)))+
  theme(legend.position = 'none')+
  scale_colour_viridis_d()+
  facet_wrap(. ~ parameter , scales = 'free_y')

## print progress marker
print('outputs checked')

#### plot predictions ####
## posterior predictive check
plot(density(as.numeric(sim_dat[1, ])),
     main = "Posterior predictive density of edge weights:\nblack = measured edge, red = predicted edge",
     ylim = c(0, 1), col = rgb(0, 0, 0, 0.25), las = 1)
for (i in 1:100) {
  j <- sample(1:1000, 1)

  mu_plot <- rep(NA, n_dyads)
  for(k in 1:n_dyads){
    mu_plot[k] <- intercept[j] + b_min[j]*sum(delta_j_min[j,(1:dyad_data$age_min_cat[k])]) + b_max[j]*sum(delta_j_max[j,(1:dyad_data$age_max_cat[k])])
  }

  sigma_plot <- dyad_data$logit_edge_cov + diag(rep(sigma[j], n_dyads))
  mv_norm <- MASS::mvrnorm(1, mu_plot, sigma_plot)

  lines(density(as.numeric(sim_dat[j, ])), col = rgb(0, 0, 0, 0.25)) # black lines for edge samples
  lines(density(mv_norm), col = rgb(1, 0, 0, 0.25))                  # red lines for predictions

}

save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

## create predictive data frame
pred <- data.frame(age_min = rep(rep(unique(cdf_all$age_min), each = length(unique(cdf_all$age_max))),
                                 n_chains*n_samples),
                   age_max = rep(rep(unique(cdf_all$age_max), length(unique(cdf_all$age_min))),
                                 n_chains*n_samples),
                   intcp = rep(intercept, each = length(unique(cdf_all$age_max))*length(unique(cdf_all$age_min))),
                   beta_min = rep(b_min, each = length(unique(cdf_all$age_max))*length(unique(cdf_all$age_min))),
                   beta_max = rep(b_max, each = length(unique(cdf_all$age_max))*length(unique(cdf_all$age_min))),
                   delta_j_min = NA, delta_j_max = NA) #%>%
#filter(age_min <= age_max)
for(i in 1:(n_samples*n_chains)){
  for(j in 1:ncol(delta_j_max)){
    pred$delta_j_min[pred$age_min == j] <- rowSums(delta_j_min[i,1:j])
    pred$delta_j_max[pred$age_max == j] <- rowSums(delta_j_max[i,1:j])
  }
}

## calculate predictions
pred$pred <- pred$intcp + pred$beta_min * pred$delta_j_min + pred$beta_max * pred$delta_j_max

## convert to unstandardised scale
pred$pred_unstd <- pred$pred*sd(cdf_all$mu) + mean(cdf_all$mu)

## summarise
pred_summary_all <- pred %>%
  group_by(age_min, age_max) %>%
  mutate(pred_lwr = quantile(pred, probs = 0.025),
         pred_mean = mean(pred),
         pred_upr = quantile(pred, probs = 0.975),
         pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
         pred_unstd_mean = mean(pred_unstd),
         pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>%
  ungroup() %>%
  select(age_min, age_max,
         pred_lwr, pred_mean, pred_upr,
         pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>%
  distinct()
pred_summary_possible <- pred %>%
  filter(age_min <= age_max) %>%
  group_by(age_min, age_max) %>%
  mutate(pred_lwr = quantile(pred, probs = 0.025),
         pred_mean = mean(pred),
         pred_upr = quantile(pred, probs = 0.975),
         pred_unstd_lwr = quantile(pred_unstd, probs = 0.025),
         pred_unstd_mean = mean(pred_unstd),
         pred_unstd_upr = quantile(pred_unstd, probs = 0.975)) %>%
  ungroup() %>%
  select(age_min, age_max,
         pred_lwr, pred_mean, pred_upr,
         pred_unstd_lwr, pred_unstd_mean, pred_unstd_upr) %>%
  distinct()

## plot
ggplot()+
  geom_ribbon(data = pred_summary_possible,
              aes(x = age_min,
                  #ymin = pred_lwr, ymax = pred_upr,
                  ymin = invlogit(pred_lwr), ymax = invlogit(pred_upr),
                  group = as.factor(age_max), fill = as.factor(age_max)),
              alpha = 0.3)+
  geom_line(data = pred_summary_possible,
            aes(x = age_min,
                #y = pred_mean,
                y = invlogit(pred_mean),
                colour = as.factor(age_max), group = as.factor(age_max)),
            linewidth = 1)+
  scale_colour_viridis_d(direction = -1)+ scale_fill_viridis_d(direction = -1)+
  geom_point(data = cdf_all, aes(x = age_min, y = invlogit(mu_std), colour = as.factor(age_max)))+
  scale_x_continuous('age of younger dyad member')+
  scale_y_continuous('mean estimated edge weight')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))+
  labs(colour = 'maximum age', fill = 'maximum age')

raw_summary <- sim %>%
  select(age_min, age_max, mu, mu_std) %>%
  group_by(age_min, age_max) %>%
  mutate(mu_mean = mean(mu),
         mu_std_mean = mean(mu_std)) %>%
  ungroup() %>%
  select(-mu, -mu_std) %>%
  distinct() %>%
  mutate(age_pair = paste0(age_min, '_', age_max))
compare <- pred_summary_possible %>%
  mutate(age_pair = paste0(age_min, '_', age_max)) %>%
  left_join(raw_summary[,c('mu_mean','mu_std_mean','age_pair')], by = 'age_pair') %>%
  rename(raw_unstd = mu_mean, raw_std = mu_std_mean) %>%
  select(age_min, age_max, raw_unstd, pred_unstd_mean, raw_std, pred_mean, pred_lwr, pred_unstd_lwr, pred_upr, pred_unstd_upr)
ggplot(compare)+
  geom_ribbon(aes(x = raw_unstd,
                  ymin = pred_unstd_lwr, ymax = pred_unstd_upr),
              alpha = 0.3)+
  geom_line(aes(x = raw_unstd, y = pred_unstd_mean),
            linewidth = 1)+
  geom_line(data = data.frame(x = c(min(cdf_all$mu),max(cdf_all$mu)),
                              y = c(min(cdf_all$mu),max(cdf_all$mu))),
            aes(x = x, y = y),
            linewidth = 0.5, colour = 'red')+
  scale_x_continuous('raw mean')+
  scale_y_continuous('predicted mean')+
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

dev.off()
save.image('step5_dyadicregression/anp_dyadic_simulation.RData')

## print progress marker
print('predictions complete')
