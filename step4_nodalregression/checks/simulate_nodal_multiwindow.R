#### set up #### 
library(LaplacesDemon)
library(tidyverse)
library(cmdstanr) #set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

theme_set(theme_bw())

#### simulate population ####
## define population parameters
min_age <- 11                                                # youngest individual
max_age <- 60                                                #  oldest  individual
n_nodes <- ((max_age+1) - min_age)                           # total nodes = 2 per age

## simulate population for first time window
sim1 <- data.frame(node = 1:n_nodes,                          # create data frame of individuals
                  age = c(sample(min_age:max_age, n_nodes,
                                 prob = 1/(min_age:max_age), replace = T)), # more at lower ages
                  mu = NA, sd = NA)

## simulate population for second time window -- half the nodes resampled from window 1, half new
sim2 <- data.frame(node = c(sample(1:n_nodes, n_nodes/2, replace = F),
                            (n_nodes+1):(n_nodes*2)),                          # create data frame of individuals
                   age = NA,
                   mu = NA, sd = NA)
for(i in 1:nrow(sim2)){
  if(sim2$node[i] %in% sim1$node){
    sim2$age[i] <- sim1$age[sim1$node == sim2$node[i]] + 2 # if present in 1st window, age = age in t1 + 2yrs
  } else {
    sim2$age[i] <- sample(min_age:max_age, 1, prob = 1/(min_age:max_age)) # if not present, randomly assign a new age, more likely to be young than old
  }
}

## simulate population for third time window -- as above, half the nodes are resampled from oervious two time windows, so some may be represented in all three windows, some twice, and some only once
n_nodes <- length(unique(sim2$node))
sim3 <- data.frame(node = c(sample(1:n_nodes, n_nodes/2, replace = F),
                            (n_nodes+1):(n_nodes*2)),                          # create data frame of individuals
                   age = NA,
                   mu = NA, sd = NA)
for(i in 1:nrow(sim3)){
  if(sim3$node[i] %in% sim2$node){
    sim3$age[i] <- sim2$age[sim2$node == sim3$node[i]] + 2   # if present in 2nd window, age = age in t2 + 2yrs
  } else {
    if(sim3$node[i] %in% sim1$node){
      sim3$age[i] <- sim1$age[sim1$node == sim3$node[i]] + 4 # if present in 1st window, age = age in t1 + 4yrs
    } else {
      sim3$age[i] <- sample(min_age:max_age, 1, prob = 1/(min_age:max_age)) # if new, randomly assign a new age, more likely to be young than old
    }
  }
}

## combine populations into single data frame
sim1$window <- 1 ; sim2$window <- 2 ; sim3$window <- 3 # create window ID variable
sim <- rbind(sim1, sim2, sim3)
sim$node_window <- paste0(sim$node, '_', sim$window)
sim$age_std <- ( sim$age - mean(sim$age) ) / sd(sim$age)
rm(sim1, sim2, sim3) ; gc()

## redefine parameters
n_data <- nrow(sim)
n_nodes <- length(unique(sim$node))
n_windows <- length(unique(sim$window))

## view ages
hist(sim$age)

#### simulate centralities ####
## simulate age effect
sim_slope <- -1

## simulate intercept
sim_intcp <- 2

## simulate random effects -- these are just random draws from a distribution centred on zero, but they do not have to have a mean of 0 between them. Should they? (i.e. if random effect of windows 1 and 2 are both positive, does window 3 need to be negative to counteract it?)
sim_window_unique <- rnorm(n_windows, mean = 0, sd = 1)
sim_node_unique <- rnorm(n_nodes, mean = 0, sd = 0.2)

## calculate mean centrality per elephant
sim$mu <- sim$age_std * sim_slope + sim_intcp + sim_window_unique[sim$window] + sim_node_unique[sim$node]    # simulate mean centrality on outcome scale -- NOTE: THIS IS PRODUCING AN EFFECT OF 1 INCREASE IN AGE STANDARD DEVIATION = -1 IN CENTRALITY OUTCOME SCALE. AGE IS STANDARDISED ALREADY, CENTRALITY IS NOT.

## plot mean centrality against age
plot(mu ~ age, data = sim[sim$window == 1,], col = 'red', pch = 19, ylim = c(0,5), las = 1)            # plot
points(mu ~ age, data = sim[sim$window == 2,], col = 'blue', pch = 19)                                 # plot
points(mu ~ age, data = sim[sim$window == 3,], col = 'green', pch = 19)                                # plot

# standardise mean centrality
sim$mu_std <- ( sim$mu - mean(sim$mu) ) / sd(sim$mu)

## plot mean standardised centrality against age
plot(mu_std ~ age_std, data = sim[sim$window == 1,], col = 'red', pch = 19, ylim = c(-3,3), las = 1)   # plot
points(mu_std ~ age_std, data = sim[sim$window == 2,], col = 'blue', pch = 19)                         # plot
points(mu_std ~ age_std, data = sim[sim$window == 3,], col = 'green', pch = 19)                        # plot

## simulate full distribution of samples per node
sim$sd <- 1         # can be made to vary amongst nodes, currently all elephants have equal variance in centrality (not realistic given that in the real data some are seen more often than others)
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_data, dimnames = list(NULL, sim$node_window))    # create matrix for full centrality distribution
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution for each elephant
}

## plot full unstandardised distribution
plot(sim_dat[1,which(sim$window == 1)] ~ sim$age[which(sim$window == 1)],
     col = 'red', pch = 19, ylim = c(-1,6))       # plot simulated values against age for window 1 (first row of simulated centralities only)
points(sim_dat[1,which(sim$window == 2)] ~ sim$age[which(sim$window == 2)],
       col = 'blue', pch = 19)                    # plot simulated values against age for window 2 (first row of simulated centralities only)
points(sim_dat[1,which(sim$window == 3)] ~ sim$age[which(sim$window == 3)],
       col = 'green', pch = 19)                   # plot simulated values against age for window 3 (first row of simulated centralities only)

## standardise full set of centralities
sim_dat_std <- sim_dat                            # create matrix to fill
for(i in 1:nrow(sim_dat_std)){
  sim_dat_std[i,] <- (sim_dat[i,] - mean(sim_dat[i,]) ) / sd(sim_dat[i,]) # standardise values within each row of the matrix (1 row = 1 network -- taken from Jordan's code)
}

## plot full standardised distribution
plot(sim_dat_std[1,which(sim$window == 1)] ~ sim$age[which(sim$window == 1)],
     col = 'red', pch = 19, ylim = c(-5,5))       # plot simulated values against age for window 1
points(sim_dat_std[1,which(sim$window == 2)] ~ sim$age[which(sim$window == 2)],
       col = 'blue', pch = 19)                    # plot simulated values against age for window 2
points(sim_dat_std[1,which(sim$window == 3)] ~ sim$age[which(sim$window == 3)],
       col = 'green', pch = 19)                   # plot simulated values against age for window 3

## visualise
df_plot <- data.frame(sim_dat_std) %>% 
  pivot_longer(cols = everything(),
               names_to = "node", values_to = "centrality") %>% 
  separate(node, into = c('X','node_window'), remove = T, sep = 1) %>% 
  dplyr::select(-X) %>% 
  separate(node_window, into = c('node','window'), sep = '_', remove = F) %>% 
  mutate(node = as.integer(node)) %>% 
  left_join(sim[,c('node_window','age')], by = 'node_window') %>% 
  filter(node %in% seq(1, n_data, by = 5))
ggplot(data = df_plot, aes(x = centrality, fill = age, group = node_window)) +
  geom_density(linewidth = 0.4, alpha = 0.6) +
  facet_grid(rows = vars(as.factor(node)), scales = "free") +
  labs(x = "Eigenvector centrality (standardised)") + 
  theme_void() + 
  theme(strip.text.y = element_text(size = 12),
        axis.text.x = element_text(angle = 0, size = 12, debug = FALSE),
        axis.title.x = element_text(size = 12),
        plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))

#### prep inputs ####
# compute normal approximation by window -- calculate means per node
sim_cent_mu_1 <- apply(sim_dat_std[,which(sim$window == 1)], 2, mean)
sim_cent_mu_2 <- apply(sim_dat_std[,which(sim$window == 2)], 2, mean)
sim_cent_mu_3 <- apply(sim_dat_std[,which(sim$window == 3)], 2, mean)

# compute normal approximation by window -- calculate covariance matrix
sim_cent_cov_1 <- cov(sim_dat_std[,which(sim$window == 1)])
sim_cent_cov_2 <- cov(sim_dat_std[,which(sim$window == 2)])
sim_cent_cov_3 <- cov(sim_dat_std[,which(sim$window == 3)])

## check normal approximation -- simulate from combined mean and covariance, plot curve against relative node ID 
sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_1, sim_cent_cov_1)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 1),1)
plot(density(sim_dat_std[,node_id_sample]), lwd = 2, las = 1,                         # plot true density curve
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_2, sim_cent_cov_2)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 2),1)
plot(density(sim_dat_std[,node_id_sample]),          # plot true density curve
     lwd = 2, las = 1,
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample - length(which(sim$window == 1))]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

sim_cent_samples <- MASS::mvrnorm(1e5, sim_cent_mu_3, sim_cent_cov_3)     # simulate from multivariate normal
node_id_sample <- sample(which(sim$window == 3),1)
plot(density(sim_dat_std[,node_id_sample]),          # plot true density curve
     lwd = 2, las = 1,
     main = "Estimated standardised centrality vs normal approximation", xlab = "Logit edge weight")
lines(density(sim_cent_samples[, node_id_sample - length(which(sim$window < 3))]), col = rgb(0,0,1,0.5), lwd = 2)      # overlay normal approximation

rm(sim_cent_samples, node_id_sample) ; gc()

## create data
eigen_list <- list(num_data = n_data,
                   num_nodes = n_nodes,
                   num_windows = length(unique(sim$window)),
                   num_nodes_window1 = length(which(sim$window == 1)),
                   num_nodes_window2 = length(which(sim$window == 2)),
                   num_nodes_window3 = length(which(sim$window == 3)),
                   centrality_mu_1 = sim_cent_mu_1,
                   centrality_mu_2 = sim_cent_mu_2,
                   centrality_mu_3 = sim_cent_mu_3,
                   centrality_cov_1 = sim_cent_cov_1,
                   centrality_cov_2 = sim_cent_cov_2,
                   centrality_cov_3 = sim_cent_cov_3,
                   node_age = sim$age_std,
                   window = sim$window,
                   nodes = sim$node,
                   nodes_window1 = sim$node[sim$window == 1],
                   nodes_window2 = sim$node[sim$window == 2],
                   nodes_window3 = sim$node[sim$window == 3])

## check inputs
plot(eigen_list$centrality_mu_1 ~ eigen_list$node_age[eigen_list$window == 1], pch = 19, col = 'red')
points(eigen_list$centrality_mu_2 ~ eigen_list$node_age[eigen_list$window == 2], pch = 19, col = 'blue')
points(eigen_list$centrality_mu_3 ~ eigen_list$node_age[eigen_list$window == 3], pch = 19, col = 'green')

#### prior predictive check ####
n <- 100
beta_age <- rnorm(n, 0, 0.8)
intercept  <- rnorm(n, 0, 0.8)
plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector (standardised)',
     ylim = c(min(c(sim_cent_mu_1, sim_cent_mu_2))-2, max(c(sim_cent_mu_1, sim_cent_mu_2))+2),
     xlim = c(min(sim$age_std), max(sim$age_std)))
abline(h = min(c(sim_cent_mu_1, sim_cent_mu_2, sim_cent_mu_3)), lty = 2)
abline(h = max(c(sim_cent_mu_1, sim_cent_mu_2, sim_cent_mu_3)), lty = 2)
for(i in 1:n){
  lines(x = seq(min(sim$age_std), max(sim$age_std), length.out = 2),
        y = intercept[i] + beta_age[i]*c(min(sim$age_std), max(sim$age_std)),
        col = rgb(0,0,1,0.4))
} # looks mad when I set true slope very small, but decent when larger -- given that I have no idea what the real version should be, I've kept them relatively weak

#### run model -- age as a continuous variable with 2 windows ####
## load model
#nodal_regression <- stan_model('models/eigen_regression_intercept.stan')
nodal_regression <- cmdstan_model('models/eigen_regression_combinewindows.stan')

## run model
n_chains <- 4
n_samples <- 1000
#fit_sim <- sampling(nodal_regression, data = eigen_list, chains = n_chains, cores = n_chains)
fit_sim <- nodal_regression$sample(data = eigen_list,
                                   chains = n_chains, parallel_chains = n_chains,
                                   iter_warmup = n_samples, iter_sampling = n_samples)

#### check outputs ####
## extract model fit -- very good!
fit_sim$summary()
summary <- fit_sim$summary()
par(mfrow = c(3,1))
hist(summary$rhat, breaks = 50)
hist(summary$ess_bulk, breaks = 50)
hist(summary$ess_tail, breaks = 50)
par(mfrow = c(1,1))

## extract posterior
#params <- rstan::extract(fit_sim)
params <- fit_sim$draws(format = 'draws_df')
colnames(params)

## separate random effects from global parameters
rand_window <- params %>% 
  dplyr::select(`rand_window[1]`,`rand_window[2]`,`rand_window[3]`)
rand_node <- params[,8:(n_nodes+7)]

## traceplot all parameters
#traceplot(fit_sim, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
params %>% 
  select(intercept,beta_age,sigma,
         `rand_node[1]`,`rand_node[50]`,`rand_node[100]`,
         `rand_window[1]`,`rand_window[2]`,
         `predictor_window1[1]`,`predictor_window1[16]`,`predictor_window1[38]`,`predictor_window1[50]`,
         `predictor_window2[1]`,`predictor_window2[25]`,`predictor_window2[50]`,`predictor_window2[75]`) %>% 
  pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
  mutate(chain_position = rep(rep(1:n_samples, each = 16,
  ), n_chains),
  chain = rep(1:n_chains, each = 16*n_samples)) %>% 
  #filter(chain == 4) %>% # inspect individual chains -- some are really bad and haven't explored at all or are very wandery
  ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
  geom_line()+
  facet_wrap(. ~ parameter, scales = 'free_y')+
  theme_bw()+
  theme(legend.position = 'none')

#### posterior predictive check ####
par(mfrow = c(3,1))
plot(density(sim_dat_std[1, which(sim$window == 1)]), las = 1, ylim = c(0,1),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data1 <- list(num_nodes_window1 = length(which(sim$window == 1)),
                    centrality_mu_1 = sim_cent_mu_1,
                    centrality_cov_1 = sim_cent_cov_1,
                    node_age = sim$age_std[sim$window == 1],
                    window = 1,
                    nodes = sim$node[sim$window == 1],
                    nodes_window1 = sim$node[sim$window == 1])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, which(sim$window == 1)]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data1$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data1$window]) + as.numeric(rand_node[j,eigen_data1$nodes[k]])
  }
  sigma <- sim_cent_cov_1 + diag(rep(params$sigma[j], eigen_list$num_nodes_window1))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

plot(density(sim_dat_std[1, which(sim$window == 2)]), las = 1, ylim = c(0,1),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data2 <- list(num_nodes_window2 = length(which(sim$window == 2)),
                    centrality_mu_2 = sim_cent_mu_2,
                    centrality_cov_2 = sim_cent_cov_2,
                    node_age = sim$age_std[sim$window == 2],
                    window = 2,
                    nodes = sim$node[sim$window == 2],
                    nodes_window2 = sim$node[sim$window == 2])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, which(sim$window == 2)]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data2$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data2$window]) + as.numeric(rand_node[j,eigen_data2$nodes[k]])
  }
  sigma <- sim_cent_cov_2 + diag(rep(params$sigma[j], eigen_list$num_nodes_window2))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}

plot(density(sim_dat_std[1, which(sim$window == 3)]), las = 1, ylim = c(0,1),
     main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
     col=rgb(0, 0, 0, 0.25))
eigen_data3 <- list(num_nodes_window3 = length(which(sim$window == 3)),
                    centrality_mu_3 = sim_cent_mu_3,
                    centrality_cov_3 = sim_cent_cov_3,
                    node_age = sim$age_std[sim$window == 3],
                    window = 3,
                    nodes = sim$node[sim$window == 3],
                    nodes_window3 = sim$node[sim$window == 3])
for (i in 1:100) {
  j <- sample(1:length(params$beta_age), 1)
  lines(density(sim_dat_std[j, which(sim$window == 3)]), col=rgb(0, 0, 0, 0.25))
  mu <- params$beta_age[j]*eigen_data3$node_age + params$intercept[j]
  for(k in 1:length(mu)) {
    mu[k] <- mu[k] + as.numeric(rand_window[j,eigen_data3$window]) + as.numeric(rand_node[j,eigen_data3$nodes[k]])
  }
  sigma <- sim_cent_cov_3 + diag(rep(params$sigma[j], eigen_list$num_nodes_window3))
  lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
}
par(mfrow = c(1,1))

#### predict from model -- predict from raw data, not from hypothetical/counterfactual data ####
## create empty matrix to fill with predictions of mean centrality per node
mu_std <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))

## populate matrix = mean centrality values per node, predicting for real data
for(i in 1:nrow(mu_std)){
  mu_std[i,] <- params$beta_age[i] * sim$age_std + params$intercept[i] + as.numeric(rand_window[i,sim$window]) + as.numeric(rand_node[i, sim$node])
}

## add mean and CI of predicted means to input data frame for comparison
sim$mu_mean_std <- apply(mu_std, 2, mean)
sim$mu_lwr_std <- NA ; sim$mu_upr_std <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr_std[i] <- rethinking::HPDI(mu_std[,i], prob = 0.95)[1]
  sim$mu_upr_std[i] <- rethinking::HPDI(mu_std[,i], prob = 0.95)[2]
}

## plot mean of model vs mean of raw data
ggplot()+
  geom_point(data = sim, aes(x = mu_std, y = mu_mean_std, colour = as.factor(window)))+
  scale_colour_viridis_d()+
  labs(colour = 'window', x = 'simulated mean (standardised)', y = 'predicted mean (standardised)')+
  geom_abline(slope = 1, intercept = 0) # add line showing where points would lie if model fit was perfect

## put together sigma arrays, separated by time window
sigma1 <- array(NA, dim = c(eigen_list$num_nodes_window1, eigen_list$num_nodes_window1, nrow(params)),
                dimnames = list(eigen_list$nodes_window1, eigen_list$nodes_window1, NULL))
sigma2 <- array(NA, dim = c(eigen_list$num_nodes_window2, eigen_list$num_nodes_window2, nrow(params)),
                dimnames = list(eigen_list$nodes_window2, eigen_list$nodes_window2, NULL))
sigma3 <- array(NA, dim = c(eigen_list$num_nodes_window3, eigen_list$num_nodes_window3, nrow(params)),
                dimnames = list(eigen_list$nodes_window3, eigen_list$nodes_window3, NULL))
for(i in 1:nrow(params)){
  sigma1[,,i] <- sim_cent_cov_1 + diag(rep(params$sigma[i], eigen_list$num_nodes_window1))
  sigma2[,,i] <- sim_cent_cov_2 + diag(rep(params$sigma[i], eigen_list$num_nodes_window2))
  sigma3[,,i] <- sim_cent_cov_3 + diag(rep(params$sigma[i], eigen_list$num_nodes_window3))
}

## create empty matrix to take full set of predicted values per elephant
predictions_std <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))

## populate matrix using mean values in matrix mu_std, and sigma values based on time window
for(i in 1:nrow(predictions)){
  predictions_std[i,sim$window == 1] <- MASS::mvrnorm(1, mu_std[i,sim$window == 1], sigma1[,,i])
  predictions_std[i,sim$window == 2] <- MASS::mvrnorm(1, mu_std[i,sim$window == 2], sigma2[,,i])
  predictions_std[i,sim$window == 3] <- MASS::mvrnorm(1, mu_std[i,sim$window == 3], sigma3[,,i])
}

## add CI of predicted data points to input data frame for comparison
sim$predict_lwr_std <- NA ; sim$predict_upr_std <- NA
for( i in 1:nrow(sim) ){
  sim$predict_lwr_std[i] <- rethinking::HPDI(predictions_std[,i], prob = 0.95)[1]
  sim$predict_upr_std[i] <- rethinking::HPDI(predictions_std[,i], prob = 0.95)[2]
}

## plot predictions
ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr_std, ymax = predict_upr_std, fill = as.factor(window)),
              alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
  geom_ribbon(aes(x = age, ymin = mu_lwr_std, ymax = mu_upr_std, fill = as.factor(window)),
              alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
  geom_line(aes(x = age, y = mu_mean_std, colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(aes(x = age, y = mu_std))+          # original data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

#### convert predictions to unstandardised scale ####
## revert predicted means and full predictions to unstandardised scale
mu_ustd <- (mu_std * sd(sim$mu)) + mean(sim$mu)
predictions_ustd <- (predictions * sd(sim$mu)) + mean(sim$mu)

## insert all means and CIs to original data frame for comparison
sim$mu_mean_ustd <- apply(mu_ustd, 2, mean)
sim$mu_lwr_ustd <- NA ; sim$mu_upr_ustd <- NA
sim$predict_lwr_ustd <- NA ; sim$predict_upr_ustd <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr_ustd[i] <- rethinking::HPDI(mu_ustd[,i], prob = 0.95)[1]
  sim$mu_upr_ustd[i] <- rethinking::HPDI(mu_ustd[,i], prob = 0.95)[2]
  sim$predict_lwr_ustd[i] <- rethinking::HPDI(predictions_ustd[,i], prob = 0.95)[1]
  sim$predict_upr_ustd[i] <- rethinking::HPDI(predictions_ustd[,i], prob = 0.95)[2]
}

## plot on unstandardised scale
ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr_ustd, ymax = predict_upr_ustd, fill = as.factor(window)),
              alpha = 0.2)+                      # background layer showing the 95% CI of all predictions
  geom_ribbon(aes(x = age, ymin = mu_lwr_ustd, ymax = mu_upr_ustd, fill = as.factor(window)),
              alpha = 0.4)+                      # mid layer showing the 95% CI of predicted means
  geom_line(aes(x = age, y = mu_mean_ustd, colour = as.factor(window)))+  # line showing mean of predicted means
  geom_point(aes(x = age, y = mu))+              # original data points (standardised centrality, actual age)
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+             # separate plots per window
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

#### extract original values from output -- simulated slope value originally used produces an effect on the unstandardised scale. The model works on the standardised scale. Convert predictions to unstandardised scale and then run contrasts to calculate the slope coefficient. ####
## predict means from model again, now using age_std + 1 stdev -- for now working with age_std+1SD not age+1yr because that's how I originally simulated the data. if change the simulation so that mean centrality = age*slope rather age_std*slope, then change this  to extract the correct slope value
mu2_std <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(mu2_std)){
  mu2_std[i,] <- params$beta_age[i] * (sim$age_std+1) + params$intercept[i] + as.numeric(rand_window[i,sim$window]) + as.numeric(rand_node[i, sim$node])
}
sim$mu_mean_plus1yr_std <- apply(mu2, 2, mean)

## full distribution of predictions using age_std + 1 stdev
predictions2_std <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(predictions2_std)){
  predictions2_std[i,sim$window == 1] <- MASS::mvrnorm(1, mu2_std[i,sim$window == 1], sigma1[,,i])
  predictions2_std[i,sim$window == 2] <- MASS::mvrnorm(1, mu2_std[i,sim$window == 2], sigma2[,,i])
  predictions2_std[i,sim$window == 3] <- MASS::mvrnorm(1, mu2_std[i,sim$window == 3], sigma3[,,i])
}

## contrast predictions on standardised scale -- check that this returns the marginal effect presented in the summary
contrast_std <- predictions2_std - predictions_std    # contrast between predicted values for raw data and all same data but add 1 stdev to age
head(contrast_std[,1:5])                              # check matrix looks right
mean(params$beta_age)                                 # output parameter direct from model
mean(contrast_std)                                    # very similar to mean of output parameter
quantile(contrast_std, prob = c(0.025, 0.975))        # very wide

## revert predictions to unstandardised scale

## contrast predictions on unstandardised scale. Does this work when sim$mu for age = x isn't the same as it would be for x+1? -- HAVE USED TWO ALTERNATIVE METHODS TO CALCULATE THE CONTRAST MATRIX. THEY GIVE TWO COMPLETELY DIFFERENT VALUES FOR SLOPE, NEITHER OF WHICH IS THE SAME AS THE ORIGINAL INPUT PARAMETER
sim_slope # value that SHOULD be produced

# METHOD 1: unstandardise predictions2 matrix, then subtract predictions_ustd from predictions2_ustd
predictions2_ustd <- predictions2_std * sd(sim$mu) + mean(sim$mu)
contrast_ustd <- predictions2_ustd - predictions_ustd
head(contrast_ustd[,1:5])                             # check matrix looks right
mean(contrast_ustd)                                   # not even close to original parameter
quantile(contrast_ustd, prob = c(0.025, 0.975))       # does include original parameter, but extremely wide

# METHOD 2: unstandardise contrast matrix directly
contrast_ustd <- contrast_std * sd(sim$mu) + mean(sim$mu)
head(contrast_ustd[,1:5])                             # check matrix looks right
mean(contrast_ustd)                                   # not even close to original parameter
quantile(contrast_ustd, prob = c(0.025, 0.975))       # does include original parameter, but extremely wide
