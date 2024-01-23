## set up #### 
library(LaplacesDemon)
library(tidyverse)
library(cmdstanr)

#set_cmdstan_path('../packages/.cmdstan/cmdstan-2.31.0/')

## simulate population ####
## define population parameters
min_age <- 11                                                # youngest individual
max_age <- 60                                                #  oldest  individual
n_nodes <- ((max_age+1) - min_age)                           # total nodes = 2 per age

## simulate population for first time window
sim <- data.frame(node = 1:n_nodes,                          # create data frame of individuals
                  age = c(sample(min_age:max_age, n_nodes, prob = 1/(min_age:max_age), # more at lower ages
                                 replace = T)),
                  mu = NA, sd = NA)

## redefine parameters
n_data <- nrow(sim)
n_nodes <- length(unique(sim$node))

## add additional time windows ####
## simulate population for first time window
sim1 <- sim

## simulate population for second time window -- half the nodes resampled from window 1, half new
sim2 <- data.frame(node = c(sample(1:n_nodes, n_nodes/2, replace = F),
                            (n_nodes+1):(n_nodes*2)),                          # create data frame of individuals
                   age = NA,
                   mu = NA, sd = NA)
for(i in 1:nrow(sim2)){
  if(sim2$node[i] %in% sim1$node){
    sim2$age[i] <- sim1$age[sim1$node == sim2$node[i]] + 2 # if present in 1st window, age = age in t1 + 2 y
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
    sim3$age[i] <- sim2$age[sim2$node == sim3$node[i]] + 2
  } else {
    if(sim3$node[i] %in% sim1$node){
      sim3$age[i] <- sim1$age[sim1$node == sim3$node[i]] + 2
    } else {
      sim3$age[i] <- sample(min_age:max_age, 1, prob = 1/(min_age:max_age)) #(min_age:max_age)[i-(n_nodes/2)]
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

## simulate centralities ####
## simulate age effect
sim_slope <- -1                     # set age effect -- smaller = bigger impact on invlogit scale as values large
sim_intcp <- 2
sim_window_unique <- rnorm(n_windows, mean = 0, sd = 1)
sim_node_unique <- rnorm(n_nodes, mean = 0, sd = 0.2)
sim$mu <- sim$age_std * sim_slope + sim_intcp + sim_window_unique[sim$window] + sim_node_unique[sim$node]    # simulate mean centrality on normal scale
plot(mu ~ age, data = sim[sim$window == 1,], col = 'red', pch = 19, ylim = c(0,5))            # plot
points(mu ~ age, data = sim[sim$window == 2,], col = 'blue', pch = 19)                        # plot
points(mu ~ age, data = sim[sim$window == 3,], col = 'green', pch = 19)                       # plot
sim$mu_std <- ( sim$mu - mean(sim$mu) ) / sd(sim$mu)
plot(mu_std ~ age_std, data = sim[sim$window == 1,], col = 'red', pch = 19, ylim = c(-5,5))   # plot
points(mu_std ~ age_std, data = sim[sim$window == 2,], col = 'blue', pch = 19)                # plot
points(mu_std ~ age_std, data = sim[sim$window == 3,], col = 'green', pch = 19)               # plot

## simulate full distribution of samples per node
sim$sd <- 1
sim_dat <- matrix(data = NA, nrow = 4000, ncol = n_data, dimnames = list(NULL, sim$node_window))    # create matrix
for(j in 1:n_data){
  sim_dat[,j] <- rnorm(n = nrow(sim_dat), mean = sim$mu[j], sd = sim$sd[j])  # simulate distribution
}
plot(sim_dat[1,1:(length(which(sim$window == 1)))] ~ sim$age[1:(length(which(sim$window == 1)))],
     col = 'red', pch = 19, ylim = c(-1,6))       # plot simulated values against age for window 1
points(sim_dat[1,(length(which(sim$window == 1))+1):(length(which(sim$window == 1 | sim$window == 2)))] ~ sim$age[(length(which(sim$window == 1))+1):(length(which(sim$window == 1 | sim$window == 2)))],
       col = 'blue', pch = 19)    # plot simulated values against age for window 2
points(sim_dat[1,(length(which(sim$window == 1 | sim$window == 2))+1):nrow(sim)] ~ sim$age[(length(which(sim$window == 1 | sim$window == 2))+1):nrow(sim)],
       col = 'green', pch = 19)   # plot simulated values against age for window 3

## standardise
sim_dat_std <- sim_dat               # create matrix to fill
for(i in 1:nrow(sim_dat_std)){
  sim_dat_std[i,] <- (sim_dat[i,] - mean(sim_dat[i,]) ) / sd(sim_dat[i,]) # standardise values
}
plot(sim_dat_std[1,which(sim$window == 1)] ~ sim$age[which(sim$window == 1)],
     col = 'red', pch = 19, ylim = c(-5,5))       # plot simulated values against age for window 1
points(sim_dat_std[1,which(sim$window == 2)] ~ sim$age[which(sim$window == 2)],
       col = 'blue', pch = 19)    # plot simulated values against age for window 2
points(sim_dat_std[1,which(sim$window == 3)] ~ sim$age[which(sim$window == 3)],
       col = 'green', pch = 19)   # plot simulated values against age for window 3

## visualise
df_wide <- data.frame(sim_dat_std)
df_plot <- df_wide %>% 
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

## normal approximation ####
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

## prep inputs ####
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

## prior predictive check ####
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

## run model -- age as a continuous variable with 2 windows ####
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

## check outputs ####
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
rand_window <- params %>% 
  dplyr::select(`rand_window[1]`,`rand_window[2]`,`rand_window[3]`)
rand_node <- params[,8:(n_nodes+7)]

## traceplot linear effect size
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

## posterior predictive check
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

## predict from model -- predict from raw data, not from hypothetical/counterfactual data ####
mu <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(mu)){
  mu[i,] <- params$beta_age[i] * sim$age_std + params$intercept[i] + as.numeric(rand_window[i,sim$window]) + as.numeric(rand_node[i, sim$node])
}
sim$mu_mean <- apply(mu, 2, mean)

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

predictions <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(predictions)){
  predictions[i,sim$window == 1] <- MASS::mvrnorm(1, mu[i,sim$window == 1], sigma1[,,i])
  predictions[i,sim$window == 2] <- MASS::mvrnorm(1, mu[i,sim$window == 2], sigma2[,,i])
  predictions[i,sim$window == 3] <- MASS::mvrnorm(1, mu[i,sim$window == 3], sigma3[,,i])
}

pred_data <- sim
pred_data$mu_model <- NA
pred_data$prediction <- NA
for(i in 1:nrow(pred_data)){
  pred_data$mu_model[i] <- list(mu[,i])
  pred_data$prediction[i] <- list(predictions[,i])
}
pred_data <- unnest(cols = prediction, data = pred_data)

sim$mu_lwr <- NA ; sim$mu_upr <- NA
sim$predict_lwr <- NA ; sim$predict_upr <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr[i] <- rethinking::HPDI(mu[,i], prob = 0.95)[1]
  sim$mu_upr[i] <- rethinking::HPDI(mu[,i], prob = 0.95)[2]
  sim$predict_lwr[i] <- rethinking::HPDI(predictions[,i], prob = 0.95)[1]
  sim$predict_upr[i] <- rethinking::HPDI(predictions[,i], prob = 0.95)[2]
}

## plot mean of model vs mean of raw data
plot(sim$mu_mean ~ sim$mu_std, pch = 19, col = rgb(0,0,1,0.2)) ; abline(a = 0, b = 1, lwd = 2)

## plot all predictions with points = mean of raw data
ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr, ymax = predict_upr, fill = as.factor(window)),
              alpha = 0.2)+
  geom_ribbon(aes(x = age, ymin = mu_lwr, ymax = mu_upr, fill = as.factor(window)),
              alpha = 0.4)+
  geom_line(aes(x = age, y = mu_mean, colour = as.factor(window)))+
  geom_point(aes(x = age, y = mu_std, colour = as.factor(window)))+
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## plot on unstandardised scale ####
## revert predictions to unstandardised scale
mu_ustd <- (mu * sd(sim$mu)) + mean(sim$mu)
predictions_ustd <- (predictions * sd(sim$mu)) + mean(sim$mu)

sim$mu_mean_ustd <- apply(mu_ustd, 2, mean)
sim$mu_lwr_ustd <- NA ; sim$mu_upr_ustd <- NA
sim$predict_lwr_ustd <- NA ; sim$predict_upr_ustd <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr_ustd[i] <- rethinking::HPDI(mu_ustd[,i], prob = 0.95)[1]
  sim$mu_upr_ustd[i] <- rethinking::HPDI(mu_ustd[,i], prob = 0.95)[2]
  sim$predict_lwr_ustd[i] <- rethinking::HPDI(predictions_ustd[,i], prob = 0.95)[1]
  sim$predict_upr_ustd[i] <- rethinking::HPDI(predictions_ustd[,i], prob = 0.95)[2]
}

ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr_ustd, ymax = predict_upr_ustd, fill = as.factor(window)),
              alpha = 0.2)+
  geom_ribbon(aes(x = age, ymin = mu_lwr_ustd, ymax = mu_upr_ustd, fill = as.factor(window)),
              alpha = 0.4)+
  geom_line(aes(x = age, y = mu_mean_ustd, colour = as.factor(window)))+
  geom_point(aes(x = age, y = mu), pch = 19, size = 0.8)+
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(as.factor(node) ~ as.factor(window))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## predict again, but first convert params to unstandardised scale
# beta_slope_ustd = beta_slope_std * ( sd(sim$mu)/sd(sim$age) )
# intcp_ustd = mean(sim$mu) - beta_slope_ustd * ( mean(sim$age)/sd(sim$age) )
# sigma_ustd = sigma * sd(sim$mu)

params_ustd <- params[,c('beta_age','intercept','sigma')]
params_ustd$beta_age <- params$beta_age * ( sd(sim$mu)/sd(sim$age) )
params_ustd$intercept <- mean(sim$mu) - params_ustd$beta_age * ( mean(sim$age)/sd(sim$age) )
params_ustd$sigma <- params_ustd$sigma * sd(sim$mu)

mu <- matrix(NA, nrow = nrow(params_ustd), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(mu)){
  mu[i,] <- params_ustd$beta_age[i] * sim$age + params_ustd$intercept[i] + as.numeric(rand_window[i,sim$window]) + as.numeric(rand_node[i, sim$node])
}
sim$mu_mean <- apply(mu, 2, mean)

sigma1 <- array(NA, dim = c(eigen_list$num_nodes_window1, eigen_list$num_nodes_window1, nrow(params_ustd)),
                dimnames = list(eigen_list$nodes_window1, eigen_list$nodes_window1, NULL))
sigma2 <- array(NA, dim = c(eigen_list$num_nodes_window2, eigen_list$num_nodes_window2, nrow(params_ustd)),
                dimnames = list(eigen_list$nodes_window2, eigen_list$nodes_window2, NULL))
sigma3 <- array(NA, dim = c(eigen_list$num_nodes_window3, eigen_list$num_nodes_window3, nrow(params_ustd)),
                dimnames = list(eigen_list$nodes_window3, eigen_list$nodes_window3, NULL))
for(i in 1:nrow(params)){
  sigma1[,,i] <- sim_cent_cov_1 + diag(rep(params_ustd$sigma[i], eigen_list$num_nodes_window1))
  sigma2[,,i] <- sim_cent_cov_2 + diag(rep(params_ustd$sigma[i], eigen_list$num_nodes_window2))
  sigma3[,,i] <- sim_cent_cov_3 + diag(rep(params_ustd$sigma[i], eigen_list$num_nodes_window3))
}

predictions <- matrix(NA, nrow = nrow(params_ustd), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(predictions)){
  predictions[i,sim$window == 1] <- MASS::mvrnorm(1, mu[i,sim$window == 1], sigma1[,,i])
  predictions[i,sim$window == 2] <- MASS::mvrnorm(1, mu[i,sim$window == 2], sigma2[,,i])
  predictions[i,sim$window == 3] <- MASS::mvrnorm(1, mu[i,sim$window == 3], sigma3[,,i])
}

pred_data <- sim
pred_data$mu_model <- NA
pred_data$prediction <- NA
for(i in 1:nrow(pred_data)){
  pred_data$mu_model[i] <- list(mu[,i])
  pred_data$prediction[i] <- list(predictions[,i])
}
pred_data <- unnest(cols = prediction, data = pred_data)

sim$mu_lwr <- NA ; sim$mu_upr <- NA
sim$predict_lwr <- NA ; sim$predict_upr <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr[i] <- rethinking::HPDI(mu[,i], prob = 0.95)[1]
  sim$mu_upr[i] <- rethinking::HPDI(mu[,i], prob = 0.95)[2]
  sim$predict_lwr[i] <- rethinking::HPDI(predictions[,i], prob = 0.95)[1]
  sim$predict_upr[i] <- rethinking::HPDI(predictions[,i], prob = 0.95)[2]
}

## plot mean of model vs mean of raw data
plot(sim$mu_mean ~ sim$mu, pch = 19, col = rgb(0,0,1,0.2)) ; abline(a = 0, b = 1, lwd = 2)

## plot all predictions with points = mean of raw data
ggplot(sim)+
  geom_ribbon(aes(x = age, ymin = predict_lwr, ymax = predict_upr, fill = as.factor(window)),
              alpha = 0.2)+
  geom_ribbon(aes(x = age, ymin = mu_lwr, ymax = mu_upr, fill = as.factor(window)),
              alpha = 0.4)+
  geom_line(aes(x = age, y = mu_mean, colour = as.factor(window)))+
  geom_point(aes(x = age, y = mu), pch = 19, size = 0.8)+
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

## extract original values from output -- simulated slope value originally used produces an effect on the unstandardised scale. The model works on the standardised scale. Convert predictions to unstandardised scale and then run contrasts to calculate the slope coefficient. ####
## predict from model again, now using age_std + 1 stdev
mu2 <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(mu2)){
  mu2[i,] <- params$beta_age[i] * (sim$age_std+1) + params$intercept[i] + as.numeric(rand_window[i,sim$window]) + as.numeric(rand_node[i, sim$node])
}
sim$mu_mean_plus1yr <- apply(mu2, 2, mean)

predictions2 <- matrix(NA, nrow = nrow(params), ncol = nrow(sim), dimnames = list(NULL, sim$node_window))
for(i in 1:nrow(predictions2)){
  predictions2[i,sim$window == 1] <- MASS::mvrnorm(1, mu2[i,sim$window == 1], sigma1[,,i])
  predictions2[i,sim$window == 2] <- MASS::mvrnorm(1, mu2[i,sim$window == 2], sigma2[,,i])
  predictions2[i,sim$window == 3] <- MASS::mvrnorm(1, mu2[i,sim$window == 3], sigma3[,,i])
}

pred_data2 <- sim
pred_data2$mu_model <- NA
pred_data2$prediction <- NA
for(i in 1:nrow(pred_data2)){
  pred_data2$mu_model[i] <- list(mu2[,i])
  pred_data2$prediction[i] <- list(predictions2[,i])
}
pred_data2 <- unnest(cols = prediction, data = pred_data2)

## revert predictions to unstandardised scale
mu2_ustd <- mu2 * sd(sim$mu) + mean(sim$mu)                   # don't know that I can do this when this wasn't the same raw data
predictions2_ustd <- predictions2 * sd(sim$mu) + mean(sim$mu) # don't know that I can do this when this wasn't the same raw data

sim$mu_mean_ustd_plus1yr <- apply(mu2_ustd, 2, mean)
sim$mu_lwr_ustd_plus1yr <- NA ; sim$mu_upr_ustd_plus1yr <- NA
sim$predict_lwr_ustd_plus1yr <- NA ; sim$predict_upr_ustd_plus1yr <- NA
for( i in 1:nrow(sim) ){
  sim$mu_lwr_ustd_plus1yr[i] <- rethinking::HPDI(mu2_ustd[,i], prob = 0.95)[1]
  sim$mu_upr_ustd_plus1yr[i] <- rethinking::HPDI(mu2_ustd[,i], prob = 0.95)[2]
  sim$predict_lwr_ustd_plus1yr[i] <- rethinking::HPDI(predictions2_ustd[,i], prob = 0.95)[1]
  sim$predict_upr_ustd_plus1yr[i] <- rethinking::HPDI(predictions2_ustd[,i], prob = 0.95)[2]
}

## contrast predictions on standardised scale -- check that this returns the marginal effect presented in the summary
contrast_std <- predictions2 - predictions
head(contrast_std[,1:5])
summary$mean[which(summary$variable == 'beta_age')]
mean(contrast_std)
quantile(contrast_std, prob = c(0.025, 0.975))

## contrast predictions on unstandardised scale -- should return marginal effect originally set?
contrast_ustd <- predictions2_ustd - predictions_ustd
head(contrast_ustd[,1:5])
mean(contrast_ustd)
sim_slope

contrast_ustd <- contrast_std * sd(sim$mu) + mean(sim$mu)


## final plots -- predict from hypothetical data to get straight line and smooth ribbons, rather than jagged ####
pred_fake <- data.frame(age_std = rep(seq(min(sim$age_std), max(sim$age_std), length.out = 5), each = nrow(sim)),
                        window = rep(sim$window,5),
                        node = rep(sim$node,5),
                        node_window = rep(sim$node_window, 5))
mu_fake <- matrix(NA, nrow = nrow(params), ncol = nrow(pred_fake), dimnames = list(NULL, pred_fake$node_window))
for(i in 1:nrow(mu_fake)){
  mu_fake[i,] <- params$beta_age[i] * pred_fake$age_std + params$intercept[i] + as.numeric(rand_window[i,pred_fake$window]) + as.numeric(rand_node[i, pred_fake$node])
}
pred_fake$mu <- apply(mu_fake, 2, mean)
pred_fake_summary <- pred_fake %>% 
  group_by(age_std, window) %>% 
  mutate(mu_mean = mean(mu),
         mu_lwr = quantile(mu, 0.025),
         mu_upr = quantile(mu, 0.975)) %>% 
  select(age_std, window, mu_mean, mu_lwr, mu_upr) %>% 
  distinct()

predictions_fake <- matrix(NA, nrow = nrow(params), ncol = nrow(pred_fake), dimnames = list(NULL, pred_fake$node_window))
for(i in 1:nrow(predictions)){
  for(age in unique(pred_fake$age_std)){
    predictions_fake[i,pred_fake$window == 1 & 
                       pred_fake$age_std == age] <- MASS::mvrnorm(1,
                                                                  mu_fake[i,pred_fake$window == 1 & 
                                                                            pred_fake$age_std == age],
                                                                  sigma1[,,i])
    predictions_fake[i,pred_fake$window == 2 & 
                       pred_fake$age_std == age] <- MASS::mvrnorm(1,
                                                                  mu_fake[i,pred_fake$window == 2 & 
                                                                            pred_fake$age_std == age],
                                                                  sigma2[,,i])
    predictions_fake[i,pred_fake$window == 3 & 
                       pred_fake$age_std == age] <- MASS::mvrnorm(1,
                                                                  mu_fake[i,pred_fake$window == 3 & 
                                                                            pred_fake$age_std == age],
                                                                  sigma3[,,i])
  }
}
pred_fake_summary$pred_mean <- apply(predictions_fake, 2, mean)

pred_fake_summary$age <- pred_fake_summary$age_std * sd(sim$age) + mean(sim$age)
pred_fake_summary$mu_mean_ustd <- pred_fake_summary$mu_mean * sd(sim$mu) + mean(sim$mu)
pred_fake_summary$mu_lwr_ustd <- pred_fake_summary$mu_lwr * sd(sim$mu) + mean(sim$mu)
pred_fake_summary$mu_upr_ustd <- pred_fake_summary$mu_upr * sd(sim$mu) + mean(sim$mu)

ggplot()+
  geom_ribbon(data = pred_fake_summary, aes(x = age, ymin = pred_lwr_ustd, ymax = pred_upr_ustd, fill = as.factor(window)),
              alpha = 0.2)+
  geom_ribbon(data = sim, aes(x = age, ymin = mu_lwr, ymax = mu_upr, fill = as.factor(window)),
              alpha = 0.4)+
  geom_line(data = pred_fake_summary, aes(x = age, y = pred_mean_ustd, colour = as.factor(window)))+
  geom_point(data = sim, aes(x = age, y = mu), pch = 19, size = 0.8)+
  scale_colour_viridis_d(begin = 0, end = 0.7)+
  scale_fill_viridis_d(begin = 0, end = 0.7)+
  facet_wrap(. ~ as.factor(window))+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(colour = 'time window', fill = 'time window',
       y = 'eigenvector centrality', x = 'age (years)')

# ## run model -- age as an ordered categorical variable with 1 window ####
## remove later windows
# sim <- sim %>% 
#   filter(window == 1) %>% 
#   dplyr::select(-window, -node_window)
# 
# ## categorise age
# sim$age_cat <- ifelse(sim$age <= 15, 1,
#                       ifelse(sim$age <= 20, 2,
#                              ifelse(sim$age <= 25, 3,
#                                     ifelse(sim$age <= 40, 4, 5))))
# 
# ## create data
# n_age_cat <- length(unique(sim$age_cat))
# eigen_list <- list(#num_data = n_data,
#   num_nodes = n_nodes,
#   num_age_cat = n_age_cat,
#   length_dirichlet = n_age_cat + 1,
#   #num_windows = length(unique(sim$window)),
#   centrality_mu = sim_cent_mu,
#   centrality_cov = sim_cent_cov,
#   node_age = sim$age_cat,
#   #window = sim$window,
#   #node = sim$node,
#   prior_age = rep(1, n_age_cat))
# 
# ## check inputs
# plot(sim_cent_mu ~ sim$age_cat)
# 
# ## prior predictive check
# n <- 100
# beta_age <- rnorm(n, 0, 1)
# intercept  <- rnorm(n, 0, 0.8)
# age_dirichlet <- rdirichlet(n, c(1,1,1,1,1))
# plot(NULL, las = 1, xlab = 'age (standardised)', ylab = 'eigenvector (standardised)',
#      ylim = c(min(sim_cent_mu)-1, max(sim_cent_mu)+1), xlim = c(min(sim$age_cat), max(sim$age_cat)))
# abline(h = min(sim_cent_mu), lty = 2) ; abline(h = max(sim_cent_mu), lty = 2)
# x <- min(sim$age_cat):max(sim$age_cat)
# for(i in 1:n){
#   y <- rep(NA, length(x))
#   for(j in 1:length(x)){
#     y[j] <- intercept[i] + beta_age[i]*sum(age_dirichlet[i,][1:x[j]])
#   }
#   lines(x = x, y = y, col = rgb(0,0,1,0.4))
# }
# rm(n, beta_age, intercept, age_dirichlet, sigma, x, y, df_plot, df_wide) ; gc()
# 
# ## load model
# nodal_regression <- cmdstan_model('models/eigen_regression_motnp.stan')
# 
# ## run model
# n_chains <- 4
# n_samples <- 1000
# #fit_sim <- sampling(nodal_regression, data = eigen_list, chains = n_chains, cores = n_chains)
# fit_sim <- nodal_regression$sample(data = eigen_list,
#                                    chains = n_chains, parallel_chains = n_chains,
#                                    iter_warmup = n_samples, iter_sampling = n_samples)
# 
# ## check outputs ####
# ## view summary
# fit_sim$summary()
# 
# ## extract posterior
# #params <- rstan::extract(fit_sim)
# params <- fit_sim$draws(format = 'draws_df')
# #rand_window <- params %>% 
# #  dplyr::select(`rand_window[1]`,`rand_window[2]`)
# #rand_node <- params[,16:(n_nodes+2)]
# delta <- params %>% 
#   select(`delta[1]`,`delta[2]`,`delta[3]`,`delta[4]`,`delta[5]`)
# delta_j <- params %>% 
#   select(`delta_j[1]`,`delta_j[2]`,`delta_j[3]`,`delta_j[4]`,`delta_j[5]`,`delta_j[6]`)
# 
# ## traceplot linear effect size
# #traceplot(fit_sim, pars = c('intercept','beta_age','sigma','predictor[1]','predictor[50]','predictor[100]'))
# parameters_to_check <- c('intercept','beta_age','sigma','predictor[1]','predictor[20]','predictor[50]')#, 'rand_window[1]','rand_window[2]','rand_node[1]','rand_node[50]','rand_node[100]')
# params %>% 
#   select(all_of(parameters_to_check)) %>%
#   pivot_longer(cols = everything(), names_to = 'parameter', values_to = 'draw') %>% 
#   mutate(chain_position = rep(rep(1:n_samples, each = length(parameters_to_check)), n_chains),
#          chain = rep(1:n_chains, each = length(parameters_to_check)*n_samples)) %>% 
#   #filter(chain == 4) %>% 
#   ggplot(aes(x = chain_position, y = draw, colour = as.factor(chain)))+
#   geom_line()+
#   facet_wrap(. ~ parameter, scales = 'free_y')+
#   theme_bw()+
#   theme(legend.position = 'none')
# 
# ## posterior predictive check
# plot(density(sim_dat_std[1, ]), las = 1, ylim = c(0,0.4),
#      main = "Posterior predictive check (standardised centrality):\nblack = data, blue = predicted",
#      col=rgb(0, 0, 0, 0.25))
# for (i in 1:100) {
#   j <- sample(1:length(params$beta_age), 1)
#   lines(density(sim_dat_std[j, ]), col=rgb(0, 0, 0, 0.25))
#   mu <- rep(NA, length(eigen_list$node_age))
#   for(k in 1:length(eigen_list$node_age)){
#     mu[k] <- params$intercept[j] + params$beta_age[j]*sum(delta_j[j,(1:eigen_list$node_age[k])]) #+ as.numeric(rand_window[j,eigen_list$window[k]]) + as.numeric(rand_node[j,eigen_list$nodes[k]])
#   }
#   sigma <- sim_cent_cov + diag(rep(params$sigma[j], n_data))
#   lines(density(MASS::mvrnorm(1, mu, sigma)), col=rgb(0, 0, 1, 0.25))
# }
# 
# ## predict from model
# new_data <- data.frame(node = 1:length(unique(sim$age_cat)),
#                        age_cat = rep(sort(unique(sim$age_cat))))
# predict <- matrix(data = NA, nrow = length(params$beta_age), ncol = nrow(new_data),
#                   dimnames = list(1:length(params$beta_age),
#                                   1:length(unique(sim$age_cat))))
# for(i in 1:nrow(predict)){
#   for(j in 1:ncol(predict)){
#     predict[i,j] <- MASS::mvnorm(n = 1,
#                                  mu = params$intercept[i] + params$beta_age[i] * sum(delta_j[i,(1:j)]),
#                                  Sigma = params$sigma[i])
#     #predict[i,j] <- params$intercept[i] + params$beta_age[i] * sum(delta_j[i,(1:j)])
#   }
# }
# new_data$mean_predict <- apply(predict, 2, mean)
# new_data$lwr_predict <- apply(predict, 2, quantile, probs = 0.025)
# new_data$upr_predict <- apply(predict, 2, quantile, probs = 0.975)
# 
# ## compare to standardised raw data
# plot(new_data$mean_predict ~ new_data$age_cat, las = 1, type = 'l', ylim = c(-2,2))        # plot against age
# points(sim$mu_std ~ sim$age_cat)                                                           # add raw points
# polygon(y = c(new_data$lwr_predict, rev(new_data$upr_predict)), x = c(new_data$age_cat, rev(new_data$age_cat)),
#         col = rgb(1,1,0,0.5), border = NA)
# 
# ## unstandardise predictions
# sim$mu
# (new_data$predict_ustd <- new_data$mean_predict * sd(sim$mu) + mean(sim$mu))
# sim <- sim %>% 
#   left_join(new_data, by = 'age_cat')
# plot(sim$predict_ustd ~ sim$mu)
# abline(a = 0, b = 1)
# 
# ## compare to unstandardised raw data
# sim$lwr_ustd <- sim$lwr_predict * sd(sim$mu) + mean(sim$mu)
# sim$upr_ustd <- sim$upr_predict * sd(sim$mu) + mean(sim$mu)
# plot(sim$predict_ustd ~ sim$age_cat, las = 1, type = 'l')        # plot against age
# points(sim$mu ~ sim$age_cat)                                     # add raw points
# polygon(y = c(sim$lwr_ustd, rev(sim$upr_ustd)), x = c(sim$age_cat, rev(sim$age_cat)),
#         col = rgb(1,1,0,0.5), border = NA)
# 
# ## convert to invlogit scale
# sim$mean_predict_invlogit <- invlogit(sim$predict_ustd)
# sim$lwr_invlogit <- invlogit(sim$lwr_ustd)
# sim$upr_invlogit <- invlogit(sim$upr_ustd)
# sim$mu_invlogit <- invlogit(sim$mu)
# 
# ## compare to invlogit raw data
# plot(sim$mu_invlogit ~ sim$age_cat, ylim = c(0,1))
# lines(sim$mean_predict_invlogit ~ sim$age_cat)
# polygon(y = c(sim$lwr_invlogit, rev(sim$upr_invlogit)), x = c(sim$age_cat, rev(sim$age_cat)),
#         col = rgb(1,1,0,0.5), border = NA)
# 
# ## extract slope estimates -- none of these match the original input, but the predictions work??!
# # original (true) slope value
# sim_slope
# # predicted (modelled) slope value = mean of differences bewteen categories, divided by the number of years each category represents
# mean( c( (new_data$predict_ustd[2] - new_data$predict_ustd[1]) / (16-10), 
#          (new_data$predict_ustd[3] - new_data$predict_ustd[2]) / (21-16), 
#          (new_data$predict_ustd[4] - new_data$predict_ustd[3]) / (26-21),
#          (new_data$predict_ustd[5] - new_data$predict_ustd[4]) / (40-26) )
# )
# 
# ## extract model fit
# summary <- fit_sim$summary()
# par(mfrow = c(3,1))
# hist(summary$rhat, breaks = 50)
# hist(summary$ess_bulk, breaks = 50)
# hist(summary$ess_tail, breaks = 50)
# par(mfrow = c(1,1))
# 

##### ignore all this -- it's old stuff that I want to be certain I don't need before I delete it! ####
# ages <- seq(min(sim$age_std), max(sim$age_std), length.out = 10)
# pred_data <- data.frame(age_std = rep(ages, each = eigen_list$num_nodes_window1),
#                         node = rep(1:eigen_list$num_nodes_window1, length(ages)),
#                         window = 1,
#                         prediction = NA)
# predictions1 <- pred_data
# for(i in 1:length(ages)){
#   mu <- params$beta_age[1]*ages[i] + params$intercept[1] + as.numeric(rand_window[1,1]) + as.numeric(rand_node[1, predictions1$node[predictions$age_std == ages[i]]])
#   sigma <- sim_cent_cov_1 + diag(rep(params$sigma[1], eigen_list$num_nodes_window1))
#   predictions1$prediction[predictions1$age_std == ages[i]] <- MASS::mvrnorm(1, mu, sigma)
# }
# for(age in 1:length(ages)){
#   for(draw in 2:nrow(params)) {
#     mu <- params$beta_age[draw]*ages[age] + params$intercept[draw] + as.numeric(rand_window[draw,1]) + as.numeric(rand_node[draw, pred_data$node[pred_data$age_std == ages[age]]])
#     sigma <- sim_cent_cov_1 + diag(rep(params$sigma[draw], eigen_list$num_nodes_window1))
#     pred_data$prediction[pred_data$age_std == ages[age]] <- MASS::mvrnorm(1, mu, sigma)
#     predictions1 <- rbind(predictions1, pred_data)
#     if(draw %% 100 == 0) { print(draw) }
#   }
# }
# save.image('step4_nodalregression/simulation_anp.RData')
# 
# pred_data <- data.frame(age_std = rep(ages, each = eigen_list$num_nodes_window2),
#                         node = rep(1:eigen_list$num_nodes_window2, length(ages)),
#                         window = 2,
#                         prediction = NA)
# predictions2 <- pred_data
# for(age in 1:length(ages)){
#   mu <- params$beta_age[1]*ages[age] + params$intercept[1] + as.numeric(rand_window[1,2]) + as.numeric(rand_node[1, predictions2$node[predictions2$age_std == ages[age]]])
#   sigma <- sim_cent_cov_2 + diag(rep(params$sigma[1], eigen_list$num_nodes_window2))
#   predictions2$prediction[predictions2$age_std == ages[age]] <- MASS::mvrnorm(1, mu, sigma)
# }
# for(age in 1:length(ages)){
#   for(draw in 2:nrow(params)) {
#     mu <- params$beta_age[draw]*ages[age] + params$intercept[draw] + as.numeric(rand_window[draw,2]) + as.numeric(rand_node[draw, pred_data$node[pred_data$age_std == ages[age]]])
#     sigma <- sim_cent_cov_2 + diag(rep(params$sigma[draw], eigen_list$num_nodes_window2))
#     pred_data$prediction[pred_data$age_std == ages[age]] <- MASS::mvrnorm(1, mu, sigma)
#     predictions2 <- rbind(predictions2, pred_data)
#     if(draw %% 100 == 0) { print(draw) }
#   }
# }
# predictions <- rbind(predictions1, predictions2)
# save.image('step4_nodalregression/simulation_anp.RData')
# rm(mu, sigma, pred_data, age, draw, i, j ,k) ; gc()

## plot predictions
pred_data <- predictions %>% 
  group_by(age_std, node, window) %>% 
  mutate(mean = mean(prediction)) %>% 
  ungroup() %>% 
  select(age_std, node, window, mean) %>% 
  distinct() %>% 
  left_join(distinct(sim[,c('age','age_std')]), by = 'age_std')
ggplot()+
  geom_point(data = sim,
             mapping = aes(x = age_std, y = mu_std,
                           colour = as.factor(node)))+
  geom_line(data = pred_data,
            mapping = aes(x = age_std, y = mean))+
  facet_wrap(. ~ window)+
  theme_classic()+
  theme(legend.position = 'none')
